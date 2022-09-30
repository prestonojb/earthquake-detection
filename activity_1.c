#include <stdio.h> 
#include <stdlib.h>
#include <stdbool.h>
#include "mpi.h"
#include <math.h>
#include <time.h>

#define pi 3.14159265358979323846
#define LATITUDE_LOWER_BOUND -90
#define LATITUDE_UPPER_BOUND 90
#define LONGITUDE_LOWER_BOUND -180
#define LONGITUDE_UPPER_BOUND 180
#define MAGNITUDE_LOWER_BOUND 0
#define MAGNITUDE_UPPER_BOUND 9
#define DEPTH_UPPER_BOUND 700
#define MAGNITUDE_UPPER_THRESHOLD 0

#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1

#define DIFF_IN_DISTANCE_THRESHOLD_IN_KM 100
#define DIFF_IN_MAGNITUDE_THRESHOLD 100

float float_rand(float min, float max) ;
bool areMatchingReadings(float* readingA, float* readingB);
bool areMatchingLocations(float latA, float longA, float latB, float longB);
bool areMatchingMagnitudes(float magnitudeA, float magnitudeB);
void getCurrentDatetime(int* currentDatetime);
void getReadings(float* readings);
void printReadings(int rank, float* readings);
double deg2rad(double);
double rad2deg(double);
float distance(float lat1, float lon1, float lat2, float lon2);

int main(int argc, char* argv[]) {
  // Initialise MPI
  int rank, total_nodes;
	MPI_Status status;
  MPI_Comm comm;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &total_nodes);

  // Seed rand() function with current timestamp and rank
  srand(time(NULL) + rank*100);

  // Get dimensions of topology from command line arguments
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);

  // Initilise MPI virtual topology (Cartesian)
  int dim[2] = {m, n};
  int wrap_around[2] = {0,0};
  int reorder = 1;
  int coord[2];

  int ierr = 0;
  ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, dim, wrap_around, reorder, &comm);
  if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);

  MPI_Cart_coords(comm, rank, 2, coord); // Save current coordinates

  // if (rank == 0) {
  //   printf("Topology dimensions: %dx%d \n", m, n);
  // }
  // printf("Current rank/Total nodes: %d/%d, Coordinate: (%d, %d) \n", rank, total_nodes, coord[0], coord[1]);

  
  /* Generate reading
     =======
      readings = [YYYY, MM, DD, HH, MM, SS, latitude, longitude, magnitude, depth]
  */
  float* readings = (float*)malloc(10 * sizeof(float));
  getReadings(readings);
  // if (rank == 0) {
  //   for (int i=0; i<10; i++) {
  //     printf("Rank %d => readings[%d]: %.2f \n", rank, i, readings[i]);
  //   }
  // }

  int left_rank, right_rank, top_rank, bottom_rank;
  MPI_Cart_shift(comm, SHIFT_ROW, DISP, &top_rank, &bottom_rank);
  MPI_Cart_shift(comm, SHIFT_COL, DISP, &left_rank, &right_rank);

  MPI_Request send_request[4];
  MPI_Request receive_request[4];
  MPI_Status send_status[4];
  MPI_Status receive_status[4];

  float magnitude = readings[8];
  int compare_readings = 0;
  if (magnitude > MAGNITUDE_UPPER_THRESHOLD) {
    compare_readings = 1;
  }
  // if (rank == 0) {
  //   compare_readings = 1;
  // }

  // Request reading from adjacent nodes
  MPI_Isend(&compare_readings, 1, MPI_INT, top_rank, 0, comm, &send_request[0]);
  MPI_Isend(&compare_readings, 1, MPI_INT, bottom_rank, 0, comm, &send_request[1]);
  MPI_Isend(&compare_readings, 1, MPI_INT, left_rank, 0, comm, &send_request[2]);
  MPI_Isend(&compare_readings, 1, MPI_INT, right_rank, 0, comm, &send_request[3]);
  
  int compare_readings_T = -1, compare_readings_B = -1, compare_readings_L = -1, compare_readings_R = -1;
  MPI_Irecv(&compare_readings_T, 1, MPI_INT, top_rank, 0, comm, &receive_request[0]);
  MPI_Irecv(&compare_readings_B, 1, MPI_INT, bottom_rank, 0, comm, &receive_request[1]);
  MPI_Irecv(&compare_readings_L, 1, MPI_INT, left_rank, 0, comm, &receive_request[2]);
  MPI_Irecv(&compare_readings_R, 1, MPI_INT, right_rank, 0, comm, &receive_request[3]);

  MPI_Waitall(4, send_request, send_status);
  MPI_Waitall(4, receive_request, receive_status);

  // if (rank == 0) {
  //   printf("Rank %d => Compare readings T: %d \n", rank, compare_readings_T);
  //   printf("Rank %d => Compare readings B: %d \n", rank, compare_readings_B);
  //   printf("Rank %d => Compare readings L: %d \n", rank, compare_readings_L);
  //   printf("Rank %d => Compare readings R: %d \n", rank, compare_readings_R);
  // }

  // Send readings to adjacent nodes (if requested)
  if (compare_readings_T == 1) {
    MPI_Isend(readings, 10, MPI_FLOAT, top_rank, 0, comm, &send_request[0]);
  }
  if (compare_readings_B == 1) {
    MPI_Isend(readings, 10, MPI_FLOAT, bottom_rank, 0, comm, &send_request[1]);
  }
  if (compare_readings_L == 1) {
    MPI_Isend(readings, 10, MPI_FLOAT, left_rank, 0, comm, &send_request[2]);
  }
  if (compare_readings_R == 1) {
    MPI_Isend(readings, 10, MPI_FLOAT, right_rank, 0, comm, &send_request[3]);
  }

  // if (rank == 2) {
  //   printf("Rank %d readings: ", rank);
  //   printReadings(rank, readings);
  //   printf("\n");
  // }

  // Receive requested readings
  if (compare_readings == 1) {
    float* readingsT = calloc(10, sizeof(float));
    float* readingsB = calloc(10, sizeof(float));
    float* readingsL = calloc(10, sizeof(float));
    float* readingsR = calloc(10, sizeof(float));
    
    MPI_Irecv(readingsT, 10, MPI_FLOAT, top_rank, 0, comm, &receive_request[0]);
    MPI_Irecv(readingsB, 10, MPI_FLOAT, bottom_rank, 0, comm, &receive_request[1]);
    MPI_Irecv(readingsL, 10, MPI_FLOAT, left_rank, 0, comm, &receive_request[2]);
    MPI_Irecv(readingsR, 10, MPI_FLOAT, right_rank, 0, comm, &receive_request[3]);
    
    MPI_Waitall(4, send_request, send_status);
    MPI_Waitall(4, receive_request, receive_status);

    // if (rank == 0) {
      // printf("Own reading: ");
      // printReadings(rank, readings);
      // printf("\n");

      // printf("readingsT: ");
      // printReadings(rank, readingsT);
      // printf("\n");
      
      // printf("readingsB: ");
      // printReadings(rank, readingsB);
      // printf("\n");
      
      // printf("readingsL: ");
      // printReadings(rank, readingsL);
      // printf("\n");
      
      // printf("readingsR: ");
      // printReadings(rank, readingsR);
      // printf("\n");
    // }

    int no_of_matches = 0;
    if(compare_readings_T == 1&& areMatchingReadings(readings, readingsT)) no_of_matches++;
    if(compare_readings_B == 1 && areMatchingReadings(readings, readingsB)) no_of_matches++;
    if(compare_readings_L == 1 && areMatchingReadings(readings, readingsL)) no_of_matches++;
    if(compare_readings_R == 1 && areMatchingReadings(readings, readingsR)) no_of_matches++;

    if(no_of_matches >= 2) {
      // Send report to base station #TODO
      printReadings(rank, readings);
    }
  }

	MPI_Finalize();
  return 0;
}

void printReadings(int rank, float* readings) {
  for(int i=0; i<10; i++) {
    printf("Rank %d => readings[%d]: %f \n", rank, i, readings[i]);
  }
}

// Readings = [YYYY, MM, DD, HH, MM, SS, latitude, longitude, magnitude, depth]
void getReadings(float* readings) {
  // Generate datetime (reading #1-6)
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);

  float_rand(0, 0); // Workaround for same first random number
  float latitude = float_rand(LATITUDE_LOWER_BOUND, LATITUDE_UPPER_BOUND);
  float longitude = float_rand(LONGITUDE_LOWER_BOUND, LONGITUDE_UPPER_BOUND);
  float magnitude = float_rand(MAGNITUDE_LOWER_BOUND, MAGNITUDE_UPPER_BOUND);
  float depth = float_rand(0, DEPTH_UPPER_BOUND);

  // Current datetime
  readings[0] = tm.tm_year + 1900;
  readings[1] = tm.tm_mon + 1;
  readings[2] = tm.tm_mday;
  readings[3] = tm.tm_hour;
  readings[4] = tm.tm_min;
  readings[5] = tm.tm_sec;

  // Sensor readings
  readings[6] = latitude;
  readings[7] = longitude;
  readings[8] = magnitude;
  readings[9] = depth;
}

bool areMatchingReadings(float* readingA, float* readingB) {
  float latA = readingA[6];
  float longA = readingA[7];
  float latB = readingA[6];
  float longB = readingA[7];

  float magnitudeA = readingA[8];
  float magnitudeB = readingB[8];

  return areMatchingLocations(latA, longA, latB, longB) && areMatchingMagnitudes(magnitudeA, magnitudeB);
}

bool areMatchingLocations(float latA, float longA, float latB, float longB) {
  return distance(latA, longA, latB, longB) < DIFF_IN_DISTANCE_THRESHOLD_IN_KM;
}

bool areMatchingMagnitudes(float magnitudeA, float magnitudeB) {
  return fabs((double) magnitudeA - magnitudeB) < DIFF_IN_MAGNITUDE_THRESHOLD;
}

/* generate a random floating point number from min to max */
float float_rand(float min, float max) 
{
    float range = (max - min); 
    float div = RAND_MAX / range;
    return min + (rand() / div);
}

float distance(float lat1, float lon1, float lat2, float lon2) {
  float theta, dist;
  if ((lat1 == lat2) && (lon1 == lon2)) {
    return 0;
  }
  else {
    theta = lon1 - lon2;
    dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
    dist = acos(dist);
    dist = rad2deg(dist);
    dist = dist * 60 * 1.1515;
    return dist * 1.609344;
  }
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts decimal degrees to radians             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double deg2rad(double deg) {
  return (deg * pi / 180);
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts radians to decimal degrees             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double rad2deg(double rad) {
  return (rad * 180 / pi);
}