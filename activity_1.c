#include <stdio.h> 
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include "mpi.h"
#include <math.h>
#include <time.h>
#include "helper.c"

#define pi 3.14159265358979323846
#define LATITUDE_LOWER_BOUND -90
#define LATITUDE_UPPER_BOUND 90
#define LONGITUDE_LOWER_BOUND -180
#define LONGITUDE_UPPER_BOUND 180
#define MAGNITUDE_LOWER_BOUND 0
#define MAGNITUDE_UPPER_BOUND 9
#define DEPTH_UPPER_BOUND 700
#define MAGNITUDE_UPPER_THRESHOLD 2.5

#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1

#define DIFF_IN_DISTANCE_THRESHOLD_IN_KM 50
#define DIFF_IN_MAGNITUDE_THRESHOLD 0.5

void generate(struct Sensor* reading);
void printReading(struct Sensor* reading);

bool areMatchingReadings(struct Sensor* readingA, struct Sensor* readingB);
bool areMatchingLocations(float latA, float longA, float latB, float longB);
bool areMatchingMagnitudes(float magnitudeA, float magnitudeB);

float float_rand(float min, float max);
float distance(float lat1, float lon1, float lat2, float lon2);

double deg2rad(double);
double rad2deg(double);

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

  if (rank == 0) {
    printf("Topology dimensions: %dx%d \n", m, n);
  }
  printf("Current rank/Total nodes: %d/%d, Coordinate: (%d, %d) \n", rank, total_nodes, coord[0], coord[1]);
  
  struct Sensor newReading;
  generate(&newReading);
  printf("Rank %d => ", rank);
  printReading(&newReading);

  int left_rank, right_rank, top_rank, bottom_rank;
  MPI_Cart_shift(comm, SHIFT_ROW, DISP, &top_rank, &bottom_rank);
  MPI_Cart_shift(comm, SHIFT_COL, DISP, &left_rank, &right_rank);

  MPI_Request send_request[4];
  MPI_Request receive_request[4];
  MPI_Status send_status[4];
  MPI_Status receive_status[4];

  float magnitude = newReading.mag;
  int compare_readings = 0;
  if (magnitude > MAGNITUDE_UPPER_THRESHOLD) {
    compare_readings = 1;
  }

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

  const int readingSize = 10;
  int blocklengths[10] = {1,1,1,1,1,1,1,1,1,1};
  MPI_Datatype MPI_READING;
  MPI_Aint offsets[10];

  offsets[0] = offsetof(struct Sensor, year);
  offsets[1] = offsetof(struct Sensor, month);
  offsets[2] = offsetof(struct Sensor, day);
  offsets[3] = offsetof(struct Sensor, hour);
  offsets[4] = offsetof(struct Sensor, minute);
  offsets[5] = offsetof(struct Sensor, second);
  offsets[6] = offsetof(struct Sensor, lat);
  offsets[7] = offsetof(struct Sensor, lon);
  offsets[8] = offsetof(struct Sensor, mag);
  offsets[9] = offsetof(struct Sensor, depth);

  MPI_Datatype types[10] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};

  MPI_Type_create_struct(readingSize, blocklengths, offsets, types, &MPI_READING);
  MPI_Type_commit(&MPI_READING);

  // Send (blocking) readings to adjacent nodes (if requested)
  if (compare_readings_T == 1) {
    MPI_Send(&newReading, 1, MPI_READING, top_rank, 0, comm);
  }
  if (compare_readings_B == 1) {
    MPI_Send(&newReading, 1, MPI_READING, bottom_rank, 0, comm);
  }
  if (compare_readings_L == 1) {
    MPI_Send(&newReading, 1, MPI_READING, left_rank, 0, comm);
  }
  if (compare_readings_R == 1) {
    MPI_Send(&newReading, 1, MPI_READING, right_rank, 0, comm);
  }

  // Receive requested readings
  if (compare_readings == 1) {
    struct Sensor readingsT, readingsB, readingsL, readingsR;
    MPI_Status status;
    
    if(top_rank >= 0) {
      MPI_Recv(&readingsT, 1, MPI_READING, top_rank, 0, comm, &status);
    }
    if(bottom_rank >= 0) {
      MPI_Recv(&readingsB, 1, MPI_READING, bottom_rank, 0, comm, &status);
    }
    if(left_rank >= 0) {
      MPI_Recv(&readingsL, 1, MPI_READING, left_rank, 0, comm, &status);
    }
    if(left_rank >= 0) {
      MPI_Recv(&readingsR, 1, MPI_READING, right_rank, 0, comm, &status);
    }
    
    int no_of_matches = 0;
    if(compare_readings_T == 1 && areMatchingReadings(&newReading, &readingsT)) no_of_matches++;
    if(compare_readings_B == 1 && areMatchingReadings(&newReading, &readingsB)) no_of_matches++;
    if(compare_readings_L == 1 && areMatchingReadings(&newReading, &readingsL)) no_of_matches++;
    if(compare_readings_R == 1 && areMatchingReadings(&newReading, &readingsR)) no_of_matches++;

    if(no_of_matches >= 2) {
      // Send report to base station #TODO
      printf("Send report to base station! \n");
    }
  }

  MPI_Type_free(&MPI_READING);
	MPI_Finalize();
  return 0;
}

void generate(struct Sensor* reading)
{
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  reading->year = tm.tm_year + 1900;
  reading->month = tm.tm_mon + 1;
  reading->day = tm.tm_mday;
  reading->hour = tm.tm_hour;
  reading->minute = tm.tm_min;
  reading->second = tm.tm_sec;

  float_rand(0, 0); // Workaround for same first random number
  reading->lat = float_rand(LATITUDE_LOWER_BOUND, LATITUDE_UPPER_BOUND);
  reading->lon = float_rand(LONGITUDE_LOWER_BOUND, LONGITUDE_UPPER_BOUND);
  reading->mag = float_rand(MAGNITUDE_LOWER_BOUND, MAGNITUDE_UPPER_BOUND);
  reading->depth = float_rand(0, DEPTH_UPPER_BOUND);
}

void printReading(struct Sensor* reading)
{
  printf("%d\t| %d\t| %d\t| %d\t| %d\t| %d\t| %.2f\t| %.2f\t| %.2f\t| %.2f\n",
          reading->year, reading->month, reading->day, reading->hour,
          reading->minute, reading->second, reading->lat, reading->lon,
          reading->mag, reading->depth);
}

bool areMatchingReadings(struct Sensor* readingA, struct Sensor* readingB) {
  float latA = readingA->lat;
  float longA = readingA->lon;
  float latB = readingA->lat;
  float longB = readingA->lon;

  float magnitudeA = readingA->mag;
  float magnitudeB = readingB->mag;

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

double deg2rad(double deg) {
  return (deg * pi / 180);
}

double rad2deg(double rad) {
  return (rad * 180 / pi);
}