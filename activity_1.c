#include <stdio.h> 
#include <stdlib.h>
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

int getIdByCoord(MPI_Comm comm, int y, int x);
float float_rand(float min, float max) ;
void getCurrentDatetime(int* currentDatetime);
void getReadings(float* readings);
double deg2rad(double);
double rad2deg(double);
double distance(double lat1, double lon1, double lat2, double lon2, char unit);

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

  /* Magnitude exceeds predefined threshold 
  */

  int left_rank, right_rank, top_rank, bottom_rank;
  MPI_Cart_shift(comm, SHIFT_ROW, DISP, &top_rank, &bottom_rank);
  MPI_Cart_shift(comm, SHIFT_COL, DISP, &left_rank, &right_rank);

  MPI_Request send_request[4];
  MPI_Request receive_request[4];
  MPI_Status send_status[4];
  MPI_Status receive_status[4];

  MPI_Isend(1, 10, MPI_FLOAT, top_rank, 0, comm, &send_request[0]);

  // Share all readings to neighbouring nodes
  MPI_Isend(readings, 10, MPI_FLOAT, top_rank, 0, comm, &send_request[0]);
  MPI_Isend(readings, 10, MPI_FLOAT, bottom_rank, 0, comm, &send_request[1]);
  MPI_Isend(readings, 10, MPI_FLOAT, left_rank, 0, comm, &send_request[2]);
  MPI_Isend(readings, 10, MPI_FLOAT, right_rank, 0, comm, &send_request[3]);

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

  float magnitude = readings[8];
  if (magnitude > MAGNITUDE_UPPER_THRESHOLD) {
    
  }

  // if(rank == 2) {
  //   // printf("Rank %d's Top sensor magnitude: %f \n", rank, recvBufT[8]);
  //   // printf("Rank %d's Bottom sensor magnitude: %f \n", rank, recvBufB[8]);
  //   // printf("Rank %d's Left sensor magnitude: %f \n", rank, recvBufL[8]);
  //   // printf("Rank %d's Right sensor magnitude: %f \n", rank, recvBufR[8]);

  //   for(int i=0; i<10; i++) {
  //     printf("Rank %d => recvBufT[%d]: %f \n", rank, i, recvBufT[i]);
  //   }
  //   for(int i=0; i<10; i++) {
  //     printf("Rank %d => recvBufB[%d]: %f \n", rank, i, recvBufB[i]);
  //   }
  //   for(int i=0; i<10; i++) {
  //     printf("Rank %d => recvBufL[%d]: %f \n", rank, i, recvBufL[i]);
  //   }
  //   for(int i=0; i<10; i++) {
  //     printf("Rank %d => recvBufR[%d]: %f \n", rank, i, recvBufR[i]);
  //   }
  // }


  // int id;
  // int x, y;
  // x = coord[0]; y = coord[1];

  // // Checks if magnitude over threshold
  // float latitude = readings[0];
  // float longitude = readings[1];
  // float magnitude = readings[2];



  // if(rank == 0)
  // {
  //   // getIdByCoord(comm, 1, 1);
  //   double d = distance(40, 100, -45, -45, 'M');
  //   printf("distance: %f", d);
  // }

	MPI_Finalize();
  return 0;
}

// Returns -2 if neighbour node out of bounds
int getIdByCoord(MPI_Comm comm, int y, int x) {
  int id;
  int coord[2] = {y, x};
  
  MPI_Cart_rank(comm, coord, &id);
  printf("The processor at position (%d, %d) has rank %d\n", coord[0], coord[1], id);
  fflush(stdout);

  return id;
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

/* generate a random floating point number from min to max */
float float_rand(float min, float max) 
{
    float range = (max - min); 
    float div = RAND_MAX / range;
    return min + (rand() / div);
}

double distance(double lat1, double lon1, double lat2, double lon2, char unit) {
  double theta, dist;
  if ((lat1 == lat2) && (lon1 == lon2)) {
    return 0;
  }
  else {
    theta = lon1 - lon2;
    dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
    dist = acos(dist);
    dist = rad2deg(dist);
    dist = dist * 60 * 1.1515;
    switch(unit) {
      case 'M':
        break;
      case 'K':
        dist = dist * 1.609344;
        break;
      case 'N':
        dist = dist * 0.8684;
        break;
    }
    return (dist);
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