#include <stdio.h> 
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include "mpi.h"
#include <math.h>
#include <time.h>
#include "helper.c"
#include <pthread.h>

#define pi 3.14159265358979323846
#define LATITUDE_LOWER_BOUND -90
#define LATITUDE_UPPER_BOUND 90
#define LONGITUDE_LOWER_BOUND -180
#define LONGITUDE_UPPER_BOUND 180
#define MAGNITUDE_LOWER_BOUND 0
#define MAGNITUDE_UPPER_BOUND 9
#define DEPTH_UPPER_BOUND 700
#define DEFAULT_MAGNITUDE_UPPER_THRESHOLD 0

#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1

#define DEFAULT_DIFF_IN_DISTANCE_THRESHOLD_IN_KM 5000
#define DEFAULT_DIFF_IN_MAGNITUDE_THRESHOLD 10

void generate(struct Sensor* reading);
void printReading(struct Sensor* reading);

void* AdjNodesCommFunc(void* pArgs);

bool areMatchingReadings(struct Sensor* readingA, struct Sensor* readingB);
bool areMatchingLocations(float latA, float longA, float latB, float longB);
bool areMatchingMagnitudes(float magnitudeA, float magnitudeB);

float float_rand(float min, float max);
float distance(float lat1, float lon1, float lat2, float lon2);

double deg2rad(double);
double rad2deg(double);

struct adj_nodes_arg_struct {
  struct Sensor* pReadingsT;
  struct Sensor* pReadingsB;
  struct Sensor* pReadingsL;
  struct Sensor* pReadingsR;
};

float MAGNITUDE_UPPER_THRESHOLD = DEFAULT_MAGNITUDE_UPPER_THRESHOLD;
float DIFF_IN_DISTANCE_THRESHOLD_IN_KM = DEFAULT_DIFF_IN_DISTANCE_THRESHOLD_IN_KM;
float DIFF_IN_MAGNITUDE_THRESHOLD = DEFAULT_DIFF_IN_MAGNITUDE_THRESHOLD;

int left_rank, right_rank, top_rank, bottom_rank;
int compare_readings = 0;
struct Sensor newReading;
int compare_readings_T = -1, compare_readings_B = -1, compare_readings_L = -1, compare_readings_R = -1;

// Initialise MPI variables
MPI_Status status;
MPI_Comm comm;

int main(int argc, char* argv[]) {
  // Initialise variables
  int rank, total_nodes, cart_rank;
  int provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided );

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &total_nodes);

  // Seed rand() function with current timestamp and rank
  srand(time(NULL) + rank*100);

  int m, n;
  int ndims = 2;
  int dims[ndims], coord[ndims], wrap_around[ndims];
  int reorder = 1;

  wrap_around[0] = wrap_around[1] = 0;

  // Get dimensions of topology from command line arguments
  if (argc != 3 && argc != 6) {
    if( rank == 0 ) printf("ERROR: Number of arguments passed must be 3 or 6."); 
    MPI_Finalize();
    return 0;
  }
  
  m = atoi(argv[1]);
  n = atoi(argv[2]);
  dims[0] = m, dims[1] = n;

  if( (m*n) != total_nodes ) {
    if( rank ==0 ) printf("ERROR: m*n != total number of processes => %d * %d = %d != %d\n", m, n, m*n,total_nodes);
    MPI_Finalize();
    return 0;
  }

  if (argc == 6) {
    MAGNITUDE_UPPER_THRESHOLD = atof(argv[3]);
    DIFF_IN_DISTANCE_THRESHOLD_IN_KM = atof(argv[4]);
    DIFF_IN_MAGNITUDE_THRESHOLD = atof(argv[5]);
  }

	// create cartesian topology for processes
  MPI_Dims_create(total_nodes, ndims, dims);
	if(rank==0) printf("Root Rank: %d. Comm Size: %d: Grid Dimension = [%d x %d] \n",rank,total_nodes,dims[0],dims[1]);

  // Initialise MPI virtual topology (Cartesian)
  int ierr = 0;
  ierr = MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm);
  if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);

  // Find my coordinates in the cartesian communicator group
  MPI_Cart_coords(comm, rank, ndims, coord);
  // Use my cartesian coordinates to find my rank in cartesian group
  MPI_Cart_rank(comm, coord, &cart_rank);

  generate(&newReading);
  printf("Rank %d => ", rank);
  printReading(&newReading);

  MPI_Cart_shift(comm, SHIFT_ROW, DISP, &top_rank, &bottom_rank);
  MPI_Cart_shift(comm, SHIFT_COL, DISP, &left_rank, &right_rank);
  
  if (newReading.mag > MAGNITUDE_UPPER_THRESHOLD) {
    compare_readings = 1;
  }

	// printf("Global rank: %d. Cart rank: %d. Coord: (%d, %d). Left: %d. Right: %d. Top: %d. Bottom: %d\n", rank, cart_rank, coord[0], coord[1], left_rank, right_rank, top_rank, bottom_rank);
  
  struct adj_nodes_arg_struct args;
  struct Sensor readingsT, readingsB, readingsL, readingsR;

  args.pReadingsT = &readingsT;
  args.pReadingsB = &readingsB;
  args.pReadingsL = &readingsL;
  args.pReadingsR = &readingsR;

  pthread_t adj_nodes_comm_t;
  pthread_create(&adj_nodes_comm_t, 0, AdjNodesCommFunc, (void*) &args);
  pthread_join(adj_nodes_comm_t, NULL);

  readingsT = *args.pReadingsT;
  readingsB = *args.pReadingsB;
  readingsL = *args.pReadingsL;
  readingsR = *args.pReadingsR;

  // if(rank == 0){
  //   if(compare_readings_T == 1) {
  //     printf("ReadingsT: \n");
  //     printReading(&readingsT);
  //     printf("\n");
  //   }
    
  //   if(compare_readings_B == 1) {
  //     printf("ReadingsB: \n");
  //     printReading(&readingsB);
  //     printf("\n");
  //   }

  //   if(compare_readings_L == 1) {
  //     printf("ReadingsL: \n");
  //     printReading(&readingsL);
  //     printf("\n");
  //   }

  //   if(compare_readings_R == 1) {
  //     printf("ReadingsR: \n");
  //     printReading(&readingsR);
  //     printf("\n");
  //   }
  // }

  // Receive requested readings
  if (compare_readings == 1) {
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

	MPI_Finalize();
  return 0;
}

/* Send/receive request to adjacent nodes
   Writes to readingsT, readingsB, readingsL, readingsR in main() scope
*/
void* AdjNodesCommFunc(void* pArguments) {
  struct adj_nodes_arg_struct *pArgs = pArguments;

  MPI_Request send_request[4];
  MPI_Request receive_request[4];
  MPI_Status send_status[4];
  MPI_Status receive_status[4];

  // Request reading from adjacent nodes
  MPI_Isend(&compare_readings, 1, MPI_INT, top_rank, 0, comm, &send_request[0]);
  MPI_Isend(&compare_readings, 1, MPI_INT, bottom_rank, 0, comm, &send_request[1]);
  MPI_Isend(&compare_readings, 1, MPI_INT, left_rank, 0, comm, &send_request[2]);
  MPI_Isend(&compare_readings, 1, MPI_INT, right_rank, 0, comm, &send_request[3]);
  
  MPI_Irecv(&compare_readings_T, 1, MPI_INT, top_rank, 0, comm, &receive_request[0]);
  MPI_Irecv(&compare_readings_B, 1, MPI_INT, bottom_rank, 0, comm, &receive_request[1]);
  MPI_Irecv(&compare_readings_L, 1, MPI_INT, left_rank, 0, comm, &receive_request[2]);
  MPI_Irecv(&compare_readings_R, 1, MPI_INT, right_rank, 0, comm, &receive_request[3]);

  MPI_Waitall(4, send_request, send_status);
  MPI_Waitall(4, receive_request, receive_status);

  // Send readings
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

  if (compare_readings == 1) {
    if(top_rank >= 0) {
      MPI_Recv(pArgs->pReadingsT, 1, MPI_READING, top_rank, 0, comm, &status);
    }
    if(bottom_rank >= 0) {
      MPI_Recv(pArgs->pReadingsB, 1, MPI_READING, bottom_rank, 0, comm, &status);
    }
    if(left_rank >= 0) {
      MPI_Recv(pArgs->pReadingsL, 1, MPI_READING, left_rank, 0, comm, &status);
    }
    if(right_rank >= 0) {
      MPI_Recv(pArgs->pReadingsR, 1, MPI_READING, right_rank, 0, comm, &status);
    }
  }

  MPI_Type_free(&MPI_READING);
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