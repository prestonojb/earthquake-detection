#include <stdio.h> 
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <unistd.h>
#include "mpi.h"
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "helper.h"
#include "activity_2.h"
#include "activity_3.h"

#define LATITUDE_LOWER_BOUND -90
#define LATITUDE_UPPER_BOUND 90
#define LONGITUDE_LOWER_BOUND -180
#define LONGITUDE_UPPER_BOUND 180
#define MAGNITUDE_LOWER_BOUND 0
#define MAGNITUDE_UPPER_BOUND 5
#define DEPTH_UPPER_BOUND 700

#define SHIFT_ROW 0
#define SHIFT_COL 1
#define DISP 1

#define READING_INTERVAL_IN_S 5

void generateNodeReading(struct Sensor* reading, int rank);
void printReading(struct Sensor* reading);

void* AdjNodesCommFunc(void* pArgs);

bool areMatchingReadings(struct Sensor* readingA, struct Sensor* readingB);
bool areMatchingLocations(float latA, float longA, float latB, float longB);
bool areMatchingMagnitudes(float magnitudeA, float magnitudeB);

float distance(float lat1, float lon1, float lat2, float lon2);

double deg2rad(double);
double rad2deg(double);

struct adj_nodes_arg_struct {
  struct Sensor* pReadingT;
  struct Sensor* pReadingB;
  struct Sensor* pReadingL;
  struct Sensor* pReadingR;
};

float MAGNITUDE_UPPER_THRESHOLD;
float DIFF_IN_DISTANCE_THRESHOLD_IN_KM;
float DIFF_IN_MAGNITUDE_THRESHOLD;

int left_rank, right_rank, top_rank, bottom_rank;
int compare_adj_readings = 0;
struct Sensor currReading;
int compare_readingT = -1, compare_readingB = -1, compare_readingL = -1, compare_readingR = -1;

// Initialise MPI variables
MPI_Status status;
MPI_Comm comm2D;

/**
 * The entry point to the Nodes
 * (previously main)
 * @param argc Argument count
 * @param argv Argument variables
 * @return
 */
int init_nodes(int m, int n, float magnitude_upper_threshold, float diff_in_distance_threshold_in_km, float diff_in_magnitude_threshold, MPI_Comm world_comm, MPI_Comm comm) {
  // Initialise variables
  int ndims = 2;
  int dims[ndims], coord[ndims], wrap_around[ndims];
  int reorder = 1;

  wrap_around[0] = wrap_around[1] = 0;
  dims[0] = m, dims[1] = n;

  MAGNITUDE_UPPER_THRESHOLD = magnitude_upper_threshold;
  DIFF_IN_DISTANCE_THRESHOLD_IN_KM = diff_in_distance_threshold_in_km;
  DIFF_IN_MAGNITUDE_THRESHOLD = diff_in_magnitude_threshold;

  int world_rank, total_processes;
  MPI_Comm_rank(world_comm, &world_rank);
  MPI_Comm_size(world_comm, &total_processes);

  int node_rank, total_nodes;
  MPI_Comm_rank(comm, &node_rank);
  MPI_Comm_size(comm, &total_nodes);

	// create cartesian topology for processes
  MPI_Dims_create(total_nodes, ndims, dims);
	// if(node_rank==0) printf("Comm Size: %d: Grid Dimension = [%d x %d] \n",total_nodes,dims[0],dims[1]);

  // Initialise MPI virtual topology (Cartesian)
  int ierr = 0;
  ierr = MPI_Cart_create(comm, ndims, dims, wrap_around, reorder, &comm2D);
  if(ierr != 0) printf("ERROR[%d] creating CART\n",ierr);

  // Find my coordinates in the cartesian communicator group
  MPI_Cart_coords(comm2D, node_rank, ndims, coord);

  MPI_Cart_shift(comm2D, SHIFT_ROW, DISP, &top_rank, &bottom_rank);
  MPI_Cart_shift(comm2D, SHIFT_COL, DISP, &left_rank, &right_rank);
  
  MPI_Status probe_status;
	int flag = 0;
  int termination_msg_buf;

  MPI_Datatype SensorType;
  defineSensorType(&SensorType);
  
  MPI_Datatype DataLogType;
  defineDataLogType(&DataLogType, SensorType);

  struct DataLog datalog;

  // Generates reading indefinitely until a termination message is received from base station
  while(1) {

    MPI_Iprobe(MPI_ANY_SOURCE, TERMINATION_TAG, world_comm, &flag, &probe_status);
    if(flag == 1) {
			MPI_Recv(&termination_msg_buf, 1, MPI_INT, probe_status.MPI_SOURCE, TERMINATION_TAG, world_comm, &status);
      printf("Node %d received termination message from base station. Terminating... \n", node_rank);
        break;
    }

    generateNodeReading(&currReading, node_rank);
    // printf("Node Rank %d => ", node_rank);
    // printReading(&currReading);
    
    if (currReading.mag > MAGNITUDE_UPPER_THRESHOLD) {
        //printf("Node %d triggered magnitude %f\n", node_rank, currReading.mag);
      compare_adj_readings = 1;
    }

    // printf("Cart rank: %d. Coord: (%d, %d). Left: %d. Right: %d. Top: %d. Bottom: %d\n", node_rank, coord[0], coord[1], left_rank, right_rank, top_rank, bottom_rank);
    
    struct adj_nodes_arg_struct args;
    struct Sensor readingT, readingB, readingL, readingR;

    args.pReadingT = &readingT;
    args.pReadingB = &readingB;
    args.pReadingL = &readingL;
    args.pReadingR = &readingR;

    pthread_t adj_nodes_comm_t;
    pthread_create(&adj_nodes_comm_t, 0, AdjNodesCommFunc, (void*) &args);
    pthread_join(adj_nodes_comm_t, NULL);

    readingT = *args.pReadingT;
    readingB = *args.pReadingB;
    readingL = *args.pReadingL;
    readingR = *args.pReadingR;

    // Receive requested readings
    if (compare_adj_readings == 1) {
      int no_of_matches = 0;
      if(compare_readingT == 1 && areMatchingReadings(&currReading, &readingT)) no_of_matches++;
      if(compare_readingB == 1 && areMatchingReadings(&currReading, &readingB)) no_of_matches++;
      if(compare_readingL == 1 && areMatchingReadings(&currReading, &readingL)) no_of_matches++;
      if(compare_readingR == 1 && areMatchingReadings(&currReading, &readingR)) no_of_matches++;

      if(no_of_matches >= 2) {
        datalog.reporterRank = node_rank;
        datalog.reporterData = currReading;
        datalog.topRank = top_rank;
        datalog.topData = readingT;
        datalog.bottomRank = bottom_rank;
        datalog.bottomData = readingB;
        datalog.leftRank = left_rank;
        datalog.leftData = readingL;
        datalog.rightRank = right_rank;
        datalog.rightData = readingR;

        MPI_Send(&datalog, 1, DataLogType, BASE_STATION, 0, world_comm);
        printf("Sensor node %d is triggered, send report to base station \n", node_rank);
        // printReading(&currReading);
      }
    }

    sleep(READING_INTERVAL_IN_S);
  }

  MPI_Type_free(&SensorType);
  MPI_Type_free(&DataLogType);
  MPI_Comm_free(&comm2D);
  return 0;
}

/**
 * @brief Send/receive request to adjacent nodes
 * Writes to readingT, readingB, readingL, readingR
 */

void* AdjNodesCommFunc(void* pArguments) {
  struct adj_nodes_arg_struct *pArgs = pArguments;

  MPI_Request send_request[4];
  MPI_Request receive_request[4];
  MPI_Status send_status[4];
  MPI_Status receive_status[4];

  // Request reading from adjacent nodes
  MPI_Isend(&compare_adj_readings, 1, MPI_INT, top_rank, 0, comm2D, &send_request[0]);
  MPI_Isend(&compare_adj_readings, 1, MPI_INT, bottom_rank, 0, comm2D, &send_request[1]);
  MPI_Isend(&compare_adj_readings, 1, MPI_INT, left_rank, 0, comm2D, &send_request[2]);
  MPI_Isend(&compare_adj_readings, 1, MPI_INT, right_rank, 0, comm2D, &send_request[3]);
  
  MPI_Irecv(&compare_readingT, 1, MPI_INT, top_rank, 0, comm2D, &receive_request[0]);
  MPI_Irecv(&compare_readingB, 1, MPI_INT, bottom_rank, 0, comm2D, &receive_request[1]);
  MPI_Irecv(&compare_readingL, 1, MPI_INT, left_rank, 0, comm2D, &receive_request[2]);
  MPI_Irecv(&compare_readingR, 1, MPI_INT, right_rank, 0, comm2D, &receive_request[3]);

  MPI_Waitall(4, send_request, send_status);
  MPI_Waitall(4, receive_request, receive_status);

  // Send readings
  MPI_Datatype SensorType;
  defineSensorType(&SensorType);

  // Send (blocking) readings to adjacent nodes (if requested)
  if (compare_readingT == 1) {
    MPI_Send(&currReading, 1, SensorType, top_rank, 0, comm2D);
  }
  if (compare_readingB == 1) {
    MPI_Send(&currReading, 1, SensorType, bottom_rank, 0, comm2D);
  }
  if (compare_readingL == 1) {
    MPI_Send(&currReading, 1, SensorType, left_rank, 0, comm2D);
  }
  if (compare_readingR == 1) {
    MPI_Send(&currReading, 1, SensorType, right_rank, 0, comm2D);
  }

  if (compare_adj_readings == 1) {
    if(top_rank >= 0) {
      MPI_Recv(pArgs->pReadingT, 1, SensorType, top_rank, 0, comm2D, &status);
    }
    if(bottom_rank >= 0) {
      MPI_Recv(pArgs->pReadingB, 1, SensorType, bottom_rank, 0, comm2D, &status);
    }
    if(left_rank >= 0) {
      MPI_Recv(pArgs->pReadingL, 1, SensorType, left_rank, 0, comm2D, &status);
    }
    if(right_rank >= 0) {
      MPI_Recv(pArgs->pReadingR, 1, SensorType, right_rank, 0, comm2D, &status);
    }
  }

  MPI_Type_free(&SensorType);
  return 0;
}

/**
 * @brief Generate sensor node reading 
 * 
 * @param reading 
 */
void generateNodeReading(struct Sensor* reading, int rank)
{
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  reading->year = tm.tm_year + 1900;
  reading->month = tm.tm_mon + 1;
  reading->day = tm.tm_mday;
  reading->hour = tm.tm_hour;
  reading->minute = tm.tm_min;
  reading->second = tm.tm_sec;

  reading->lat = float_rand_seed(LATITUDE_LOWER_BOUND, LATITUDE_UPPER_BOUND, rank);
  reading->lon = float_rand_seed(LONGITUDE_LOWER_BOUND, LONGITUDE_UPPER_BOUND, rank);
  reading->mag = float_rand_seed(MAGNITUDE_LOWER_BOUND, MAGNITUDE_UPPER_BOUND, rank);
  reading->depth = float_rand_seed(0, DEPTH_UPPER_BOUND, rank);
    //printf("Node %d: %f, %f, %f, %f\n", rank, reading->lat, reading->lon, reading->mag, reading->depth);
}

/**
 * @brief Returns true if both readings match, else false
 * 
 * @param readingA 
 * @param readingB 
 * @return true 
 * @return false 
 */
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