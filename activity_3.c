/*
 * This file serves as the main entry to the program
 * Contains the code for the base station
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <pthread.h>
#include <math.h>
#include "activity_1.h"
#include "activity_2.h"
#include "helper.h"
#include <unistd.h>

#define INTERVAL 5

#define DEFAULT_MAGNITUDE_UPPER_THRESHOLD 0
#define DEFAULT_DIFF_IN_DISTANCE_THRESHOLD_IN_KM 100000
#define DEFAULT_DIFF_IN_MAGNITUDE_THRESHOLD 1000

#define TERMINATION_TAG 10

void update();
void createBalloonPosix();
void defineSensorType(MPI_Datatype* SensorType);
void defineDataLogType(MPI_Datatype* DataLogType, MPI_Datatype SensorType);
int saveLog(int conclusion, int intervalCount, struct DataLog n, struct Sensor b);
void exitBase(MPI_Comm world_comm);
int checkSentinel();

float MAGNITUDE_UPPER_THRESHOLD = DEFAULT_MAGNITUDE_UPPER_THRESHOLD;
float DIFF_IN_DISTANCE_THRESHOLD_IN_KM = DEFAULT_DIFF_IN_DISTANCE_THRESHOLD_IN_KM;
float DIFF_IN_MAGNITUDE_THRESHOLD = DEFAULT_DIFF_IN_MAGNITUDE_THRESHOLD;

// Shared Balloon sensor readings
#define SIZE 10
struct Sensor readings[SIZE];

int main(int argc, char* argv[]) {
    // Initialise variables
    int world_rank, total_processes, provided;
    
    // Initialise MPI variables
    MPI_Status status;

    // Seed rand() function with current timestamp and rank
    srand(time(NULL) + world_rank*100);

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided );
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_processes);

    // Get dimensions of topology from command line arguments
    if (argc != 3 && argc != 6) {
        if( world_rank == 0 ) printf("ERROR: Number of arguments passed must be 3 or 6."); 
        MPI_Finalize();
        return 0;
    }

    int m = atoi(argv[1]);
    int n = atoi(argv[2]);

    if( (m*n) != total_processes - 1 ) {
        if( world_rank == 0 ) printf("ERROR: m*n != total number of processes - 1 (base station) => %d * %d = %d != %d\n", m, n, m*n,total_processes-1);
        MPI_Finalize();
        return 0;
    }

    if (argc == 6) {
        MAGNITUDE_UPPER_THRESHOLD = atof(argv[3]);
        DIFF_IN_DISTANCE_THRESHOLD_IN_KM = atof(argv[4]);
        DIFF_IN_MAGNITUDE_THRESHOLD = atof(argv[5]);
    }

    MPI_Comm nodes_comm;
    // Split into base station and sensor nodes
    MPI_Comm_split( MPI_COMM_WORLD, world_rank == 0, 0, &nodes_comm);

    if (world_rank == 0) {
        createBalloonPosix(readings);
        update(MPI_COMM_WORLD);
    } else {
        init_nodes(m, n, MAGNITUDE_UPPER_THRESHOLD, DIFF_IN_DISTANCE_THRESHOLD_IN_KM, DIFF_IN_MAGNITUDE_THRESHOLD, MPI_COMM_WORLD, nodes_comm);
    }

    // if (rank == 0) update();
    // else if (rank == total_nodes - 1) init_balloon();
    // else init_nodes(argc, argv, rank, total_nodes - 2);

    // if(world_rank != 0) {
    // }

    // int errorcode;
    // MPI_Abort(MPI_COMM_WORLD, errorcode);
    printf("byee from rank %d! \n", world_rank);

    MPI_Finalize();
    return 0;
}

/**
 * Create a POSIX thread to start the balloon sensor
 */
void createBalloonPosix() {
    pthread_t balloon_comm;
    pthread_create(&balloon_comm, 0, startBalloon, (void*) readings);
}

/**
 * Runs continuously with a fixed delay.
 * Stops when encountered a sentinel value.
 */
void update(MPI_Comm world_comm) {
    int intervalCount = 0;

    struct Sensor sensor;
    MPI_Datatype SensorType;
    defineSensorType(&SensorType);

    struct DataLog dataLog;
    MPI_Datatype DataLogType;
    defineDataLogType(&DataLogType, SensorType);

    // Loop until sentinel value encountered
    int sentinelVal = checkSentinel();
    while (sentinelVal == 0) {

        intervalCount++;

        printf("Listening for seismic activity...\n");
        
        MPI_Recv(&dataLog, 1, DataLogType, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Received seismic data from node %d! Comparing with balloon sensor.\n", dataLog.reporterRank);
        
        
        // Todo: use POSIX to send and receive


        // Retrieving last value from shared balloon readings array
        int finalVal = queue_head;  // From "helper.h"
        struct Sensor balloonReading = readings[finalVal];

        // Todo: Compare balloon and nodes reading

        int conclusive = 1;
        saveLog(conclusive, intervalCount, dataLog, sensor);

        sleep(INTERVAL);
        sentinelVal = checkSentinel();
    }

    /* Clean up the type */
    MPI_Type_free( &SensorType );

    exitBase(MPI_COMM_WORLD);
}

/**
 * Defines an MPI datatype of type Sensor
 * @param SensorType
 */
void defineSensorType(MPI_Datatype* SensorType) {
    const int readingSize = 10;
    int blocklengths[10] = {1,1,1,1,1,1,1,1,1,1};
    MPI_Datatype types[10] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
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

    MPI_Type_create_struct(readingSize, blocklengths, offsets, types, SensorType);
    MPI_Type_commit(SensorType);
}

void defineDataLogType(MPI_Datatype* DataLogType, MPI_Datatype SensorType) {
    const int readingSize = 10;
    int blocklengths[10] = {1,1,1,1,1,1,1,1,1,1};
    MPI_Datatype types[10] = {MPI_INT, SensorType,
                              MPI_INT, SensorType,
                              MPI_INT, SensorType,
                              MPI_INT, SensorType,
                              MPI_INT, SensorType};
    MPI_Aint offsets[10];

    offsets[0] = offsetof(struct DataLog, reporterRank);
    offsets[1] = offsetof(struct DataLog, reporterData);
    offsets[2] = offsetof(struct DataLog, topRank);
    offsets[3] = offsetof(struct DataLog, topData);
    offsets[4] = offsetof(struct DataLog, bottomRank);
    offsets[5] = offsetof(struct DataLog, bottomData);
    offsets[6] = offsetof(struct DataLog, leftRank);
    offsets[7] = offsetof(struct DataLog, leftData);
    offsets[8] = offsetof(struct DataLog, rightRank);
    offsets[9] = offsetof(struct DataLog, rightData);

    MPI_Type_create_struct(readingSize, blocklengths, offsets, types, DataLogType);
    MPI_Type_commit(DataLogType);
}


/**
 * Save the log in a text file
 * @param conclusive Whether the reporting is a conclusive alert (1) or not (0)
 * @param intervalCount Number of intervals finished when reported
 * @param n The nodes sensor data
 * @param b The balloon sensor data
 * @return 0: Success
 *        -1: Error
 */
int saveLog(int conclusion, int intervalCount, struct DataLog n, struct Sensor b) {
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char* line = "------------------------------------";
    FILE *f;
    f=fopen("log.txt","a");
    if (!f) return 1;
    fprintf(f,"%s\n%s\n", line, line);

    fprintf(f, "Iterations: %d\n", intervalCount);

    fprintf(f, "Logged Time: %d:%d:%d  %d-%d-%d\n",
            tm.tm_sec, tm.tm_min, tm.tm_hour,
            tm.tm_mday, tm.tm_mon + 1, tm.tm_year + 1900);

    struct Sensor s = n.reporterData;
    fprintf(f, "Alerted Time: %d:%d:%d  %d-%d-%d\n",
            s.second, s.minute, s.hour, s.day, s.month, s.year);

    char* c_result[2] = {"Inconclusive", "Conclusive"};
    fprintf(f, "Alert is %s\n\n", c_result[conclusion]);

    fprintf(f, "NODE SEISMIC REPORT:\n");
    fprintf(f, "ID\tCoordinates\t\t\tMagnitude\n");
    fprintf(f, "%d\t(%f, %f)\t\t\t%f\n\n", n.reporterRank, s.lat, s.lon, s.mag);

    // Todo: Fill in Distance
    fprintf(f, "ID\tCoordinates\tDistance\tMagnitude\n");

    struct Sensor s_t = n.topData;
    fprintf(f, "%d\t(%f, %f)\t\t\t\t%f\n", n.topRank, s_t.lat, s_t.lon, s_t.mag);

    struct Sensor s_b = n.bottomData;
    fprintf(f, "%d\t(%f, %f)\t\t\t\t%f\n", n.bottomRank, s_b.lat, s_b.lon, s_b.mag);

    struct Sensor s_l = n.leftData;
    fprintf(f, "%d\t(%f, %f)\t\t\t\t%f\n", n.leftRank, s_l.lat, s_l.lon, s_l.mag);

    struct Sensor s_r = n.rightData;
    fprintf(f, "%d\t(%f, %f)\t\t\t\t%f\n\n", n.rightRank, s_r.lat, s_r.lon, s_r.mag);

    fprintf(f, "BALLOON SEISMIC REPORT:\n");
    fprintf(f, "Report Time: %d:%d:%d  %d-%d-%d\n",
            b.second, b.minute, b.hour, b.day, b.month, b.year);

    fprintf(f, "Coordinates: (%f, %f)\n", b.lat, b.lon);
    // Todo: Distance from Balloon to Reporting Node
    fprintf(f, "Distance from Balloon to Reporting Node: -\n");
    fprintf(f, "Magnitude: %f\n", b.mag);
    fprintf(f, "Magnitude difference with Reporting Node: %f\n\n", fabs(s.mag - b.mag));

    // Todo: Communication Time
    fprintf(f, "Communication Time (in seconds): - s\n");
    fprintf(f, "Total messages sent by Node to Base: 1\n");
    fprintf(f, "Coordinate Threshold: 20\n"); // Todo: hardcode this value
    fprintf(f, "Magnitude Difference Threshold: 0.5\n"); // Todo: hardcode this value
    fprintf(f, "Earthquake Magnitude Threshold: 2.5\n"); // Todo: hardcode this value

    fprintf(f,"%s\n%s\n\n", line, line);
    fclose(f);

    return 0;
}

/**
 * Exit the program.
 * Send an exit signal to all nodes before quitting.
 * @param sentinelValue Exit status
 */
void exitBase(MPI_Comm world_comm) {

    // Comm Balloon to quit
    pthread_t balloon_comm;
    int message = 1;
    pthread_create(&balloon_comm, 0, receiveMessage, &message);
    pthread_join(balloon_comm, NULL);

    // Send termination message to sensor nodes
    int total_processes;
    MPI_Comm_size(world_comm, &total_processes);
    int total_nodes = total_processes - 1;
    
    int termination_msg = 0;

    MPI_Request send_request[total_nodes];
    MPI_Status send_status[total_nodes];

    for(int i=0; i < total_nodes; i++) {
        int node_rank = i+1; // offset base station rank
        
        printf("Send termination message to node %d \n", node_rank);
        MPI_Isend(&termination_msg, 1, MPI_INT, node_rank, TERMINATION_TAG, world_comm, &send_request[i]);
    }
    
    printf("Waiting to receive send ACK of all termination messages from sensor nodes... \n");
    MPI_Waitall(total_nodes, send_request, send_status);
}

/**
 * Check if program should stop.
 * This func will check for a sentinel value
 * in a file to determine exit.
 * @return 0: continue,
 *         1: exit,
 *        -1: file error,
 */
int checkSentinel() {
    FILE *f;
    f = fopen("sentinel.txt", "r");
    if (!f) return -1;
    char s[2];
    fgets(s,2,f);
    fclose(f);
    if (!strcmp(s, "0"))
        return 0;
    else return 1;
}
