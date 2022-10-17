/*
 * This file serves as the main entry to the program
 * Contains the code for the base station
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
#include <pthread.h>
#include <math.h>
#include <unistd.h>
#include "helper.h"
#include "activity_1.h"
#include "activity_2.h"

#define READING_INTERVAL_IN_S 1

void update();
void createBalloonPosix();
void defineSensorType(MPI_Datatype* SensorType);
void defineDataLogType(MPI_Datatype* DataLogType, MPI_Datatype SensorType);
int saveLog(int conclusion, int intervalCount, struct DataLog n, struct Sensor b);
void exitBase(MPI_Comm world_comm);
int checkSentinel();
void printColNode(FILE* f, int id, float lat, float lon, float dist, float mag, float depth);
void* terminationToNodesComm();
void* recvDataLogFromNodesCommFunc(void* pArg);

float MAGNITUDE_UPPER_THRESHOLD = DEFAULT_MAGNITUDE_UPPER_THRESHOLD;
float DIFF_IN_DISTANCE_THRESHOLD_IN_KM = DEFAULT_DIFF_IN_DISTANCE_THRESHOLD_IN_KM;
float DIFF_IN_MAGNITUDE_THRESHOLD = DEFAULT_DIFF_IN_MAGNITUDE_THRESHOLD;

// Shared Balloon sensor readings
#define SIZE 10
struct Sensor readings[SIZE];
pthread_t balloon_comm[2];

int main(int argc, char* argv[]) {
    // Create empty log file
    FILE *f = fopen("log.txt","w");
    fclose(f);

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

    return 0;
}

/**
 * Create a POSIX thread to start the balloon sensor
 */
void createBalloonPosix() {
    pthread_create(&balloon_comm[0], 0, startBalloon, (void*) readings);
}

/**
 * @brief Generate balloon reading, based on set magnitude threshold
 * 
 * @param reading 
 */
void generateBalloonReading(struct Sensor* reading)
{
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    reading->year = tm.tm_year + 1900;
    reading->month = tm.tm_mon + 1;
    reading->day = tm.tm_mday;
    reading->hour = tm.tm_hour;
    reading->minute = tm.tm_min;
    reading->second = tm.tm_sec;

    reading->lat = float_rand(-20, -10);
    reading->lon = float_rand(150, 170);
    reading->mag = float_rand(MAGNITUDE_UPPER_THRESHOLD, 9);
    reading->depth = float_rand(4, 10);
}

/**
 * Runs continuously with a fixed delay.
 * Stops when encountered a sentinel value.
 */
void update(MPI_Comm world_comm) {
    int intervalCount = 0;
    struct DataLog dataLog;

    // Loop until sentinel value encountered
    int sentinelVal = checkSentinel();
    while (sentinelVal == 0) {
        intervalCount++;

        printf("Listening for seismic activity...\n");

        pthread_t recv_datalog_from_nodes_comm_t;
        pthread_create(&recv_datalog_from_nodes_comm_t, 0, recvDataLogFromNodesCommFunc, &dataLog);
        pthread_join(recv_datalog_from_nodes_comm_t, NULL);

        // Retrieving last value from shared balloon readings array
        int finalVal = queue_head;  // From "helper.h"
        struct Sensor balloonReading = readings[finalVal-1];

        int conclusive = areMatchingMagnitudes(dataLog.reporterData.mag, balloonReading.mag);

        saveLog(conclusive, intervalCount, dataLog, balloonReading);
        printf("Base station logs alert from node %d to log.txt file. \n", dataLog.reporterRank);

        sleep(READING_INTERVAL_IN_S);
        sentinelVal = checkSentinel();
    }

    exitBase(MPI_COMM_WORLD);
}

/**
 * @brief Base station to receive alert report from sensor nodes
 * 
 * @param pArg 
 * @return void* 
 */
void* recvDataLogFromNodesCommFunc(void* pArg) {
    struct DataLog *pDataLog = pArg;

    struct Sensor sensor;
    MPI_Datatype SensorType;
    defineSensorType(&SensorType);

    MPI_Datatype DataLogType;
    defineDataLogType(&DataLogType, SensorType);

    MPI_Recv(pDataLog, 1, DataLogType, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Base station received seismic data from node %d, comparing with balloon sensor... \n", pDataLog->reporterRank);

    /* Clean up the type */
    MPI_Type_free( &SensorType );
    MPI_Type_free( &DataLogType );
    return 0;
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

    fprintf(f, "Logged Time:\t%d:%d:%d  %d-%d-%d\n",
            tm.tm_hour, tm.tm_min, tm.tm_sec,
            tm.tm_mday, tm.tm_mon + 1, tm.tm_year + 1900);

    struct Sensor s = n.reporterData;
    fprintf(f, "Alerted Time:\t%d:%d:%d  %d-%d-%d\n",
            s.hour, s.minute, s.second, s.day, s.month, s.year);

    char* c_result[2] = {"Inconclusive", "Conclusive"};
    fprintf(f, "Alert is %s\n\n", c_result[conclusion]);

    fprintf(f, "NODE SEISMIC REPORT:\n");
    fprintf(f, "ID\tCoordinates\t\t\t\t\t\tMagnitude\tDepth\n");
    fprintf(f, "%d\t(%.2f, %.2f)  \t\t\t\t%.2f\t\t%.2f\n\n", n.reporterRank, s.lat, s.lon, s.mag, s.depth);

    fprintf(f, "ID\tCoordinates\t\t\tDistance\tMagnitude\tDepth\n");

    struct Sensor s_t = n.topData;
    float distSensor = distance(s.lat, s.lon, s_t.lat, s_t.lon);
    printColNode(f, n.topRank, s_t.lat, s_t.lon, distSensor, s_t.mag, s_t.depth);

    struct Sensor s_b = n.bottomData;
    distSensor = distance(s.lat, s.lon, s_b.lat, s_b.lon);
    printColNode(f, n.bottomRank, s_b.lat, s_b.lon, distSensor, s_b.mag, s_b.depth);

    struct Sensor s_l = n.leftData;
    distSensor = distance(s.lat, s.lon, s_l.lat, s_l.lon);
    printColNode(f, n.leftRank, s_l.lat, s_l.lon, distSensor, s_l.mag, s_l.depth);

    struct Sensor s_r = n.rightData;
    distSensor = distance(s.lat, s.lon, s_r.lat, s_r.lon);
    printColNode(f, n.rightRank, s_r.lat, s_r.lon, distSensor, s_r.mag, s_r.depth);

    fprintf(f, "\nBALLOON SEISMIC REPORT:\n");
    fprintf(f, "Report Time: %d:%d:%d  %d-%d-%d\n",
            b.hour, b.minute, b.second, b.day, b.month, b.year);

    float distBalloon2Node = distance(b.lat, b.lon, s.lat, s.lon);
    fprintf(f, "Coordinates: (%.2f, %.2f)\n", b.lat, b.lon);
    fprintf(f, "Distance from Balloon to Reporting Node: %.2f\n", distBalloon2Node);
    fprintf(f, "Magnitude: %.2f\n", b.mag);
    fprintf(f, "Magnitude difference with Reporting Node: %.2f\n\n", fabs(s.mag - b.mag));

    float time = tm.tm_sec - s.second;
    fprintf(f, "Communication Time (in seconds): %.3fs\n", time);
    fprintf(f, "Total messages sent by Node to Base: 1\n");
    fprintf(f, "Coordinate Threshold: %.2f\n", DIFF_IN_DISTANCE_THRESHOLD_IN_KM);
    fprintf(f, "Magnitude Difference Threshold: %.2f\n", DIFF_IN_MAGNITUDE_THRESHOLD);
    fprintf(f, "Earthquake Magnitude Threshold: %.2f\n", MAGNITUDE_UPPER_THRESHOLD);

    fprintf(f,"%s\n%s\n\n", line, line);
    fclose(f);

    return 0;
}

void printColNode(FILE* f, int id, float lat, float lon, float dist, float mag, float depth) {
    if (id == -2) {
        fprintf(f, "-\t(----, ----)\t\t-\t\t\t-\t\t\t-\n");
    }
    else {
        fprintf(f, "%d\t(%.2f, %.2f)  \t%.2f\t\t%.2f\t\t%.2f\n", id, lat, lon, dist, mag, depth);
    }
}
struct termination_to_nodes_arg_struct {
    int total_nodes;
};

/**
 * Exit the program.
 * Send an exit signal to all nodes before quitting.
 * @param sentinelValue Exit status
 */
void exitBase(MPI_Comm world_comm) {
    printf("Balloon sensor node exiting...");
    // Comm Balloon to quit
    int message = 1;
    pthread_create(&balloon_comm[1], 0, receiveMessage, &message);
    pthread_join(balloon_comm[0], NULL);
    pthread_join(balloon_comm[1], NULL);

    printf("Sensor Nodes exiting...");
    // Send termination message to sensor nodes
    int total_processes;
    MPI_Comm_size(world_comm, &total_processes);
    int total_nodes = total_processes - 1;
    
    struct termination_to_nodes_arg_struct args;

    args.total_nodes = total_nodes;

    pthread_t termination_to_nodes_comm_t;
    pthread_create(&termination_to_nodes_comm_t, 0, terminationToNodesComm, (void*) &args);
    pthread_join(termination_to_nodes_comm_t, NULL);

    printf("Base has quit successfully!");
    MPI_Finalize();
}

/**
 * @brief Sends termination message to sensor nodes
 * 
 * @param pArguments 
 * @return void* 
 */
void* terminationToNodesComm(void *pArguments) {
    struct termination_to_nodes_arg_struct *pArgs = pArguments;
    int total_nodes = pArgs->total_nodes;

    MPI_Request send_request[total_nodes];
    MPI_Status send_status[total_nodes];
    
    int termination_msg = 0;

    for(int i=0; i < total_nodes; i++) {
        int node_rank = i+1; // offset base station rank
        
        printf("Send termination message to node %d \n", node_rank);
        MPI_Isend(&termination_msg, 1, MPI_INT, node_rank, TERMINATION_TAG, MPI_COMM_WORLD, &send_request[i]);
    }
    
    // printf("Waiting to receive send ACK of all termination messages from sensor nodes... \n");
    MPI_Waitall(total_nodes, send_request, send_status);

    return 0;
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
