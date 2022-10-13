/*
 * This file serves as the main entry to the program
 * Contains the code for the base station
 */

#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <pthread.h>
#include <math.h>
#include "activity_1.h"
#include "activity_2.h"
#include "helper.h"
#include <unistd.h>

#define INTERVAL 5

void update();
void defineSensorType(struct Sensor sensor, MPI_Datatype* SensorType);
void defineDataLogType(struct DataLog dataLog, MPI_Datatype* DataLogType, MPI_Datatype SensorType);
int saveLog(int conclusion, int intervalCount, struct DataLog n, struct Sensor b);
void exitBase(int sentinelValue);
int checkSentinel();



int main(int argc, char* argv[]) {

    // Initialise MPI variables
    MPI_Status status;
    MPI_Comm comm;

    // Initialise variables
    int rank, total_nodes, provided;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided );
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &total_nodes);

    if (rank == 0) update();
    else init_nodes(argc, argv, rank, total_nodes - 1);

    return 0;
}

/**
 * Runs continuously with a fixed delay.
 * Stops when encountered a sentinel value.
 */
void update() {

    int intervalCount = 0;

    struct Sensor sensor;
    MPI_Datatype SensorType;
    defineSensorType(sensor, &SensorType);

    struct DataLog dataLog;
    MPI_Datatype DataLogType;
    defineDataLogType(dataLog, &DataLogType, SensorType);

    // Loop until sentinel value encountered
    int sentinelVal = checkSentinel();
    while (sentinelVal == 0) {

        intervalCount++;

        printf("Listening for seismic activity...\n");
        // Todo: use POSIX to send and receive
        MPI_Recv(&dataLog, 1, DataLogType, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Received seismic data! Comparing with balloon sensor.\n");

        // Todo: check balloon shared array

        int conclusive = 1;
        saveLog(conclusive, intervalCount, dataLog, sensor);

        sleep(INTERVAL);
        sentinelVal = checkSentinel();
    }

    /* Clean up the type */
    MPI_Type_free( &SensorType );
    MPI_Finalize( );

    exitBase(sentinelVal);
}

/**
 * Defines an MPI datatype of type Sensor
 * @param sensor
 * @param SensorType
 */
void defineSensorType(struct Sensor sensor, MPI_Datatype* SensorType) {
    const int vars = 10;
    MPI_Datatype type[2] = { MPI_INT, MPI_FLOAT };
    int blockLen[2] = {6, 4 };
    MPI_Aint disp[vars];

    // Get address for all variables
    MPI_Get_address(&sensor.year, &disp[0]);
    MPI_Get_address(&sensor.month, &disp[1]);
    MPI_Get_address(&sensor.day, &disp[2]);
    MPI_Get_address(&sensor.hour, &disp[3]);
    MPI_Get_address(&sensor.minute, &disp[4]);
    MPI_Get_address(&sensor.second, &disp[5]);
    MPI_Get_address(&sensor.lat, &disp[6]);
    MPI_Get_address(&sensor.lon, &disp[7]);
    MPI_Get_address(&sensor.mag, &disp[8]);
    MPI_Get_address(&sensor.depth, &disp[9]);

    // Fix displacement
    for (int i = 1; i < vars; i++)
        disp[i] = disp[i] - disp[i-1];
    disp[0] = 0;

    // Create MPI datatype
    MPI_Type_create_struct(vars, blockLen, disp, type, SensorType);
    MPI_Type_commit(SensorType);
}

void defineDataLogType(struct DataLog dataLog, MPI_Datatype* DataLogType, MPI_Datatype SensorType) {
    const int vars = 10;
    MPI_Datatype type[2] = { MPI_INT, SensorType };
    int blockLen[2] = {5, 5 };
    MPI_Aint disp[vars];

    // Get address for all variables
    MPI_Get_address(&dataLog.reporterRank, &disp[0]);
    MPI_Get_address(&dataLog.reporterData, &disp[1]);
    MPI_Get_address(&dataLog.topRank, &disp[2]);
    MPI_Get_address(&dataLog.topData, &disp[3]);
    MPI_Get_address(&dataLog.bottomRank, &disp[4]);
    MPI_Get_address(&dataLog.bottomData, &disp[5]);
    MPI_Get_address(&dataLog.leftRank, &disp[6]);
    MPI_Get_address(&dataLog.leftData, &disp[7]);
    MPI_Get_address(&dataLog.rightRank, &disp[8]);
    MPI_Get_address(&dataLog.rightData, &disp[9]);

    // Fix displacement
    for (int i = 1; i < vars; i++)
        disp[i] = disp[i] - disp[i-1];
    disp[0] = 0;

    // Create MPI datatype
    MPI_Type_create_struct(vars, blockLen, disp, type, DataLogType);
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
void exitBase(int sentinelValue) {


    if (sentinelValue == 1)
        printf("Sentinel value detected. Quitting!\n");
    else
        printf("Encountered error reading sentinel file");

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
