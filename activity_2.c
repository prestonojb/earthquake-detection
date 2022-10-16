#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include "helper.h"

#define READING_INTERVAL_IN_S 2

void generateBalloonReading(struct Sensor* reading);
void printReading(struct Sensor* reading);
void printQueue(struct Sensor arr[], int length);

int shutdown = 0;

/**
 * Started on a thread by the Base
 */
void* startBalloon(void* pArg) {
    struct Sensor* sharedReadings = pArg;
    while (shutdown == 0) {
        struct Sensor newReading;
        generateBalloonReading(&newReading);
        enqueue(sharedReadings, QUEUE_SIZE, newReading);

        // printReading(&newReading);
        // printQueue(sharedReadings, QUEUE_SIZE);
        sleep(READING_INTERVAL_IN_S);
    }
    return 0;
}

/**
 * Receive a message from Base
 * @param pArg
 * @return
 */
void* receiveMessage(void *pArg) {
    int* p = (int*) pArg;
    int message = *p;

    if (message == 1) {
        shutdown = 1;
    }
    return 0;
}

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
    reading->mag = float_rand(1, 5);
    reading->depth = float_rand(4, 10);
}

void printQueue(struct Sensor arr[], int length)
{
    printf(" YYYY\t| MM\t| DD\t| HH\t| MM\t| SS\t| Latitude\t| Longitude\t| Mag\t|Depth\n");
    for (int i = 0; i < length; ++i) {
        printReading(&arr[i]);
    }
    printf("\n");
}