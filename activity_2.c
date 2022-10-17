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

void printQueue(struct Sensor arr[], int length)
{
    printf(" YYYY\t| MM\t| DD\t| HH\t| MM\t| SS\t| Latitude\t| Longitude\t| Mag\t|Depth\n");
    for (int i = 0; i < length; ++i) {
        printReading(&arr[i]);
    }
    printf("\n");
}