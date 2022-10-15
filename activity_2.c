#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>
#include "helper.h"
#include "activity_2.h"

#define SIZE 10


void getGrid(int *x, int *y);
void generateReading2(struct Sensor* reading);
void printReading2(struct Sensor* reading);
void printQueue(struct Sensor arr[], int length);

int shutdown = 0;

/**
 * Started on a thread by the Base
 */
void* startBalloon(void* pArg) {
    struct Sensor* sharedReadings = pArg;
    while (shutdown == 0) {
        struct Sensor newReading;
        generateReading2(&newReading);
        enqueue(sharedReadings, SIZE, newReading);

        // printReading2(&newReading);
        printQueue(sharedReadings, SIZE);
        sleep(5);
    }
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
}

void generateReading2(struct Sensor* reading)
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

void printReading2(struct Sensor* reading)
{
    printf(" %d\t| %d\t| %d\t| %d\t| %d\t| %d\t| %.2f\t| %.2f\t| %.2f\t| %.2f\n",
           reading->year, reading->month, reading->day, reading->hour,
           reading->minute, reading->second, reading->lat, reading->lon,
           reading->mag, reading->depth);
}

void printQueue(struct Sensor arr[], int length)
{
    printf(" YYYY\t| MM\t| DD\t| HH\t| MM\t| SS\t| Latitude\t| Longitude\t| Mag\t|Depth\n");
    for (int i = 0; i < length; ++i) {
        printReading2(&arr[i]);
    }
    printf("\n");
}

void getGrid(int *x, int *y) {
    printf("Enter the Grid size");
    printf("\n\tx: ");
    scanf("%d", x);
    printf("\ty: ");
    scanf("%d", y);
}


