#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "helper.c"

#define SIZE 3

void delay(int seconds);
void getGrid(int *x, int *y);
void generate(struct Sensor* reading);
void printReading(struct Sensor* reading);
void printQueue(struct Sensor arr[], int length);
float float_rand(float min, float max);

int main() {
    struct Sensor readings[SIZE];
    struct Sensor newReading;
    generate(&newReading);
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

    reading->lat = float_rand(-20, -10);
    reading->lon = float_rand(150, 170);
    reading->mag = float_rand(1, 5);
    reading->depth = float_rand(4, 10);
}

void printReading(struct Sensor* reading)
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
        printReading(&arr[i]);
    }
    printf("\n");
}

/**
 * Adapted from answer:
 * https://stackoverflow.com/questions/13408990/how-to-generate-random-float-number-in-c
 * @param min
 * @param max
 * @return
 */
float float_rand(float min, float max)
{
    float scale = rand() / (float) RAND_MAX;
    return min + scale * ( max - min );
}

void getGrid(int *x, int *y) {
    printf("Enter the Grid size");
    printf("\n\tx: ");
    scanf("%d", x);
    printf("\ty: ");
    scanf("%d", y);
}

/**
 * Delay the program
 * @param seconds Time to delay
 */
void delay(int seconds) {
    int ms = 1000 * seconds;
    clock_t start_time = clock();
    // looping till required time is not achieved
    while (clock() < start_time + ms);
}

