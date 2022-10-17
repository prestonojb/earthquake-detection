#ifndef _HELPER_H_
#define _HELPER_H_

#include <stdio.h>

#define pi 3.14159265358979323846

#define DEFAULT_MAGNITUDE_UPPER_THRESHOLD 2.5
#define DEFAULT_DIFF_IN_DISTANCE_THRESHOLD_IN_KM 20
#define DEFAULT_DIFF_IN_MAGNITUDE_THRESHOLD 1

#define QUEUE_SIZE 10

#define BASE_STATION 0
#define TERMINATION_TAG 10

struct Sensor {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second;
    float lat;
    float lon;
    float mag;
    float depth;
};

struct DataLog {
    int reporterRank;
    struct Sensor reporterData;

    int topRank;
    struct Sensor topData;
    int bottomRank;
    struct Sensor bottomData;
    int leftRank;
    struct Sensor leftData;
    int rightRank;
    struct Sensor rightData;
};

extern int queue_head;
void enqueue(struct Sensor arr[], int length, struct Sensor element);
float float_rand(float min, float max);
float float_rand_seed(float min, float max, int seed);
void printReading(struct Sensor* reading);

#endif