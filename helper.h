#ifndef _HELPER_H_
#define _HELPER_H_

#include <stdio.h> 

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
void printReading(struct Sensor* reading);

#endif