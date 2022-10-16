#include <stdlib.h>
#include "helper.h"

int queue_head = 0;

void enqueue(struct Sensor arr[], int length, struct Sensor element)
{
    if (queue_head == length) {
        for (int i = 1; i < length; ++i) {
            arr[i-1] = arr[i];
        }
        queue_head--;
    }

    arr[queue_head] = element;

    queue_head++;
}

/**
 * Generate a random floating point number from min to max
 * @param min
 * @param max
 * @return
 */
float float_rand(float min, float max)
{
    float range = (max - min);
    float div = RAND_MAX / range;
    return min + (rand() / div);
}

void printReading(struct Sensor* reading)
{
  printf("%d\t| %d\t| %d\t| %d\t| %d\t| %d\t| %.2f\t| %.2f\t| %.2f\t| %.2f\n",
          reading->year, reading->month, reading->day, reading->hour,
          reading->minute, reading->second, reading->lat, reading->lon,
          reading->mag, reading->depth);
}