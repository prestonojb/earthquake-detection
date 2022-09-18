#include <stdio.h>
#include <stdlib.h>
#include <time.h>

float float_rand(float min, float max) ;
void getCurrentDatetime(int* currentDatetime);
void getReadings(float* readings);

int main() {
  srand (time(NULL));

  int* currentDatetime = (int*)malloc(6 * sizeof(int));
  float* readings = (float*)malloc(4 * sizeof(float)); // Heap memory

  getCurrentDatetime(currentDatetime);
  getReadings(readings); 

  // for (int i=0; i<6; i++) {
  //   printf("currentDatetime[%d]: %d \n", i, currentDatetime[i]);
  // }

  // for (int i=0; i<4; i++) {
  //   printf("readings[%d]: %.2f \n", i, readings[i]);
  // }

  return 0;
}

// currentDatetime = [YYYY, MM, DD, HH, MM, SS]
void getCurrentDatetime(int* currentDatetime) {
  // Generate datetime (reading #1-6)
  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  int datetimeArr[6] = {tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec};

  for(int i =0; i<6; i++) {
    currentDatetime[i] = datetimeArr[i];
  }
}

// Readings = [latitude, longitude, magnitude, depth]
void getReadings(float* readings) {

  float_rand(0, 0); // Workaround for same first random number
  float latitude = float_rand(-90, 90);
  float longitude = float_rand(-180, 180);
  float magnitude = float_rand(0, 9);
  float depth = float_rand(0, 700);

  float readingsArr[4] = {latitude, longitude, magnitude, depth};

  for(int i=0; i<4; i++) {
    readings[i] = readingsArr[i];
  }
}

/* generate a random floating point number from min to max */
float float_rand(float min, float max) 
{
    float range = (max - min); 
    float div = RAND_MAX / range;
    return min + (rand() / div);
}