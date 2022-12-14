#ifndef _ACTIVITY_1_H_
#define _ACTIVITY_1_H_

int init_nodes(int m, int n, float magnitude_upper_threshold, float diff_in_distance_threshold_in_km, float diff_in_magnitude_threshold, MPI_Comm world_comm, MPI_Comm comm);
bool areMatchingReadings(struct Sensor* readingA, struct Sensor* readingB);
bool areMatchingMagnitudes(float magnitudeA, float magnitudeB);
float distance(float lat1, float lon1, float lat2, float lon2);
#endif