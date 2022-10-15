#ifndef _ACTIVITY_3_H_
#define _ACTIVITY_3_H_

void defineSensorType(MPI_Datatype* SensorType);
void defineDataLogType(MPI_Datatype* DataLogType, MPI_Datatype SensorType);
#endif