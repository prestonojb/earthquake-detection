compile: activity_1.c activity_2.c activity_3.c helper.c
	mpicc activity_1.c activity_2.c activity_3.c helper.c -lm -o earthquake_detection.o
clean:
	rm earthquake_detection.o