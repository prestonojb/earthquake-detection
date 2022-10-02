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

int queue_head = 0;
void enqueue(struct Sensor arr[], int length, struct Sensor element);

void enqueue(struct Sensor arr[], int length, struct Sensor element)
{
    queue_head++;
    if (queue_head == length) {
        for (int i = 1; i < length; ++i) {
            arr[i-1] = arr[i];
        }
        queue_head--;
    }

    arr[queue_head] = element;
}
