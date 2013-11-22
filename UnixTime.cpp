#include "time.h"

double getUnixTime()
{
    /* Function to get accurate time, prescission to about some microseconds.
     *
     * Important: to get the function clock_gettime() to work link to library
     * -lrt when linking.
     * */
    struct timespec tv;

    // Get time:
    if(clock_gettime(CLOCK_REALTIME, &tv) != 0) return 0;

    // Returns time in seconds:
    return (((double) tv.tv_sec) + (double) (tv.tv_nsec / 1000000000.0));
}
