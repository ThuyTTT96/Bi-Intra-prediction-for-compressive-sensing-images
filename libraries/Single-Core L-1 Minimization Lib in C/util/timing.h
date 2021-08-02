// Functions that approximate MATLAB's handy tic and toc timing functions.
#ifndef __TIMING_H__
#define __TIMING_H__

#ifdef __APPLE__
#define USE_GETTIMEOFDAY
#include <sys/time.h>
#elif __linux
#include <time.h>
#define USE_CLOCK_GETTIME
#else
#error "Platform not recognized!  Not sure what clock to use!"
#endif

#ifdef USE_GETTIMEOFDAY

#define TIMINGINIT struct timeval timeValue1, timeValue2;\
			struct timezone timeZone;\
			double timeDifference;

#define TIC gettimeofday(&timeValue1, &timeZone);

#define TOC gettimeofday(&timeValue2, &timeZone);\
			timeDifference = (double)(timeValue2.tv_sec-timeValue1.tv_sec) + \
					((double)(timeValue2.tv_usec - timeValue1.tv_usec))*1e-6;
#endif


#ifdef USE_CLOCK_GETTIME
#define TIMINGINIT struct timespec timeValue1, timeValue2;\
			double timeDifference;

#define TIC clock_gettime(CLOCK_REALTIME,&timeValue1);

#define TOC clock_gettime(CLOCK_REALTIME,&timeValue2);\
			timeDifference = (double)(timeValue2.tv_sec-timeValue1.tv_sec) + \
					((double)(timeValue2.tv_nsec - timeValue1.tv_nsec))*1e-9;
#endif

#endif
