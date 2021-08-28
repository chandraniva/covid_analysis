#ifndef _COV_READ_DATA
#define _COV_READ_DATA

struct daydata {
	int I;
	int R;
	struct daydata* next;
} ;

struct daydata* readdata();

#endif
