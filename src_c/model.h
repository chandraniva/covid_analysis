#ifndef _COV_MODEL
#define _COV_MODEL
#include "readdata.h"

typedef struct {
	double S;
	double I;
	double R;
	double Q;
} Y;

typedef struct {
	float beta;
	float gamma;
	float phi;
} PARAMS;

typedef int TIME;

double square_err(TIME t0, TIME t, PARAMS *params, struct daydata*);

#endif
