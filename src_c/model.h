#ifndef _COV_MODEL
#define _COV_MODEL

typedef struct {
	float S;
	float I;
	float R;
	float Q;
} Y;

typedef struct {
	float beta;
	float gamma;
	float phi;
} PARAMS;

typedef float TIME;

Y integrator(TIME t0, TIME t,  Y y, PARAMS *params);

#endif
