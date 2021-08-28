#include<stdio.h>
#include<math.h>

// this is clone of model.py in c
// refer model.py
// need to optimize with data

//globals

float N = 1.3526e9;

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

float fear_function(TIME t, Y y, PARAMS *params) {
	return params->phi * sqrt(y.I);
}

Y dydt_fear(TIME t, Y y, PARAMS *params) {
	Y dydt;
	float betaISN = params->beta * y.I * y.S / N;
	float gamaI = params->gamma * y.I;
	float fear_factor = fear_function(t,y,params);
	dydt.I = betaISN - gamaI ;
	dydt.R = gamaI;
	dydt.S = - betaISN - y.S * fear_factor / N ;
	dydt.Q = y.S * fear_factor / N ;
	return dydt;
}

Y y_sprod(float h, Y y) {
	y.I *= h; y.Q *= h; y.S *= h; y.R *= h;
	return y;
}

Y y_sum(Y y1, Y y2) {
	y1.I += y2.I; y1.Q += y2.Q; y1.S += y2.S; y1.R += y2.R;
	return y1;
}

Y integrator(TIME t0, TIME t,  Y y, PARAMS *params) {
	// y has the initial values in the start
	float h = 1e-1;
	// t0 goes to t

	Y delta;
	Y k[4];

	int i;
	while( t0<t ) {

		k[0] = y_sprod(h,dydt_fear(t0, y, params));
		k[1] = y_sprod(h,dydt_fear(t0+h/2, y_sum(y,y_sprod(0.5,k[1])), params));
		k[2] = y_sprod(h,dydt_fear(t0+h/2, y_sum(y,y_sprod(0.5,k[2])), params));
		k[3] = y_sprod(h,dydt_fear(t0+h, y_sum(y,k[3]), params));

		delta = y_sprod(1.0/6.0,y_sum(y_sum(k[0],y_sprod(2,k[1])),y_sum(y_sprod(2,k[2]),k[3])));
		y = y_sum(y,delta);
		printf("%f\t%f\n",t0,y.I);
		t0+=h;
	}
	return y;
}


int main() {
	PARAMS params = {
		.beta = 0.08,
		.gamma = 0.02,
		.phi = 0.002
	};
	Y y = {
		.S = 1e10,
		.I = 1 ,
		.R = 3 ,
		.Q = 0 
	};
	y = integrator(1,300,y,&params);
}
