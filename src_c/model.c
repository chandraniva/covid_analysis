#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include "model.h"
#include "readdata.h"

// this is clone of model.py in c
// refer model.py
// need to optimize with data

//globals

double N = 1.3526e9;

double fear_function(TIME t, Y y, PARAMS *params) {
	return params->phi * pow(y.I,2);
	//return 0;
}

Y dydt_fear(TIME t, Y y, PARAMS *params) {
	Y dydt;
	double fear_factor = fear_function(t,y,params);
	dydt.I = params->beta * y.I * y.S / N - params->gamma * y.I ;
	dydt.R = params->gamma * y.I;
	dydt.S = - params->beta * y.I * y.S / N - y.S * fear_factor / N ;
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
double square_err(TIME t0, TIME t, PARAMS *params , struct daydata* today) {
	Y y;

	double serror = 0;
	// y has the initial values in the start
	double h = 0.1;
	// t0 goes to t

	Y delta;
	Y k[4];

	int i;
	double time;

	for(i=0; i<76; i++) {
		today = today->next;
	}
	y.I = today->I;
	y.R = today->R;
	y.Q = 0;
	y.S = N - y.I - y.R - y.Q;

	while( t0<t && today ) {
		time = t0;
		for(i=0;i<10;i++) {
			time+=h;
			k[0] = y_sprod(h,dydt_fear(time, y, params));
			k[1] = y_sprod(h,dydt_fear(time+h/2, y_sum(y,y_sprod(0.5,k[0])), params));
			k[2] = y_sprod(h,dydt_fear(time+h/2, y_sum(y,y_sprod(0.5,k[1])), params));
			k[3] = y_sprod(h,dydt_fear(time+h, y_sum(y,k[2]), params));
			//printf("k.I : %g %g %g %g\n",k[0].I,k[1].I,k[2].I,k[3].I);

			delta = y_sprod(1.0/6.0,y_sum(y_sum(k[0],y_sprod(2,k[1])),y_sum(y_sprod(2,k[2]),k[3])));
			//printf("delta I : %g\n",delta.I);
			y = y_sum(y,delta);
		}
		//printf("%d\t%g\n",t0,y.I);
		printf("%d\t%g\t%d\n",t0,y.I,today->I);

		serror += pow((y.I - today->I),2);
		t0++;
		today = today->next;
	}
	printf("%g\n",y.I);
	printf("error: %g\n",serror);
	return serror;
}


