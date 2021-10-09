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
	return params->phi * pow(y.I,params->pw);
		//* (int)(1/(1 + exp(-1e308 * (y.I - 5*1e5) ) ));
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

Y rk_increment_day(TIME t, Y y,PARAMS *params){
	Y k[4];
	double time = t;
	Y delta;
	int i;

	double h=0.1;
	for(i=0;i<10;i++) {
		time+=h;
		k[0] = y_sprod(h,dydt_fear(time, y, params));
		k[1] = y_sprod(h,dydt_fear(time+h/2, y_sum(y,y_sprod(0.5,k[0])), params));
		k[2] = y_sprod(h,dydt_fear(time+h/2, y_sum(y,y_sprod(0.5,k[1])), params));
		k[3] = y_sprod(h,dydt_fear(time+h, y_sum(y,k[2]), params));

		delta = y_sprod(1.0/6.0,y_sum(y_sum(k[0],y_sprod(2,k[1])),y_sum(y_sprod(2,k[2]),k[3])));
		y = y_sum(y,delta);
	}
	return y;
}
void simulate(Y y, int total_days, PARAMS *params) {
	for(int day=0; day<total_days ; day++) {
		y = rk_increment_day(day,y,params);
		printf("%d\t%g\n",day,y.I);
	}
}


// STR LINE 
//
double square_err_st_line(TIME t0, TIME t, PARAMS *params , struct daydata* today) {

	FILE *oftr;
	oftr = fopen("c.dat","w");

	double serror = 0;
	// y has the initial values in the start
	int i;
	for(i=0; i<t0; i++) {
		today = today->next;
	}

	// x day
	double temp;

	while( t0<t && today ) {
		temp = - params->beta * exp(-t0);
		serror += pow( temp - today->I  ,2);
		fprintf(oftr,"%d\t%g\t%d\n",t0,temp,today->I);
		t0++;
		today = today->next;
	}
	fclose(oftr);
	return serror;
}

//
double square_err(TIME t0, TIME t, PARAMS *params , struct daydata* today) {

	Y y;
	double serror = 0;
	// y has the initial values in the start
	int i;
	for(i=0; i<t0; i++) {
		today = today->next;
	}

	// initial values for the model
	y.I = today->I;
	y.R = today->R;
	y.Q = 0;
	y.S = N - y.I - y.R - y.Q;

	while( t0<t && today ) {
		y = rk_increment_day(t0,y,params);
		serror += pow((y.I - today->I),2);
		t0++;
		today = today->next;
	}
	return serror;
}


double square_err_verbose(TIME t0, TIME t, PARAMS *params , struct daydata* today, char* file_name) {

	FILE *oftr;
	oftr = fopen(file_name,"w");
	Y y;

	double serror = 0;
	// y has the initial values in the start
	int i;
	for(i=0; i<t0; i++) {
		today = today->next;
	}

	// initial values for the model
	y.I = today->I;
	y.R = today->R;
	y.Q = 0;
	y.S = N - y.I - y.R - y.Q;

	while( t0<t && today ) {
		y = rk_increment_day(t0,y,params);
		fprintf(oftr,"%d\t%g\t%d\n",t0,y.I,today->I);
		serror += pow((y.I - today->I),2);
		t0++;
		today = today->next;
	}
	//printf("%g\n",y.I);
	//printf("error: %g\n",serror);
	fclose(oftr);
	return serror;
}


