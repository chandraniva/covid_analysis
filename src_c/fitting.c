#include<stdio.h>
#include<nlopt.h>
#include "model.h"
#include "readdata.h"

struct daydata* today;


double x[] = { 
	0.0486891,0.0119584,5.18946e-05
	};

void print_actual_data() {
	struct daydata* today;
	today = readdata();
	// now all the day data is in this link list
	int count = 0;
	while(today) {
		printf("%d\t%d\t%d\n",count,today->I,today->R);
		count++;
		today=today->next;
	}
}

double objective_fun_wrapper(unsigned n, const double* x, double* grad, void* f_data) {
	PARAMS params;
	params.beta=x[0];
	params.gamma=x[1];
	params.phi=x[2];
	// first wave 76th day to 377 day
	return square_err(76,377,&params,today);
}

void test_current_params_1st() {
	today = readdata();
	PARAMS params;
	params.beta=x[0];
	params.gamma=x[1];
	params.phi=x[2];
	// first wave 76th day to 377 day
	square_err(76,377,&params,today);
}

void test_current_params_2nd() {
	today = readdata();
	PARAMS params;
	params.beta=x[0];
	params.gamma=x[1];
	params.phi=x[2];
	// first wave 76th day to 377 day
	square_err(378,576,&params,today);
}


void optimize() {
	double lb[3] = { 0 , 0 , 0};
	today = readdata();

	nlopt_opt opt;

	opt = nlopt_create(
			NLOPT_LN_BOBYQA
			, 3);
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective (opt, objective_fun_wrapper, NULL);
	nlopt_set_xtol_rel(opt, 1e-4);

	double minf;
	printf("%d\n",nlopt_optimize(opt, x, &minf));
	printf("%g,%g,%g\n",x[0],x[1],x[2]);
}


int main () {
	//optimize();
	test_current_params_2nd();
	//print_actual_data();
}
