#include<stdio.h>
#include<nlopt.h>
#include "model.h"
#include "readdata.h"

struct daydata* today;


double opt[] = { 
			0.0495045,0.0123377,0.00050104,1.83065
	};

double *x = opt;
TIME object_edate;
TIME object_sdate;

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
	params.pw=x[3];
	// first wave 76th day to 377 day
	// FIT
	return square_err(object_sdate,object_edate,&params,today);
}

void test_current_params(int start, int end) {
	today = readdata();
	PARAMS params;
	// TODO make func
	params.beta=x[0];
	params.gamma=x[1];
	params.phi=x[2];
	params.pw=x[3];
	// first wave 76th day to 377 day
	printf("%g\n",square_err_verbose(start,end,&params,today));
}



double* optimize(TIME sdate, TIME edate) {

	object_edate = edate;
	object_sdate = sdate;

	double lb[] = { 0 , 0.0123 , 0 , 0};
	today = readdata();

	nlopt_opt opt;

	opt = nlopt_create( NLOPT_LN_BOBYQA , sizeof(lb)/sizeof(lb[0]) );
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_min_objective (opt, objective_fun_wrapper, NULL);
	nlopt_set_xtol_rel(opt, 1e-10);

	double minf;
	printf("%d\n",nlopt_optimize(opt, x, &minf));
	printf("%g,%g,%g,%g\n",x[0],x[1],x[2],x[3]);
	printf("optimization error: %g\n",minf);
	return x;
}


int main () {
	// first wave 378th day to 576 day
	//test_current_params(76, 377);
	//test_current_params(378, 576);
	x = optimize(420,500);
	test_current_params(420, 500);
}
