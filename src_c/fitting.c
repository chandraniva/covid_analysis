#include<stdio.h>
#include<nlopt.h>
#include "model.h"
#include "readdata.h"

struct daydata* today;


double opt[] = { 0.0498588,0.0123458,0.00178522,1.73658 };

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
	// FIT
	return square_err(object_sdate,object_edate,&params,today);
}

void test_current_params(int start, int end, char* file_name) {
	today = readdata();
	PARAMS params;
	// TODO make func
	params.beta=x[0];
	params.gamma=x[1];
	params.phi=x[2];
	params.pw=x[3];
	// first wave 76th day to 377 day
	printf("Fitting error: %g\n",square_err_verbose(start,end,&params,today,file_name));
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
	printf("nlopt exit status(+ve is success): %d\n",nlopt_optimize(opt, x, &minf));
	printf("fitted parameters (beta,gama,phi,pw): %g,%g,%g,%g\n",x[0],x[1],x[2],x[3]);
	printf("optimization error: %g\n",minf);
	return x;
}


int main () {
	// second wave 378th day to 576 day
	// test_current_params(76, 377);
	// test_current_params(378, 576);
	
	//430-510
	x = optimize(76,377);
	test_current_params(76,377,"data-1-fit-on-1.dat");
	test_current_params(430,576,"data-2-fit-on-1.dat");
	x = optimize(430,510);
	test_current_params(430,576,"data-2-fit-on-2.dat");
}
