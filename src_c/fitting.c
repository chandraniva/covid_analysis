#include<stdio.h>
#include "model.h"
#include "readdata.h"


int main () {
	struct daydata* today;
	today = readdata();
	// now all the day data is in this link list
	int count = 0;
	while(today) {
		printf("%d\t%d\n",count,today->I);
		count++;
		today=today->next;
	}
}

