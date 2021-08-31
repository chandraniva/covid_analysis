#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "readdata.h"

struct daydata* readdata() {
	int a; //temp variable

	FILE* f1;
	f1 = fopen("case_time_series.csv","r");

	char buffer[120];
	char *key;

	struct daydata* HEAD = NULL;
	struct daydata* current = NULL;
	
	int line=0;
	while(fgets(buffer,120,f1)){

		line++;
		if ( !(line > 1) ) {
			continue;
		}

		if ( ! HEAD )  {
			HEAD = (struct daydata*) malloc(sizeof(int));
			current = HEAD;
		} else {
			current->next = (struct daydata*) malloc(sizeof(int));
			current = current->next;
		}
		current->next = NULL;
		strtok(buffer,",");
		strtok(NULL,",");
		strtok(NULL,",");
		key = strtok(NULL,",");
		//printf("%ld\t",atol(key));
		a = atol(key);
		strtok(NULL,",");
		key = strtok(NULL,",");
		//printf("%ld\n",atol(key));
		current->R = atol(key);
		current->I = a - current->R ;
	}
	fclose(f1);

	return HEAD;
}



