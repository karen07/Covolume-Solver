// Exemple_solution.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdlib.h>
#include "include/Covolume_solver.h"

#pragma comment(lib, "Covolume_solver")

int main()
{
	double *d_arr;
	double *v_arr;
	double *p_arr;
	double *x;

	FILE *file = NULL;
	int i;

	riman_solver_covolume(100.0, 0.0, 100.0e+6, 1.0, 0.0, 0.1e+6, 0.001, 1.3, 1, 0.4, 0.0002, 1000, 10e-5, &d_arr, &v_arr, &p_arr, &x);

	fopen_s(&file, "C:/Users/Karen/Desktop/karen.dat", "w+");
	if (file == NULL) {
		return 0;
	}
	fprintf(file, "#X              #D            #V              #P\n");
	for (i = 0; i < 1000; i++) {
		fprintf(file, "%f %f %f %f\n", x[i], d_arr[i], v_arr[i], p_arr[i]);
	}
	printf("Success\n");
	system("pause");
	return 0;
}

