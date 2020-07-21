/**
 * \Created on: Mar 18, 2020
 * \Author: Asena Durukan, Asli Altiparmak, Elif Ozbay, Hasan Ozan Sogukpinar, Nuri Furkan Pala
 * \file main.cu
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "test.h"
#include "ecm.h"


__device__ void timeTest1(int *a){
	int t_index = threadIdx.x + (blockIdx.x * blockDim.x);

	if (t_index < SIZE) {
		*a +=5;
	}

}
__global__ void timeTest() {

	int t_index = threadIdx.x + (blockIdx.x * blockDim.x);

	if (t_index < SIZE) {

		int a = 0;

		for(int i = 0; i < 10000000; i++){
			timeTest1(&a);
		}

	}
}
int main() {
//     time_t start, end;
//     start = time(NULL);
	   pro_curve_point_gmp_test(1000);
//     aff_curve_point_gmp_test();
//     pro_add_gmp_test();
//     pro_add_magma_test();
//     pro_dbl_magma_test();
//     pro_ladder_magma_test();
//     pro_ladder_gmp_test(1000);
//		ecm_gmp_test(1000);
//    end = time(NULL);
//    printf("Time spent: %ld\n", end - start);

//	trace();

	return 0;
}

