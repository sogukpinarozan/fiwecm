/**
 * \Created on: Mar 18, 2020
 * \Author: Asena Durukan, Asli Altiparmak, Elif Ozbay, Hasan Ozan Sogukpinar, Nuri Furkan Pala
 * \file ecm.cu
 * \brief Implementation of ecm.h library.
 */
//TODO: change all function names to fiweX
#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "montgomery.h"
#include "mplib.h"
#include "ecm.h"

__managed__ int primes[21] = {5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83};
__managed__  int prime_powers[5] = {1,2,3,4,5};


void generate_B_smooth(ui z, ui_t l) {
    int i, j, k, prime, power, exceed = 0;
    ui_t z_[2 * l], factor[l];

    z_[0] = 1L;
    for(i = 1; i < l; i++) {
        z_[i] = 0L;
        factor[i] = 0L;
    }
    while(exceed != 1) {
        big_cpy(z, z_, 0, l);
        prime = primes[rand() % PRIMESL];
        power = prime_powers[rand() % POWERSL];
        factor[0] = 1L;
        for(j = 0; j < power; j++) {
            factor[0] *= prime;
        }
        big_mul(z_, z, l, factor, l);
        for(k = l; k < 2 * l; k++) {
            if(z_[k] != 0L) {
                k = 2 * l;
                exceed = 1;
            }
        }
    }
}

__device__ void generateBSmooth(ui z, ui_t l) {

	int t_index = threadIdx.x + (blockIdx.x * blockDim.x);

	if(t_index < SIZE){

		int i, j, k, prime, power, exceed = 0, primeIndex, powerIndex;

		ui z_ = new ui_t[2 * l];
		ui factor = new ui_t[l];

		z_[0] = 1L;
		for(i = 1; i < l; i++) {
			z_[i] = 0L;
			factor[i] = 0L;
		}

		curandState globalState;
		curand_init((int)clock() - t_index, t_index, 0, &globalState);


		while(exceed != 1) {
			bigCopy(z, 0, l, z_, 0);

			primeIndex = gpuGenerate(&globalState) * PRIMESL;
			powerIndex = gpuGenerate(&globalState) * POWERSL;
			prime = primes[primeIndex];
			power = prime_powers[powerIndex];

			factor[0] = 1L;
			for(j = 0; j < power; j++) {
				factor[0] *= prime;
			}

			bigMul(z_, z, l, factor, l);

			for(k = l; k < 2 * l; k++) {
				if(z_[k] != 0L) {
					k = 2 * l;
					exceed = 1;
				}
			}
		}

	}
}

__global__ void fiwEcm(ui d, ui n, ui_t nl, ui mu, ui_t mul, ui returnArr){

	int t_index = threadIdx.x + (blockIdx.x * blockDim.x);

	if(t_index < SIZE){

		MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
		c->A = (ui) malloc(sizeof(ui_t) * nl);
		c->B = (ui) malloc(sizeof(ui_t) * nl);

		PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
		p->X = (ui) malloc(sizeof(ui_t) * nl);
		p->Y = (ui) malloc(sizeof(ui_t) * nl);
		p->Z = (ui) malloc(sizeof(ui_t) * nl);

		PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
		p1->X = (ui) malloc(sizeof(ui_t) * nl);
		p1->Y = (ui) malloc(sizeof(ui_t) * nl);
		p1->Z = (ui) malloc(sizeof(ui_t) * nl);

		ui A24 = new ui_t[nl];

		int i, j, kl, flag = 0;
		ui_t is_zero = 0, is_n = 0, is_one = 0;

		kl = (nl >> 1) + 1; // For the choice of appropriate B. May change.
		ui k = new ui_t[kl];

		//TODO: if n is prime, barret loops infinitely

		PRO_POINT R0, R0_, R1, R1_;

		R0 = (PRO_POINT)malloc(sizeof(PRO_POINT_t));
		R0_ = (PRO_POINT)malloc(sizeof(PRO_POINT_t));
		R1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t));
		R1_ = (PRO_POINT)malloc(sizeof(PRO_POINT_t));

		R0->X = (ui)malloc(sizeof(ui_t) * nl);
		R0->Z = (ui)malloc(sizeof(ui_t) * nl);

		R0_->X = (ui)malloc(sizeof(ui_t) * nl);
		R0_->Z = (ui)malloc(sizeof(ui_t) * nl);

		R1->X = (ui)malloc(sizeof(ui_t) * nl);
		R1->Z = (ui)malloc(sizeof(ui_t) * nl);

		R1_->X = (ui)malloc(sizeof(ui_t) * nl);
		R1_->Z = (ui)malloc(sizeof(ui_t) * nl);

		for(i = 0; i < CRV_THRESHOLD; i++) {

			do {
				//proCurvePoint(d, c, p1, n, nl, mu, mul, &flag);
			} while(flag == -1);
			if(flag == 0) {
				returnArr[t_index] = 1;
				return;
			} else {
				bigGetA24(A24, c->A, n, nl, mu, mul, &flag);
				if(flag == 1) {
					generateBSmooth(k, kl);
					j = 0;
					while(j < k_THRESHOLD) {

						proLadder(p, p1, R0, R0_, R1, R1_, A24, k, kl, n, nl, mu, mul);
						bigIsEqualUi(&is_zero, p->Z, nl, 0L);
						printf("%d: is_zero :%u\n", j, is_zero);

						if(is_zero == 1) {
							binaryGCD(d, k, nl, n, nl);
							bigIsEqualUi(&is_one, d, nl, 1L);
							bigIsEqual(&is_n, d, n, nl);
							if(!is_one && !is_n) {
								returnArr[t_index] = 2;
								return;
							}
						}
						generateBSmooth(k, kl);
						j++;
					}
				} else {
					returnArr[t_index] = 3;
					return;
				}
			}
		}

		printf("%d\n", i);
		returnArr[t_index] = 0;
		return;
	}
}

int ecm(ui d, ui n, ui_t nl) {
    ui A24 = (ui)malloc(sizeof(ui_t) * nl);

    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);

    ui_t mu[nl + 1], B[nl];
    int i, j, coeff, is_one, is_n, flag;

    big_get_mu(mu, n, nl);
    for(i = 0; i <nl; i++) {
        d[i] = 0L;
    }
    for(i = 0; i < CRV_THRESHOLD; i++) {
        do {
            //pro_curve_point(d, c, p1, n, nl, mu, nl + 1, &flag);
        } while(flag == -1);
        if(flag == 0) {
            return 1;
        } else {
            big_get_A24(A24, c->A, n, nl, mu, nl + 1, &flag);
            if(flag == 1) {
                coeff = 1;
                for(j = 1; j < 11; j++) {
                    coeff *= j;
                }
                B[0] = (ui_t)coeff;
                for(j = 1; j < nl; j++) {
                    B[j] = 0L;
                }
                // TODO: How to choose B?
                // TODO: How big is B?
                for(j = 0; j < k_THRESHOLD; j++) {
                    pro_ladder(p, p1, A24, B, nl, n, nl, mu, nl + 1);
                    if(p->Z == 0) {
                        big_gcd(d, nl, B, nl, n, nl);
                        big_is_equal_ui(&is_one, d, nl, 1L);
                        big_is_equal(&is_n, d, n, nl);
                        if(!is_one && !is_n) {
                            return 2;
                            // FIXME: Does not find in ladder!
                        }
                        coeff++;
                        // TODO: How to increment B?
                    }
                }
            } else {
                return 3;
                // FIXME: Does not find in (A+2)/4!
            }
        }
    }
    return 0;
}
