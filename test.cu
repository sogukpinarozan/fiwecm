/**
 * \file test.c
 * \brief Implementation of test.h library.
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mplib.h"
#include "montgomery.h"
#include "ecm.h"
#include "test.h"
#include <cudaProfiler.h>


void pro_curve_point_gmp_test(int THRESHOLD) {
    mpz_t mp_n, mp_a, mp_b, mp_x, mp_y, mp_z, mp_y2, mp_x2, mp_x3, mp_by2, mp_by2z, mp_ax2, mp_ax2z, mp_z2, mp_xz2, mp_d, mp_mod, mp_dif;
    ui_t nl;


    int i, j, res = 0, trues = 0, falses = 0, singulars = 0, *d_flag;
    int flag[SIZE];

    mpz_init(mp_n);
    mpz_init(mp_a);
    mpz_init(mp_b);
    mpz_init(mp_x);
    mpz_init(mp_y);
    mpz_init(mp_z);
    mpz_init(mp_y2);
    mpz_init(mp_x2);
    mpz_init(mp_x3);
    mpz_init(mp_by2);
    mpz_init(mp_by2z);
    mpz_init(mp_ax2);
    mpz_init(mp_ax2z);
    mpz_init(mp_z2);
    mpz_init(mp_xz2);
    mpz_init(mp_d);
    mpz_init(mp_mod);
    mpz_init(mp_dif);
    mpz_set_ui(mp_mod, 0L);
    mpz_set_ui(mp_dif, 0L);

    srand(time(NULL));

    nl = 10;
    ui d_n, d_mu, d_d, d_controlArr, d;
    ui_t n[nl], mu[nl + 1];

    d = (ui)malloc(sizeof(ui_t) * SIZE * nl);

    cudaMalloc(&d_n, nl * sizeof(ui_t));
    cudaMalloc(&d_mu, (nl + 1) * sizeof(ui_t));
    cudaMalloc(&d_flag, SIZE * sizeof(ui_t));
    cudaMalloc(&d_controlArr, SIZE * sizeof(ui_t));
    cudaMalloc(&d_d, SIZE * nl * sizeof(ui_t));

    big_rand(n, nl);
    n[0] -= 37;
    big_print(stdout, n, nl, "n", NULL);
    cudaMemcpy(d_n, n, sizeof(ui_t) * nl, cudaMemcpyHostToDevice);

    big_get_mu(mu, n, nl);

    cudaMemcpy(d_mu, mu, sizeof(ui_t) * (nl + 1), cudaMemcpyHostToDevice);

    size_t pValue;
	cudaError_t error;

	cudaDeviceSetLimit(cudaLimitStackSize, 10000 * sizeof(ui_t));
	cudaDeviceGetLimit(&pValue, cudaLimitStackSize);

	cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, 10000 * sizeof(ui_t));
	cudaDeviceGetLimit(&pValue, cudaLimitDevRuntimePendingLaunchCount);

	cudaDeviceGetLimit(&pValue, cudaLimitMallocHeapSize);
	cudaDeviceSetLimit(cudaLimitMallocHeapSize , pValue * 50);
	cudaDeviceGetLimit(&pValue, cudaLimitMallocHeapSize);

    //TODO: static or dynamic?
    MONTG_CURVE c;
    PRO_POINT p;


	cudaMallocManaged(&c, SIZE * sizeof(MONTG_CURVE_t));
	cudaMallocManaged(&p, SIZE * sizeof(PRO_POINT_t));

    for(i = 0; i < SIZE; i++){
		cudaMallocManaged(&(c[i].A), nl * sizeof(MONTG_CURVE_t));
		cudaMallocManaged(&(c[i].B), nl * sizeof(MONTG_CURVE_t));
		cudaMallocManaged(&(p[i].X), nl * sizeof(MONTG_CURVE_t));
		cudaMallocManaged(&(p[i].Y), nl * sizeof(MONTG_CURVE_t));
		cudaMallocManaged(&(p[i].Z), nl * sizeof(MONTG_CURVE_t));
    }

	proCurvePoint<<<BLOCKNUM, THREADNUM>>>(d_d, c, p, d_n, nl, d_mu, nl + 1, d_flag);
    cudaDeviceSynchronize();

    error = cudaGetLastError();
    if(error != cudaSuccess)
    	fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(error));

    cudaMemcpy(d, d_d, sizeof(ui_t) * SIZE * nl, cudaMemcpyDeviceToHost);

    cudaMemcpy(flag, d_flag, sizeof(ui_t) * SIZE, cudaMemcpyDeviceToHost);

//    int count = 0;
//    for(i = 0; i < SIZE; i++){
//    	//printf("%d: %d\n", i, flag[i]);
//    	if(flag[i] == 1)
//    		count++;
//    }

//    printf("count : %d\n", count);


    for(i = 0; i < SIZE; i++){
		 mpz_import(mp_n, nl, -1, 4, 0, 0, n);
		 // n mod d = 0 check;
		if(flag[i] == 0) {
			mpz_import(mp_d, nl, -1, 4, 0, 0, d);
			mpz_mod(mp_mod, mp_n, mp_d);
			res = mpz_cmp_ui(mp_mod, 0);
			if(res == 0) {
				trues++;
			} else {
				falses++;
				printf("False at index: %d\n", i);
			}
			//(BY2Z - (X3 + AX2Z + XZ2)) mod n eq 0;
		}
		else if(flag[i] == 1){
			mpz_import(mp_a, nl, -1, 4, 0, 0, c[i].A);
			mpz_import(mp_b, nl, -1, 4, 0, 0, c[i].B);
			mpz_import(mp_x, nl, -1, 4, 0, 0, p[i].X);
			mpz_import(mp_y, nl, -1, 4, 0, 0, p[i].Y);
			mpz_import(mp_z, nl, -1, 4, 0, 0, p[i].Z);
			mpz_mul(mp_y2, mp_y, mp_y);
			mpz_mul(mp_by2, mp_b, mp_y2);
			mpz_mul(mp_by2z, mp_by2, mp_z);
			mpz_mul(mp_x2, mp_x, mp_x);
			mpz_mul(mp_x3, mp_x2, mp_x);
			mpz_mul(mp_ax2, mp_a, mp_x2);
			mpz_mul(mp_ax2z, mp_ax2, mp_z);
			mpz_mul(mp_z2, mp_z, mp_z);
			mpz_mul(mp_xz2, mp_x, mp_z2);
			mpz_sub(mp_dif, mp_by2z, mp_x3);
			mpz_sub(mp_dif, mp_dif, mp_ax2z);
			mpz_sub(mp_dif, mp_dif, mp_xz2);
			mpz_mod(mp_dif, mp_dif, mp_n);
			res = mpz_cmp_ui(mp_dif, 0);
			if(res == 0) {
				trues++;
			}
			else{
				falses++;
				printf("False at index: %d\n", i);
			}

		}
		else{
			singulars++;
			printf("Singularity at index: %d\n", i);
		}
    }

    printf("TRUE: %d\n", trues);
    printf("FALSE: %d\n", falses);
    printf("SINGULAR: %d\n", singulars);


    i = 0;
    int t_index = 0;
//
    while(i < SIZE * nl){
    	printf("t_index: %d: ", t_index);
    	t_index++;

    	j = i;

    	while(j < (i + nl)){
    		printf("%u ", d[j]);
    		j++;
    	}

    	printf("\n");
    	i += nl;

    }

	cudaDeviceReset();
}

void aff_curve_point_gmp_test(int THRESHOLD) {
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    AFF_POINT p = (AFF_POINT)malloc(sizeof(AFF_POINT_t) * 1);
    mpz_t mp_n, mp_a, mp_b, mp_x, mp_y, mp_y2, mp_x2, mp_x3, mp_by2, mp_ax2, mp_d, mp_mod, mp_dif;
    ui_t nl;
    int i, j, res = 0, flag = 0, trues = 0, falses = 0, singulars = 0;

    mpz_init(mp_n);
    mpz_init(mp_a);
    mpz_init(mp_b);
    mpz_init(mp_x);
    mpz_init(mp_y);
    mpz_init(mp_y2);
    mpz_init(mp_x2);
    mpz_init(mp_x3);
    mpz_init(mp_by2);
    mpz_init(mp_ax2);
    mpz_init(mp_d);
    mpz_init(mp_mod);
    mpz_init(mp_dif);
    mpz_set_ui(mp_mod, 0L);
    mpz_set_ui(mp_dif, 0L);

    for (i = 0; i < THRESHOLD; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], d[nl];

        flag = 0;
        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        for (j = 0; j < nl; j++) {
            d[j] = 0;
        }

        aff_curve_point(d, c, p, n, nl, mu, nl + 1, &flag);

        mpz_import(mp_n, nl, -1, 4, 0, 0, n);
        if(flag == 0) {
            mpz_import(mp_d, nl, -1, 4, 0, 0, d);
            mpz_mod(mp_mod, mp_n, mp_d);
            res = mpz_cmp_ui(mp_mod, 0);
            if(res == 0) {
                trues++;
            } else {
                falses++;
                printf("False at index: %d\n", i);
            }
        } else if(flag == 1){
            mpz_import(mp_a, nl, -1, 4, 0, 0, c->A);
            mpz_import(mp_b, nl, -1, 4, 0, 0, c->B);
            mpz_import(mp_x, nl, -1, 4, 0, 0, p->x);
            mpz_import(mp_y, nl, -1, 4, 0, 0, p->y);
            mpz_mul(mp_y2, mp_y, mp_y);
            mpz_mul(mp_by2, mp_b, mp_y2);
            mpz_mul(mp_x2, mp_x, mp_x);
            mpz_mul(mp_x3, mp_x2, mp_x);
            mpz_mul(mp_ax2, mp_a, mp_x2);
            mpz_sub(mp_dif, mp_by2, mp_x3);
            mpz_sub(mp_dif, mp_dif, mp_ax2);
            mpz_sub(mp_dif, mp_dif, mp_x);
            mpz_mod(mp_dif, mp_dif, mp_n);
            res = mpz_cmp_ui(mp_dif, 0);
            if(res == 0) {
                trues++;
            } else {
                falses++;
                printf("False at index: %d\n", i);
            }
        } else {
            singulars++;
            printf("Singularity at index: %d\n", i);
        }
    }
    printf("TRUE: %d\n", trues);
    printf("FALSE: %d\n", falses);
    printf("SINGULAR: %d\n", singulars);
}

void pro_add_gmp_test(int THRESHOLD) {
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p2 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT pd = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    mpz_t  mp_n, mp_x1, mp_x2, mp_xd, mp_z1, mp_z2, mp_zd, mp_a, mp_b, mp_c, mp_d, mp_da, mp_cb, mp_e, mp_f, mp_e2, mp_f2, mp_g, mp_h, mp_X, mp_Z;
    ui_t nl;
    int i, trues = 0, falses = 0;

    mpz_init(mp_n);
    mpz_init(mp_x1);
    mpz_init(mp_x2);
    mpz_init(mp_xd);
    mpz_init(mp_z1);
    mpz_init(mp_z2);
    mpz_init(mp_zd);
    mpz_init(mp_a);
    mpz_init(mp_b);
    mpz_init(mp_c);
    mpz_init(mp_d);
    mpz_init(mp_da);
    mpz_init(mp_cb);
    mpz_init(mp_e);
    mpz_init(mp_f);
    mpz_init(mp_e2);
    mpz_init(mp_f2);
    mpz_init(mp_g);
    mpz_init(mp_h);
    mpz_init(mp_X);
    mpz_init(mp_Z);

    for (i = 0; i < THRESHOLD; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], X1[nl], X2[nl], Xd[nl], Z1[nl], Z2[nl], Zd[nl];

        mpz_set_ui(mp_X, 0L);
        mpz_set_ui(mp_Z, 0L);

        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        big_mod_rand(X1, nl, n, nl, mu, nl + 1);
        big_mod_rand(X2, nl, n, nl, mu, nl + 1);
        big_mod_rand(Xd, nl, n, nl, mu, nl + 1);
        big_mod_rand(Z1, nl, n, nl, mu, nl + 1);
        big_mod_rand(Z2, nl, n, nl, mu, nl + 1);
        big_mod_rand(Zd, nl, n, nl, mu, nl + 1);

        mpz_import(mp_n, nl, -1, 4, 0, 0, n);
        mpz_import(mp_x1, nl, -1, 4, 0, 0, X1);
        mpz_import(mp_x2, nl, -1, 4, 0, 0, X2);
        mpz_import(mp_xd, nl, -1, 4, 0, 0, Xd);
        mpz_import(mp_z1, nl, -1, 4, 0, 0, Z1);
        mpz_import(mp_z2, nl, -1, 4, 0, 0, Z2);
        mpz_import(mp_zd, nl, -1, 4, 0, 0, Zd);

        p1->X = X1;
        p2->X = X2;
        pd->X = Xd;
        p1->Z = Z1;
        p2->Z = Z2;
        pd->Z = Zd;

        pro_add(p, p1, p2, pd, n, nl, mu, nl + 1);
        mpz_import(mp_X, nl, -1, 4, 0, 0, p->X);
        mpz_import(mp_Z, nl, -1, 4, 0, 0, p->Z);

        mpz_add(mp_a, mp_x2, mp_z2);
        mpz_mod(mp_a, mp_a, mp_n);
        mpz_sub(mp_b, mp_x2, mp_z2);
        while(mpz_sgn(mp_b) == -1) {
            mpz_add(mp_b, mp_b, mp_n);
        }
        mpz_add(mp_c, mp_x1, mp_z1);
        mpz_mod(mp_c, mp_c, mp_n);
        mpz_sub(mp_d, mp_x1, mp_z1);
        while(mpz_sgn(mp_d) == -1) {
            mpz_add(mp_d, mp_d, mp_n);
        }
        mpz_mul(mp_da, mp_d, mp_a);
        mpz_mod(mp_da, mp_da, mp_n);
        mpz_mul(mp_cb, mp_c, mp_b);
        mpz_mod(mp_cb, mp_cb, mp_n);
        mpz_add(mp_e, mp_da, mp_cb);
        mpz_mod(mp_e, mp_e, mp_n);
        mpz_sub(mp_f, mp_da, mp_cb);
        while(mpz_sgn(mp_f) == -1) {
            mpz_add(mp_f, mp_f, mp_n);
        }
        mpz_mul(mp_e2, mp_e, mp_e);
        mpz_mod(mp_e2, mp_e2, mp_n);
        mpz_mul(mp_f2, mp_f, mp_f);
        mpz_mod(mp_f2, mp_f2, mp_n);
        mpz_mul(mp_g, mp_zd, mp_e2);
        mpz_mod(mp_g, mp_g, mp_n);
        mpz_mul(mp_h, mp_xd, mp_f2);
        mpz_mod(mp_h, mp_h, mp_n);

        int res = mpz_cmp(mp_g, mp_X) + mpz_cmp(mp_h, mp_Z);
        if(res == 0) {
            trues++;
        } else {
            falses++;
        }
    }
    printf("TRUE: %d\n", trues);
    printf("FALSE: %d\n", falses);
}

void pro_add_magma_test(int THRESHOLD) {
    FILE *fp = fopen("/home/ozbayelif/Development/FIWE/ecm/pro_add_test.magma", "a");
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p2 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT pd = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    ui_t nl;
    int i;

    fprintf(fp, "clear;\n");
    fprintf(fp, "/****************************************************************************/\n");
    fprintf(fp, "ADDM:=function(X1, Z1, X2, Z2, Xd, Zd, n)\n");
    fprintf(fp, "return (Zd * (((X1 - Z1) * (X2 + Z2)) + ((X1 + Z1) * (X2 - Z2)))^2) mod n, (Xd * (((X1 - Z1) * (X2 + Z2)) - ((X1 + Z1) * (X2 - Z2)))^2) mod n;\n");
    fprintf(fp, "end function;\n");
    fprintf(fp, "trues := 0;\n");
    fprintf(fp, "falses := 0;\n");

    for (i = 0; i < THRESHOLD; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], X1[nl], X2[nl], Xd[nl], Z1[nl], Z2[nl], Zd[nl];

        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        big_mod_rand(X1, nl, n, nl, mu, nl + 1);
        big_mod_rand(X2, nl, n, nl, mu, nl + 1);
        big_mod_rand(Xd, nl, n, nl, mu, nl + 1);
        big_mod_rand(Z1, nl, n, nl, mu, nl + 1);
        big_mod_rand(Z2, nl, n, nl, mu, nl + 1);
        big_mod_rand(Zd, nl, n, nl, mu, nl + 1);

        p1->X = X1;
        p2->X = X2;
        pd->X = Xd;
        p1->Z = Z1;
        p2->Z = Z2;
        pd->Z = Zd;

        pro_add(p, p1, p2, pd, n, nl, mu, nl + 1);

        big_print(fp, n, nl, "n", NULL);
        big_print(fp, p1->X, nl, "X1", NULL);
        big_print(fp, p1->Z, nl, "Z1", NULL);
        big_print(fp, p2->X, nl, "X2", NULL);
        big_print(fp, p2->Z, nl, "Z2", NULL);
        big_print(fp, pd->X, nl, "Xd", NULL);
        big_print(fp, pd->Z, nl, "Zd", NULL);

        big_print(fp, p->X, nl, "pX", NULL);
        big_print(fp, p->Z, nl, "pZ", NULL);

        fprintf(fp, "qX, qZ := ADDM(X1, Z1, X2, Z2, Xd, Zd, n);\n");
        fprintf(fp, "if qX eq pX then\n");
        fprintf(fp, "trues := trues + 1;\n");
        fprintf(fp, "else\n");
        fprintf(fp, "falses := falses + 1;\n");
        fprintf(fp, "end if;\n");

    }
    fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_add_results.magma\", trues);\n");
    fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_add_results.magma\", falses);\n");

    fclose(fp);
}

void pro_dbl_magma_test(int THRESHOLD) {
    FILE *fp = fopen("/home/ozbayelif/Development/FIWE/ecm/pro_dbl_test.magma", "a");
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    ui_t nl;
    int i;

    fprintf(fp, "clear;\n");
    fprintf(fp, "/****************************************************************************/\n");
    fprintf(fp, "DBLM:=function(X1, Z1, A24)\n");
    fprintf(fp, "return ((X1+Z1)*(X1+Z1))*((X1-Z1)*(X1-Z1)),(((X1+Z1)*(X1+Z1))-((X1-Z1)*(X1-Z1)))*(((X1-Z1)*(X1-Z1))+A24*(((X1+Z1)*(X1+Z1)) - ((X1-Z1)*(X1-Z1))));\n");
    fprintf(fp, "end function;\n");
    fprintf(fp, "trues := 0;\n");
    fprintf(fp, "falses := 0;\n");

    for (i = 0; i < THRESHOLD; i++) {
        nl = (ui_t)(rand() % 100 + 1);
        ui_t n[nl], mu[nl + 1], X1[nl], Z1[nl], A24[nl];

        big_rand(n, nl);
        big_get_mu(mu, n, nl);
        big_mod_rand(X1, nl, n, nl, mu, nl + 1);
        big_mod_rand(Z1, nl, n, nl, mu, nl + 1);
        big_mod_rand(A24, nl, n, nl, mu, nl + 1);

        p1->X = X1;
        p1->Z = Z1;

        pro_dbl(p, p1, A24, n, nl, mu, nl + 1);

        big_print(fp, n, nl, "n", NULL);
        fprintf(fp, "R:=Integers(n);\n");
        big_print(fp, p1->X, nl, "X1", "R");
        big_print(fp, p1->Z, nl, "Z1", "R");
        big_print(fp, A24, nl, "A24", "R");

        big_print(fp, p->X, nl, "pX", "R");
        big_print(fp, p->Z, nl, "pZ", "R");

        fprintf(fp, "qX, qZ := DBLM(X1, Z1, A24);\n");
        fprintf(fp, "if qX eq pX then\n");
        fprintf(fp, "trues := trues + 1;\n");
        fprintf(fp, "else\n");
        fprintf(fp, "falses := falses + 1;\n");
        fprintf(fp, "end if;\n");

    }
    fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_dbl_results.magma\", trues);\n");
    fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_dbl_results.magma\", falses);\n");

    fclose(fp);
}

void pro_ladder_gmp_test(int THRESHOLD) {


    int i, nl, kl, ll, *flag, trues = 0, falses = 0;
    ui n, mu, A24, d, k, l;
    nl = (ui_t)(rand() % 10 + 1), kl = 1, ll = 1;

    cudaMallocManaged(&flag, SIZE * sizeof(int));
    cudaMallocManaged(&n, nl * sizeof(ui_t));
    cudaMallocManaged(&mu, (nl + 1) * sizeof(ui_t));
    cudaMallocManaged(&A24, nl * sizeof(ui_t));
    cudaMallocManaged(&d, SIZE * nl * sizeof(ui_t));
    cudaMallocManaged(&k, SIZE * kl * sizeof(ui_t));
    cudaMallocManaged(&l, SIZE * ll * sizeof(ui_t));


    mpz_t mp_n, mp_p3X, mp_p3Z, mp_p5X, mp_p5Z, mp_Xlk, mp_Xkl;

    MONTG_CURVE c;

    cudaMallocManaged(&c, sizeof(MONTG_CURVE_t));
    cudaMallocManaged(&c->A, sizeof(ui_t) * nl);
    cudaMallocManaged(&c->B, sizeof(ui_t) * nl);

    PRO_POINT p1, p2, p3, p4, p5;

    cudaMallocManaged(&p1, sizeof(PRO_POINT_t));
    cudaMallocManaged(&p1->X, sizeof(ui_t) * nl);
    cudaMallocManaged(&p1->Y, sizeof(ui_t) * nl);
    cudaMallocManaged(&p1->Z, sizeof(ui_t) * nl);

    cudaMallocManaged(&p2, sizeof(PRO_POINT_t));
    cudaMallocManaged(&p2->X, sizeof(ui_t) * nl);
    cudaMallocManaged(&p2->Y, sizeof(ui_t) * nl);
    cudaMallocManaged(&p2->Z, sizeof(ui_t) * nl);

    cudaMallocManaged(&p3, sizeof(PRO_POINT_t));
	cudaMallocManaged(&p3->X, sizeof(ui_t) * nl);
	cudaMallocManaged(&p3->Y, sizeof(ui_t) * nl);
	cudaMallocManaged(&p3->Z, sizeof(ui_t) * nl);

	cudaMallocManaged(&p4, sizeof(PRO_POINT_t));
	cudaMallocManaged(&p4->X, sizeof(ui_t) * nl);
	cudaMallocManaged(&p4->Y, sizeof(ui_t) * nl);
	cudaMallocManaged(&p4->Z, sizeof(ui_t) * nl);

	cudaMallocManaged(&p5, sizeof(PRO_POINT_t));
	cudaMallocManaged(&p5->X, sizeof(ui_t) * nl);
	cudaMallocManaged(&p5->Y, sizeof(ui_t) * nl);
	cudaMallocManaged(&p5->Z, sizeof(ui_t) * nl);

    mpz_init(mp_n);
    mpz_init(mp_p3X);
    mpz_init(mp_p3Z);
    mpz_init(mp_p5X);
    mpz_init(mp_p5Z);
    mpz_init(mp_Xlk);
    mpz_init(mp_Xkl);

    mpz_set_ui(mp_n, 0L);
    mpz_set_ui(mp_p3X, 0L);
    mpz_set_ui(mp_p3Z, 0L);
    mpz_set_ui(mp_p5X, 0L);
    mpz_set_ui(mp_p5Z, 0L);
    mpz_set_ui(mp_Xlk, 0L);
    mpz_set_ui(mp_Xkl, 0L);

	PRO_POINT R0, R0_, R1, R1_;

	cudaMallocManaged(&R0, sizeof(PRO_POINT));
	cudaMallocManaged(&R0_, sizeof(PRO_POINT));
	cudaMallocManaged(&R1, sizeof(PRO_POINT));
	cudaMallocManaged(&R1_, sizeof(PRO_POINT));

	cudaMallocManaged(&R0->X, sizeof(ui_t) * nl);
	cudaMallocManaged(&R0_->X, sizeof(ui_t) * nl);

	cudaMallocManaged(&R1->X, sizeof(ui_t) * nl);
	cudaMallocManaged(&R1_->X, sizeof(ui_t) * nl);

	cudaMallocManaged(&R0->Z, sizeof(ui_t) * nl);
	cudaMallocManaged(&R0_->Z, sizeof(ui_t) * nl);

	cudaMallocManaged(&R1->Z, sizeof(ui_t) * nl);
	cudaMallocManaged(&R1_->Z, sizeof(ui_t) * nl);

    for (i = 0; i < THRESHOLD; i++) {
        big_rand(n, nl);
        n[0]--;

        big_get_mu(mu, n, nl);
        big_rand(k, kl);
        big_rand(l, ll);

        pro_curve_point(d, c, p1, n, nl, mu, nl + 1, flag);
        if(flag[0] != 1){
            i--;
            continue;
        }
        big_get_A24(A24, c->A, n, nl, mu, nl + 1, flag);
        if(flag[0] != 1) {
            i--;
            continue;
        };

        //big_print(stdout, p->X, nl, "X", NULL);

//        proLadder<<<1, 1>>>(p2, p1, R0, R0_, R1, R1, A24, k, kl, n, nl, mu, nl + 1);
//        cudaDeviceSynchronize();
//
//        //big_print(stdout, p2->X, nl, "X", NULL);
//
//        proLadder<<<1, 1>>>(p3, p2, R0, R0_, R1, R1, A24, k, kl, n, nl, mu, nl + 1);
//		cudaDeviceSynchronize();
//
//		//big_print(stdout, p3->X, nl, "X", NULL);
//
//		proLadder<<<1, 1>>>(p4, p1, R0, R0_, R1, R1, A24, k, kl, n, nl, mu, nl + 1);
//		cudaDeviceSynchronize();
//
//		//big_print(stdout, p4->X, nl, "X", NULL);
//
//		proLadder<<<1, 1>>>(p5, p4, R0, R0_, R1, R1, A24, k, kl, n, nl, mu, nl + 1);
//		cudaDeviceSynchronize();

		//big_print(stdout, p5->X, nl, "X", NULL);

        mpz_import(mp_n, nl, -1, 4, 0, 0, n);
        mpz_import(mp_p3X, nl, -1, 4, 0, 0, p3->X);
        mpz_import(mp_p5X, nl, -1, 4, 0, 0, p5->X);
        mpz_import(mp_p3Z, nl, -1, 4, 0, 0, p3->Z);
        mpz_import(mp_p5Z, nl, -1, 4, 0, 0, p5->Z);

        mpz_mul(mp_Xlk, mp_p3X, mp_p5Z);                    // X1*Z2
        mpz_mod(mp_Xlk, mp_Xlk, mp_n);                      // X1*Z2 mod n
        mpz_mul(mp_Xkl, mp_p5X, mp_p3Z);                    // X2*Z1
        mpz_mod(mp_Xkl, mp_Xkl, mp_n);                      // X2*Z1 mod n

        if(mpz_cmp(mp_Xlk, mp_Xkl) == 0) {
            trues++;
        } else {
            falses++;
        }
    }
    printf("TRUE: %d\n", trues);
    printf("FALSE: %d\n", falses);
}

void pro_ladder_magma_test(int THRESHOLD) {
    FILE *fp = fopen("/home/ozbayelif/Development/FIWE/ecm/pro_ladder_test.magma", "a");
    // FILE *fp = stdout;
    MONTG_CURVE c = (MONTG_CURVE)malloc(sizeof(MONTG_CURVE_t) * 1);
    PRO_POINT p = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p1 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p2 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p3 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p4 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    PRO_POINT p5 = (PRO_POINT)malloc(sizeof(PRO_POINT_t) * 1);
    int i, nl, kl, ll, flag;
    nl = (ui_t)5, kl = 1, ll = 1;
    ui_t mu[nl + 1], A24[nl], d[nl], k[kl], l[ll];
    ui_t n[] = {3411243619, 3283606458, 2946840869, 2642350139, 82690173}; // Prime
    big_get_mu(mu, n, nl);

    fprintf(fp, "clear;\n\n");
    fprintf(fp, "load \"/home/ozbayelif/Development/FIWE/ecm/montgomery.m\";\n\n");
    fprintf(fp, "trues:=0;\n\n");
    fprintf(fp, "falses:=0;\n\n");

    big_print(fp, n, nl, "n", NULL);
    fprintf(fp, "F:=GF(n);\n\n");

    for (i = 0; i < THRESHOLD; i++) {
        big_rand(k, kl);
        big_rand(l, ll);
        big_print(fp, k, kl, "k", NULL);
        big_print(fp, l, ll, "l", NULL);

        //pro_curve_point(d, c, p1, n, nl, mu, nl + 1, &flag);
        if(flag != 1) {
            i--;
            continue;
        }
        big_print(fp, c->A, nl, "A", "F");
        big_print(fp, c->B, nl, "B", "F");
        fprintf(fp, "S<X,Y,Z>:=ProjectiveSpace(F,2);\n\n");
        fprintf(fp, "C<X,Y,Z>:=Curve(S,[B*Y^2*Z-(X^3+A*X^2*Z+X*Z^2)]);\n\n");
        fprintf(fp, "W,MtoW:=EllipticCurve(C,C![0,1,0]);\n\n");
        fprintf(fp, "WtoM:=Inverse(MtoW);\n\n");
        big_print(fp, p1->X, nl, "X1", "F");
        big_print(fp, p1->Y, nl, "Y1", "F");
        big_print(fp, p1->Z, nl, "Z1", "F");
        fprintf(fp, "P1:=C![X1,Y1,Z1];\n\n");

        big_get_A24(A24, c->A, n, nl, mu, nl + 1, &flag);
        if(flag != 1) {
            i--;
            continue;
        }
        big_print(fp, A24, nl, "A24", NULL);
        fprintf(fp, "assert A24 eq (A+2)/4;\n\n");

        pro_ladder(p2, p1, A24, k, kl, n, nl, mu, nl + 1);  // p2 = k*P
        big_print(fp, p2->X, nl, "Xk_", NULL);
        big_print(fp, p2->Z, nl, "Zk_", NULL);
        fprintf(fp, "Xk,Zk:=LADDM(X1,Z1,k,A24);\n\n");
        fprintf(fp, "assert (Xk_ eq Xk and Zk_ eq Zk);\n\n");

        pro_ladder(p3, p2, A24, l, ll, n, nl, mu, nl + 1);  // p3 = l*(k*P)
        big_print(fp, p3->X, nl, "Xlk_", NULL);
        big_print(fp, p3->Z, nl, "Zlk_", NULL);
        fprintf(fp, "Xlk,Zlk:=LADDM(Xk,Zk,l,A24);\n\n");
        fprintf(fp, "assert (Xlk_ eq Xlk and Zlk_ eq Zlk);\n\n");

        pro_ladder(p4, p1, A24, l, ll, n, nl, mu, nl + 1);  // p4 = l*P
        big_print(fp, p4->X, nl, "Xl_", NULL);
        big_print(fp, p4->Z, nl, "Zl_", NULL);
        fprintf(fp, "Xl,Zl:=LADDM(X1,Z1,l,A24);\n\n");
        fprintf(fp, "assert (Xl_ eq Xl and Zl_ eq Zl);\n\n");

        pro_ladder(p5, p4, A24, k, kl, n, nl, mu, nl + 1);  // p5 = k*(l*P)
        big_print(fp, p5->X, nl, "Xkl_", NULL);
        big_print(fp, p5->Z, nl, "Zkl_", NULL);
        fprintf(fp, "Xkl,Zkl:=LADDM(Xl,Zl,k,A24);\n\n");
        fprintf(fp, "assert (Xkl_ eq Xkl and Zkl_ eq Zkl);\n\n");

        fprintf(fp, "Xlk_:=F!Xlk_/Zlk_;\n\n");
        fprintf(fp, "Xkl_:=F!Xkl_/Zkl_;\n\n");

        fprintf(fp, "if Xlk_ eq Xkl_ then\n");
        fprintf(fp, "trues:=trues + 1;\n");
        fprintf(fp, "else\n");
        fprintf(fp, "falses:=falses + 1;\n");
        fprintf(fp, "end if;\n");
    }
    fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_ladder_results.magma\", trues);\n");
    fprintf(fp, "Write(\"/home/ozbayelif/Development/FIWE/ecm/pro_ladder_results.magma\", falses);\n");

    fclose(fp);
}

void ecm_gmp_test(int THRESHOLD) {

    int i, nl, res, trues = 0, falses = 0, success = 0, fail = 0, success_type[3] = {0};

    mpz_t mp_n, mp_d, mp_mod;

    mpz_init(mp_n);
    mpz_init(mp_mod);
    mpz_init(mp_d);
    mpz_set_ui(mp_d, 0L);

    ui n, d, returnArr, mu, h_returnArr[0];

    nl = 10;

    ui_t mul = nl + 1;

    ui_t h_n[nl], h_d[SIZE * nl], h_mu[mul], bl;

   // generate_B_smooth(h_n, nl);
    big_rand(h_n, nl);
    h_n[0]--;
    big_get_mu(h_mu, h_n, nl);

    bl = (nl >> 1) + 1;

    cudaMalloc(&n, sizeof(ui_t) * nl);
    cudaMalloc(&mu, sizeof(ui_t) * mul);

    cudaMemcpy(n, h_n, sizeof(ui_t) * nl, cudaMemcpyHostToDevice);
    cudaMemcpy(mu, h_mu, sizeof(ui_t) * mul, cudaMemcpyHostToDevice);

    for(i = 0; i < SIZE * nl; i++)
    	h_d[i] = 0;

    cudaMalloc(&d, sizeof(ui_t) * nl * SIZE);
    cudaMemcpy(d, h_d, sizeof(ui_t) * nl, cudaMemcpyHostToDevice);

    cudaMalloc(&returnArr, sizeof(ui_t) * SIZE);

    size_t pValue;

    cudaDeviceSetLimit(cudaLimitStackSize, 10000 * sizeof(ui_t));
	cudaDeviceGetLimit(&pValue, cudaLimitStackSize);

	cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, 10000 * sizeof(ui_t));
	cudaDeviceGetLimit(&pValue, cudaLimitDevRuntimePendingLaunchCount);

	cudaDeviceGetLimit(&pValue, cudaLimitMallocHeapSize);
	cudaDeviceSetLimit(cudaLimitMallocHeapSize , pValue * 25);
	cudaDeviceGetLimit(&pValue, cudaLimitMallocHeapSize);

    fiwEcm<<<BLOCKNUM, THREADNUM>>>(d, n, nl, mu, mul, returnArr);
    cudaDeviceSynchronize();

    cudaMemcpy(h_returnArr, returnArr, sizeof(ui_t) * SIZE, cudaMemcpyDeviceToHost);

    printf("%u\t", h_returnArr[0]);
    printf("%u\n", h_returnArr[1]);

    cudaError_t error = cudaGetLastError();
	if(error != cudaSuccess)
		fprintf(stderr, "ERROR: %s\n", cudaGetErrorString(error));

    big_print(stdout, h_n, nl, "n", NULL);
    cudaMemcpy(h_d, d, sizeof(ui_t) * SIZE * nl, cudaMemcpyDeviceToHost);

    i = 0;
 	int j;
    int t_index = 0;
//
    while(i < SIZE * nl){
    	printf("t_index: %d: ", t_index);
    	t_index++;

    	j = i;

    	while(j < (i + nl)){
    		printf("%u ", h_d[j]);
    		j++;
    	}

    	printf("\n");
    	i += nl;

    }


//    for(i = 0; i < THRESHOLD; i++) {
//
//        ui_t n[nl], d[nl];
//
//        //big_print(stdout, n, nl, "n", NULL);
//
//        int ret = ecm(d, n, nl);
//        if(ret != 0) {
//            success++;
//            success_type[ret - 1]++;
//            mpz_import(mp_n, nl, -1, 4, 0, 0, n);
//            mpz_import(mp_d, nl, -1, 4, 0, 0, d);
//            mpz_mod(mp_mod, mp_n, mp_d);
//            res = mpz_cmp_ui(mp_mod, 0L);
//            if(res == 0) {
//                big_print(stdout, d, nl, "d", NULL);
//                true++;
//            } else {
//                false++;
//                printf("False at index: %d\n", i);
//            }
//        } else {
//            printf("fail \n");
//            fail++;
//        }
//    }
//    printf("TRUE: %d\n", trues);
//    printf("FALSE: %d\n", falses);
//    printf("SUCCESS: %d\n", success);
//    printf("FAIL: %d\n", fail);
//    printf("FOUND IN pro_curve_point: %d\n", success_type[0]);
//    printf("FOUND IN ladder: %d\n", success_type[1]);
//    printf("FOUND IN A24: %d\n", success_type[2]);
}
