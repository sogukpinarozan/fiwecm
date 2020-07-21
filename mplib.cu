/**
 * \Created on: Mar 18, 2020
 * \Author: Asena Durukan, Asli Altiparmak, Elif Ozbay, Hasan Ozan Sogukpinar, Nuri Furkan Pala
 * \file mplib.cu
 * \brief Implementation of mplib.h library.
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "mplib.h"


#define ufvg(u, f, v, g, c, h, l) __asm__( \
	"imulq %3;" \
	"movq %%rax, %%rsi;" \
	"movq %%rdx, %%rcx;" \
	"movq %4, %%rax;" \
	"imulq %5;" \
	"addq %%rsi, %%rax;" \
	"adcq %%rcx, %%rdx;" \
	\
	"movq %6, %%r10;" \
	"btq $63, %%r10;" \
	"movq $0, %%r11;" \
	"sbbq $0, %%r11;" \
	"addq %%r10, %%rax;" \
	"adcq %%r11, %%rdx;" \
	\
	"movq $0x3fffffffffffffff, %%rdi;" \
	"movq $0xc000000000000000, %%r11;" \
	"movq %%rdx, %%r10;" \
	"andq %%r11, %%r10;" \
	"shld $2, %%rax, %%rdx;" \
	"andq %%rdi, %%rax;" \
	"orq %%r10, %%rax;" \
	: "=a"((l)), "=d"((h)) \
	: "a"((u)), "m"((f)), "m"((v)), "m"((g)), "m"((c)) \
)


void big_rand(ui z, ui_t l) {
	int i;

	for (i = 0; i < l; i++) {
		z[i] = ((ui_t) rand()) * ((ui_t) rand()) * ((ui_t) rand())
				* ((ui_t) rand());
	}
}

void big_mod_rand(ui z, ui_t l, ui n, ui_t nl, ui mu, ui_t mul) {
	ui_t z_[2 * nl];
	int i;

	for (i = 0; i < l; i++) {
		z_[i] = ((ui_t) rand()) * ((ui_t) rand()) * ((ui_t) rand())
				* ((ui_t) rand());
	}
	for (i = l; i < 2 * nl; i++) {
		z_[i] = 0L;
	}
	barret_reduction(z, z_, 2 * nl, n, nl, mu, mul);
}

void big_print(FILE *fp, ui a, ui_t al, const char *s, const char *R) {
	if (R != NULL) {
		fprintf(fp, "%s := %s!(", s, R);
	} else {
		fprintf(fp, "%s := (", s);
	}
	fprintf(fp, "%u", a[0]);
	for (int i = 1; i < al; i++) {
		fprintf(fp, " + %u * (2^%d)^%d", a[i], W, i);
	}
	fprintf(fp, ");\n\n");
}

void big_is_equal(int *z, ui a, ui b, ui_t l) {
	int i;
	*z = 1;
	for (i = 0; i < l; i++) {
		if (a[i] != b[i]) {
			*z = 0;
		}
	}
}

void big_is_equal_ui(int *z, ui a, ui_t al, ui_t b) {
	int i;
	*z = 1;
	if (a[0] != b) {
		*z = 0;
	}
	for (i = 1; i < al; i++) {
		if (a[i] != 0) {
			*z = 0;
		}
	}
}

void big_add(ui z, ui a, ui_t al, ui b, ui_t bl) {
	int i;
	ui_t carry_bit = 0;

	for (i = 0; i < al; i++) {
		z[i] = a[i] + b[i] + carry_bit;
		if (z[i] < a[i]) {
			carry_bit = 1;
		} else if (z[i] > a[i]) {
			carry_bit = 0;
		}
	}
}

void big_mod_add(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl, ui mu,
		ui_t mul) {
	int i;
	ui_t z_[2 * nl], carry_bit = 0;

	for (i = al + 1; i < 2 * nl; i++) {
		z_[i] = 0L;
	}
	for (i = 0; i < al; i++) {
		z_[i] = a[i] + b[i] + carry_bit;
		if (z_[i] < a[i]) {
			carry_bit = 1;
		} else if (z_[i] > a[i]) {
			carry_bit = 0;
		}
	}
	z_[al] = carry_bit;
	barret_reduction(z, z_, 2 * nl, n, nl, mu, mul);
}

void big_sub(ui z, int *d, ui a, ui_t al, ui b, ui_t bl) {
	int i;
	ui_t borrow_bit = 0;

	for (i = 0; i < al; i++) {
		z[i] = a[i] - b[i] - borrow_bit;
		if (z[i] < a[i]) {
			borrow_bit = 0;
		} else if (z[i] > a[i]) {
			borrow_bit = 1;
		}
	}
	*d = borrow_bit;
}

void big_mod_sub(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl) {
	int i;
	ui_t z_[nl], borrow_bit = 0;

	for (i = 0; i < al; i++) {
		z_[i] = a[i] - b[i] - borrow_bit;
		if (z_[i] < a[i]) {
			borrow_bit = 0;
		} else if (z_[i] > a[i]) {
			borrow_bit = 1;
		}
	}
	if (borrow_bit) {
		big_add(z, z_, nl, n, nl);
	} else {
		big_cpy(z, z_, 0, nl);
	}
}

void big_mul(ui z, ui a, ui_t al, ui b, ui_t bl) {
	int i, j;
	ui_t u, v;
	uni_t uv;

	for (i = 0; i <= al; i++) {
		z[i] = 0;
	}
	for (i = 0; i < al; i++) {
		u = 0;
		for (j = 0; j < bl; j++) {
			uv = (uni_t) z[i + j] + (uni_t) a[i] * (uni_t) b[j] + (uni_t) u;
			u = uv >> W;
			v = uv & 0xFFFFFFFF;  // TODO: W != 32?
			z[i + j] = v;
		}
		z[i + bl] = u;
	}
}

void big_mod_mul(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl, ui mu,
		ui_t mul) {
	int i, j;
	ui_t u, v, z_[2 * nl];
	uni_t uv;

	for (i = al + bl; i < 2 * nl; i++) {
		z_[i] = 0;
	}
	for (i = 0; i <= al; i++) {
		z_[i] = 0;
	}
	for (i = 0; i < al; i++) {
		u = 0;
		for (j = 0; j < bl; j++) {
			uv = (uni_t) z_[i + j] + (uni_t) a[i] * (uni_t) b[j] + (uni_t) u;
			u = uv >> W;
			v = uv & 0xFFFFFFFF; // TODO: W != 32?
			z_[i + j] = v;
		}
		z_[i + bl] = u;
	}
	barret_reduction(z, z_, 2 * nl, n, nl, mu, mul);
}

void big_get_mu(ui mu, ui n, ui_t nl) {
	mpz_t mp_n, mp_b2k, mp_mu;

	mpz_init(mp_n);
	mpz_init(mp_b2k);
	mpz_init(mp_mu);

	mpz_set_ui(mp_b2k, 0L);
	mpz_set_ui(mp_mu, 0L);

	mpz_import(mp_n, nl, -1, 4, 0, 0, n);
	mpz_add_ui(mp_b2k, mp_b2k, 1);
	mpz_mul_2exp(mp_b2k, mp_b2k, W * (2 * nl));
	mpz_fdiv_q(mp_mu, mp_b2k, mp_n);
	mpz_export(mu, NULL, -1, 4, 0, 0, mp_mu);
}

__device__ void bigGetA24(ui z, ui A, ui n, ui_t nl, ui mu, ui_t mul, int *flag){

	int t_index = threadIdx.x + (blockIdx.x * blockDim.x);

	if(t_index < SIZE){

		ui c_2 = new ui_t[nl];
		ui c_4 = new ui_t[nl];
		ui A2 = new ui_t[nl];
		ui ic_4 = new ui_t[nl];
		ui tempMultResult = new ui_t[nl];
		ui_t isEqual = 0;

		int i;

		c_2[0] = 2L;
		c_4[0] = 4L;

		for(i = 1; i < nl; i++) {
			c_2[i] = 0L;
			c_4[i] = 0L;
		}

		bigModAdd(A2, A, nl, c_2, nl, n, nl, mu, mul);
		bigInvert(ic_4, c_4, nl, n, nl, mu, mul);
		bigMul(tempMultResult, ic_4, nl, c_4, nl);
		bigIsEqual(&isEqual, tempMultResult, n, nl);

//		bigPrint(ic_4, nl, "ic4");
//		bigPrint(c_4, nl, "c4");
//		bigPrint(tempMultResult, nl, "tempMultResult");
//		bigPrint(n, nl, "n");
//		printf("isEqual: %d\n", isEqual);

		if(isEqual != 1){ //Inverse exists
			bigModMul(z, A2, nl, ic_4, nl, n, nl, mu, mul);
			*flag = 1;
		}
		else{
			binaryGCD(z, c_4, nl, n, nl);
			*flag = 0;
		}
	}
}

void big_get_A24(ui z, ui A, ui n, ui_t nl, ui mu, ui_t mul, int *flag) {
	ui_t c_2[nl], c_4[nl], A2[nl], ic_4[nl];
	int i, ret;

	c_2[0] = 2L;
	c_4[0] = 4L;
	for (i = 1; i < nl; i++) {
		c_2[i] = 0L;
		c_4[i] = 0L;
	}
	big_mod_add(A2, A, nl, c_2, nl, n, nl, mu, mul);
	ret = big_invert(ic_4, c_4, nl, n, nl);
	if (ret) { // Inverse exists
		big_mod_mul(z, A2, nl, ic_4, nl, n, nl, mu, mul);
		*flag = 1;
	} else { // Inverse does not exist
		big_gcd(z, nl, c_4, nl, n, nl);
		*flag = 0;
	}
}

// ml = 2 * nl
void barret_reduction(ui z, ui m, ui_t ml, ui n, ui_t nl, ui mu, ui_t mul) { // Calculate m mod n
	ui_t k = nl, md[k + 1], mdmu[mul + k + 1], q[mul], mm[k + 1], qn[mul + nl],
			qnm[k + 1], r2[k + 1], r3[k + 1];
	int i, b;

	big_cpy(md, m, k - 1, k + 1); // md = m / b^(k - 1)
	big_mul(mdmu, md, k + 1, mu, mul); // mdmu = md * mu
	big_cpy(q, mdmu, k + 1, mul); // q = (m / b^(k - 1) * mu) / b^(k + 1)
	big_cpy(mm, m, 0, k + 1); // mm = m mod b^(k + 1)
	big_mul(qn, q, mul, n, nl); // qn = q * n
	big_cpy(qnm, qn, 0, k + 1); // qnm = qn mod b^(k + 1)
	big_sub(r3, &i, mm, k + 1, qnm, k + 1); // r3 = mm - qnm
	big_cpy(z, r3, 0, k);
	big_sub(r2, &b, r3, nl, n, nl); // while r >= n do: r <- r - n
	r2[nl] = r3[nl] - b;
	while (!(r2[nl] >> (W - 1))) {
		big_cpy(z, r2, 0, k);
		big_cpy(r3, r2, 0, k + 1);
		big_sub(r2, &b, r3, nl, n, nl);
		r2[nl] = r3[nl] - b;
	}
}

void big_gcd(ui d, ui_t dl, ui a, ui_t al, ui b, ui_t bl) {
	mpz_t mp_a, mp_b, mp_d;
	int i;

	for (i = 0; i < dl; i++) {
		d[i] = 0L;
	}
	mpz_init(mp_a);
	mpz_init(mp_b);
	mpz_init(mp_d);

	mpz_set_ui(mp_a, 0L);
	mpz_set_ui(mp_b, 0L);
	mpz_set_ui(mp_d, 0L);

	mpz_import(mp_a, al, -1, 4, 0, 0, a);
	mpz_import(mp_b, bl, -1, 4, 0, 0, b);
	mpz_gcd(mp_d, mp_a, mp_b);
	mpz_export(d, NULL, -1, 4, 0, 0, mp_d);
}

int big_invert(ui z, ui a, ui_t al, ui b, ui_t bl) {

	int i, ret;
	mpz_t mp_z, mp_a, mp_b;

	mpz_init(mp_z);
	mpz_init(mp_a);
	mpz_init(mp_b);

	mpz_import(mp_a, al, -1, 4, 0, 0, a);
	mpz_import(mp_b, bl, -1, 4, 0, 0, b);
	mpz_set_ui(mp_z, 0L);

	for (i = 0; i < bl; i++) {
		z[i] = 0L;
	}
	ret = mpz_invert(mp_z, mp_a, mp_b);
	mpz_export(z, NULL, -1, 4, 0, 0, mp_z);       // iY2Z = Inv(Y^2Z)

	return ret;
}

//montgomery simultaneous inversion algorithm

__device__ void bigPrint(ui a, ui_t al, const char *s) {
	printf("%s := (", s);
	printf("%u", a[0]);
	for (int i = 1; i < al; i++) {
		printf(" + %u * (2^%d)^%d", a[i], W, i);
	}
	printf(");\n\n");
}

__global__ void bigPrintG(ui a, ui_t al, char *s) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {
	printf("%s := (", s);
		printf("%u", a[0]);
		for (int i = 1; i < al; i++) {
			printf(" + %u * (2^%d)^%d", a[i], W, i);
		}
		printf(");\n\n");
	}
}

void memoryAllocationGPU(ui *deviceArray, ui_t arraySize) {
	//Memory Allocation for GPU
	cudaMalloc(deviceArray, arraySize * sizeof(ui_t));
	if (deviceArray == NULL)
		printf("deviceArray no space");
}

__device__ void bigCpy(ui z, ui a, ui_t start, ui_t end) {

	//md, m, k - 1, k + 1, 2 * nl
	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {

//		int firstIndexZ = t_index * end, lastIndexZ = firstIndexZ + end;
//		int firstIndexA = t_index * arraySize + start;


		bigCopy(z, 0, end, a, start);

//		for(i = (firstIndexZ), j = (firstIndexA); i < (lastIndexZ); i++, j++) { \
//				z[i] = a[j]; \
//			} \
	}
}
}

__device__ void bigMul(ui z, ui a, ui_t al, ui b, ui_t bl) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {

		ui_t u, low, high;
		int i, j;

		//mdmu, md, k + 1, mu, mul

//		int firstIndexA = t_index * al, lastIndexA = firstIndexA + al;
//		int firstIndexB = t_index * bl, lastIndexB = firstIndexB + bl;
//		int firstIndexZ = t_index * (al + bl);

		for (i = 0; i <= al; i++) {
			z[i] = 0L;
		}

		for (i = 0; i < al; i++) {
			u = 0;
			for (j = 0; j < bl; j++) {
				low = 0;
				high = 0;
				__mul_lo(low, a[i], b[j]);
				__mul_hi(high, a[i], b[j]);
				__add_cc(low, z[i+j], low);
				__addcy2(high);
				__add_cc(low, u, low);
				__addcy2(high);

				z[i + j] = low;
				u = high;
			}
			z[i + bl] = u;
		}
	}

}

//TODO: macros
__device__ void bigSub(ui z, ui controlBit, ui a, ui_t al, ui b, ui_t bl) {


	int t_index = threadIdx.x + (blockIdx.x * blockDim.x);

	if (t_index < SIZE) {

		//int firstIndexA = t_index * al, lastIndexA = firstIndexA + al;
		ui_t borrow_bit = 0;
		int i;

		for (i = 0; i < al; i++) {
			z[i] = a[i] - b[i] - borrow_bit;
			if (z[i] < a[i]) {
				borrow_bit = 0;
			} else if (z[i] > a[i]) {
				borrow_bit = 1;
			}
		}

		*controlBit = borrow_bit;
	}

}

__device__ void bigAdd(ui z, ui a, ui_t al, ui b, ui_t bl) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {

		int i;

		//int firstIndexA = t_index * al, lastIndexA = firstIndexA + al;

		__add_cc(z[0], a[0], b[0]);
		for (i = 1; i < al; i++) {
			__addc_cc(z[i], a[i], b[i]);
		}
		__addcy(z[al]);
	}
}

__device__ void bigModSub(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl) {

	int t_index = threadIdx.x + (blockIdx.x * blockDim.x);

	if (t_index < SIZE) {

		ui z_ = new ui_t[2 * nl];
		ui_t borrow_bit = 0;
		int firstIndexA = t_index * al, lastIndexA = firstIndexA + al;
		int i;

		for (i = firstIndexA; i < lastIndexA; i++) {
			z_[i] = a[i] - b[i] - borrow_bit;
			if (z_[i] < a[i]) {
				borrow_bit = 0;
			} else if (z_[i] > a[i]) {
				borrow_bit = 1;
			}
		}

		if (borrow_bit) {
			bigAdd(z, z_, nl, n, nl);
		} else {
			bigCopy(z, 0, nl, z_, 0);
		}
	}

}

__device__ float gpuGenerate(curandState* globalState) {

	int t_index = threadIdx.x + (blockIdx.x * blockDim.x);

	curandState localState = *globalState;
	float RANDOM = curand_uniform(&localState);
	*globalState = localState;
	return RANDOM;

}



__device__ void bigModRand(ui z, ui_t l, ui n, ui_t nl, ui mu, ui_t mul, curandState* globalState) {

	int t_index = threadIdx.x + (blockIdx.x * blockDim.x);

	if (t_index < SIZE) {
		int i;
		ui z_ = new ui_t[2 * nl];

		curand_init((int)clock() - t_index, t_index, 0, globalState);

		for (i = 0; i < l; i++) {
			z_[i] = (ui_t) (gpuGenerate(globalState) * 4294967295);
		}
		for (i = l; i < nl + l; i++) {
			z_[i] = 0L;
		}

		barretReduction(z, z_, 2 * nl, n, nl, mu, mul);

	}
}

__device__ void bigModMul(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl,
		ui mu, ui_t mul) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {
		int i, j;
		ui z_ = new ui_t[2 * nl];
		ui_t u, low, high;



		for (i = al + bl; i < 2 * nl; i++) {
			z_[i] = 0L;
		}

		for (i = 0; i <= al; i++) {
			z_[i] = 0L;
		}

		for (i = 0; i < al; i++) {
			u = 0;
			for (j = 0; j < bl; j++) {
				low = 0L;
				high = 0L;
				__mul_lo(low, a[i], b[j]);
				__mul_hi(high, a[i], b[j]);
				__add_cc(low, z_[i+j], low);
				__addcy2(high);
				__add_cc(low, u, low);
				__addcy2(high);

				z_[i + j] = low;
				u = high;
			}
			z_[i + bl] = u;
		}
		barretReduction(z, z_, 2 * nl, n, nl, mu, mul);
	}
}

__device__ void bigModMulD(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl,
		ui mu, ui_t mul) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {
		int i, j;
		ui z_ = new ui_t[2 * nl];
		ui_t u, low, high;

		int firstIndexA = t_index * al, lastIndexA = firstIndexA + al;
		int firstIndexB = t_index * bl, lastIndexB = firstIndexB + bl;
		int firstIndexAB = t_index * (al + bl);

		if(t_index == 0){
			firstIndexAB += (al + bl);
		}

		int firstIndexZ_ = t_index * 2 * nl, lastIndexZ_ = firstIndexZ_	+ (2 * nl);

		for (i = firstIndexAB; i < lastIndexZ_; i++) {
			z_[i] = 0L;
		}

		for (i = firstIndexZ_; i <= lastIndexA; i++) {
			z_[i] = 0;
		}

		for (i = firstIndexA; i < lastIndexA; i++) {
			u = 0;
			for (j = firstIndexB; j < lastIndexB; j++) {
				low = 0;
				high = 0;
				__mul_lo(low, a[i], b[j]);
				__mul_hi(high, a[i], b[j]);
				__add_cc(low, z_[i+j], low);
				__addcy2(high);
				__add_cc(low, u, low);
				__addcy2(high);

				z_[i + j] = low;
				u = high;
			}
			z_[i + lastIndexB] = u;
		}

		barretReduction(z, z_, 2 * nl, n, nl, mu, mul);
	}
}

__device__ void bigModAdd(ui z, ui a, ui_t al, ui b, ui_t bl, ui n, ui_t nl,
		ui mu, ui_t mul) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {



		int i;

		//TODO: Why start from 1?
		ui z_ = new ui_t[2 * nl];
		for (i = 1; i < 2 * nl; i++) {
			z_[i] = 0L;
		}

		__add_cc(z_[0], a[0], b[0]);
		for (i = 1; i < al; i++) {
			__addc_cc(z_[i], a[i], b[i]);
		}
		__addcy(z_[al]);



		barretReduction(z, z_, 2 * nl, n, nl, mu, mul);

	}
}


__device__ void bigIsEqualUi(ui z, ui a, ui_t al, ui_t b) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {
		int i;
		//int firstIndexA = t_index * al, lastIndexA = firstIndexA + al;
		*z = 1;

		if (a[0] != b) {
			*z = 0;
		} else {
			if (al > 1) {
				for (i = 1; i < al; i++) {
					if (a[i] != 0) {
						*z = 0;
					}
				}
			}
		}
	}
}

__device__ void bigIsEqualUiD(ui z, ui a, ui_t al, ui_t b) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {

		int i;
		int firstIndexA = t_index * al, lastIndexA = firstIndexA + al;
		*z = 1;

		if (a[firstIndexA] != b) {
			*z = 0;
		} else {
			if (al > 1) {
				for (i = firstIndexA + 1; i < lastIndexA; i++) {
					if (a[i] != 0L) {
						*z = 0;
						break;
					}
				}
			}
		}

	}
}


__device__ void bigIsEqual(ui z, ui a, ui b, ui_t l) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {

		//int firstIndexA = t_index * l, lastIndexA = firstIndexA + l;
		int i;
		*z = 1;

		for (i = 0; i < l; i++) {
			if (a[i] != b[i]) {
				*z = 0;
				break;
			}
		}
	}

}

__device__ void bigIsEqualD(ui z, ui a, ui b, ui_t l) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {

		int firstIndexA = t_index * l, lastIndexA = firstIndexA + l;
		int i, indexB = 0;
		*z = 1;

		for (i = firstIndexA; i < lastIndexA; i++) {
			if (a[i] != b[indexB]) {
				*z = 0;
				break;
			}
			indexB++;
		}
	}

}

__device__ void barretReduction(ui z, ui m, ui_t ml, ui n, ui_t nl, ui mu,
		ui_t mul) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {
		ui_t k = nl;
		ui md = new ui_t[k + 1];

		ui mdmu = new ui_t[mul + k + 1];
		ui q = new ui_t[mul];
		ui mm = new ui_t[k + 1];
		ui qn = new ui_t[mul + nl];
		ui qnm = new ui_t[k + 1];
		ui r2 = new ui_t[k + 1];
		ui r3 = new ui_t[k + 1];
		ui_t controlBit;

		bigCopy(md, 0, k + 1, m, k - 1);
		bigMul(mdmu, md, k + 1, mu, mul); // mdmu = md * mu
		bigCopy(q, 0, mul, mdmu, k + 1);
		bigCopy(mm, 0, k + 1, m, 0);
		bigMul(qn, q, mul, n, nl); // qn = q * n
		bigCopy(qnm, 0, k + 1, qn, 0);
		bigSub(r3, &controlBit, mm, k + 1, qnm, k + 1); // r3 = mm - qnm
		bigCopy(z, 0, k, r3, 0);
		bigSub(r2, &controlBit, r3, nl, n, nl); // while r >= n do: r <- r - n
		r2[nl] = r3[nl] - controlBit;

		while ((!(r2[nl] >> (W - 1)))) { //pozitifse giriyor

			bigCopy(z, 0, k, r2, 0);
			bigCopy(r3, 0, k + 1, r2, 0);
			bigSub(r2, &controlBit, r3, nl, n, nl);
			r2[nl] = r3[nl]- controlBit;
		}
	}

}

__device__ void binaryGCD(ui d, ui a, ui_t al, ui n, ui_t nl) {

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {

		ui u = new ui_t[al];
		ui v = new ui_t[nl];

		ui u_temp = new ui_t[al];
		ui v_temp = new ui_t[nl];
		ui sub_temp = new ui_t[al];
		ui_t is_zero;
		ui_t controlBit;

		//TODO: check either of them is zero

		int firstIndexD = t_index * nl, lastIndexD = firstIndexD + nl;
		int i, lastBit, shift = 0;

		bigCopy(u, 0, al, a, 0); 			// u <- a
		bigCopy(v, 0, nl, n, 0); 			// v <- n


		bigCopy(u_temp, 0, al, a, 0); 		// u <- a
		bigCopy(v_temp, 0, nl, n, 0); 		// v <- n

		//Find number of 2 factors
	    while (((u_temp[0] | v_temp[0]) & 1) == 0) {

	    	shift++;

	        u_temp[0] = u_temp[0] >> 1;
			for(i = 1; i < al; i++){
				lastBit = (u_temp[i] << (W - 1)) >> (W - 1);
				u_temp[i - 1] |= lastBit << (W - 1);
				u_temp[i] = u_temp[i] >> 1;

			}
		    v_temp[0] = v_temp[0] >> 1;
			for(i = 1; i < nl; i++){
				lastBit = (v_temp[i] << (W - 1)) >> (W - 1);
				v_temp[i - 1] |= lastBit << (W - 1);
				v_temp[i] = v_temp[i] >> 1;

			}
	    }

		bigIsEqualUi(&is_zero, u, al, 0L);
		while(!is_zero){
			while(u[0] % 2 == 0){

				u[0] = u[0] >> 1;
				for(i = 1; i < al; i++){

					lastBit = (u[i] << (W - 1)) >> (W - 1);
					u[i - 1] |= lastBit << (W - 1);
					u[i] = u[i] >> 1;

				}
			}
			bigSub(sub_temp, &controlBit, u, al, v, nl);
			if(controlBit == 1){ 	// u < v
				// [u, v] <- [v, u]
				bigCopy(u_temp, 0, al, u, 0);
				bigCopy(u, 0, al, v, 0);
				bigCopy(v, 0, nl, u_temp, 0);
			}
			// u <- u - v
			bigSub(sub_temp, &controlBit, u, al, v, nl);
			bigCopy(u, 0, al, sub_temp, 0);

			bigIsEqualUi(&is_zero, u, al, 0L);
		}

		bigCopy(d, firstIndexD, lastIndexD, v, 0); 	// return v
		//printf("t_index(%u %d %d) %d: %u %u %u %u && %u %u %u %u\n", nl, firstIndexD, lastIndexD, t_index, d[firstIndexD], d[firstIndexD + 1], d[firstIndexD + 2], d[firstIndexD + 3], v[0], v[1], v[2], v[3]);

		//TODO: what if overflow?
		d[firstIndexD] <<= shift;


	}
}

//Binary algorithm for inversion in mod n
__device__ void bigInvert(ui z, ui a, ui_t al, ui n, ui_t nl, ui mu, ui_t mul){

	int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	if (t_index < SIZE) {

		int i;
		ui n_temp = new ui_t[nl + 1];
		ui u = new ui_t[al + 1];
		ui v = new ui_t[nl + 1];
		ui x1 = new ui_t[nl + 1];
		ui x2 = new ui_t[nl + 1];
		ui x1_temp = new ui_t[nl + 2];
		ui x2_temp = new ui_t[nl + 2];
		ui zero_num = new ui_t[nl + 1];
		ui sub_temp = new ui_t[nl + 1];
		ui sub_temp2 = new ui_t[nl];

		int addControl = 0;

		al+= 1; nl += 1;

		ui_t is_one_u = 0, is_one_v = 0, controlBit = 0, lastBit = 0, is_negative_x1 = 0, is_negative_x2 = 0;

		bigCopy(n_temp, 0, nl - 1, n, 0);
		n_temp[nl - 1] = 0;

		//u <- a
		bigCopy(u, 0, al - 1, a, 0);
		a[al - 1] = 0;
		//v <- p
		bigCopy(v, 0, nl - 1, n_temp, 0);
		v[nl - 1] = 0;

		//x1 <- 1, x2 <- 0
		x1[0] = 1;	x2[0] = 0;

		for(i = 1; i < nl; i++){
			x1[i] = 0;
			x2[i] = 0;
		}

		for(i = 0; i < nl; i++)
			zero_num[i] = 0;

		bigIsEqualUi(&is_one_u, u, al, 1L);
		bigIsEqualUi(&is_one_v, v, nl, 1L);


		while(is_one_u != 1 && is_one_v != 1){

			//While u is even
			while(u[0] % 2 == 0){

				// u <- u / 2
				u[0] = u[0] >> 1;
				for(i = 1; i < al; i++){
					lastBit = (u[i] << (W - 1)) >> (W - 1);
					u[i - 1] |= lastBit << (W - 1);
					u[i] = u[i] >> 1;
				}

				// if x1 is even
				if(x1[0] % 2 == 0){

					// x1 <- x1/2
					x1[0] = x1[0] >> 1;
					for(i = 1; i < al; i++){
						lastBit = (x1[i] << (W - 1)) >> (W - 1);
						x1[i - 1] |= lastBit << (W - 1);
						x1[i] = x1[i] >> 1;
					}


				}
				else {
					// x1 <- (x1 + p) / 2

					if(is_negative_x1 == 1){
						bigSub(sub_temp, &controlBit, n_temp, nl, x1, nl);
						if(controlBit == 1){
							bigSub(sub_temp2, &controlBit, zero_num, nl, sub_temp, nl);
							bigCopy(x1_temp, 0, nl, sub_temp2, 0);
						}
						else{
							is_negative_x1 = 0;
							bigCopy(x1_temp, 0, nl, sub_temp, 0);
						}
					}
					else{
						bigAdd(x1_temp, x1, nl, n_temp, nl);
						is_negative_x1 = 0;
					}

					x1_temp[0] = x1_temp[0] >> 1;
					for(i = 1; i < al + 1; i++){
						lastBit = (x1_temp[i] << (W - 1)) >> (W - 1);
						x1_temp[i - 1] |= lastBit << (W - 1);
						x1_temp[i] = x1_temp[i] >> 1;
					}

					x1_temp[nl] = 0;
					bigCopy(x1, 0, nl, x1_temp, 0);

				}

			}

			//while v is even
			while(v[0] % 2 == 0){

				// v <- v / 2
				v[0] = v[0] >> 1;
				for(i = 1; i < al; i++){
					lastBit = (v[i] << (W - 1)) >> (W - 1);
					v[i - 1] |= lastBit << (W - 1);
					v[i] = v[i] >> 1;
				}


				//if x2 is even
				if(x2[0] % 2 == 0){

					// x2 <- x2 / 2

					// -8 / 2 => 4294967296 - (4294967288)
					x2[0] = x2[0] >> 1;
					for(i = 1; i < al; i++){
						lastBit = (x2[i] << (W - 1)) >> (W - 1);
						x2[i - 1] |= lastBit << (W - 1);
						x2[i] = x2[i] >> 1;
					}

				}
				else {
					// x2 <- (x2 + p) / 2

					if(is_negative_x2 == 1){
						bigSub(sub_temp, &controlBit, n_temp, nl, x2, nl);
						if(controlBit == 1){
							bigSub(sub_temp2, &controlBit, zero_num, nl, sub_temp, nl);
							bigCopy(x2_temp, 0, nl, sub_temp2, 0);
						}
						else{
							bigCopy(x2_temp, 0, nl, sub_temp, 0);
							is_negative_x2 = 0;
						}

					}
					else{
						bigAdd(x2_temp, x2, nl, n_temp, nl);
						is_negative_x2 = 0;
					}


					x2_temp[0] = x2_temp[0] >> 1;
					for(i = 1; i < al + 1; i++){
						lastBit = (x2_temp[i] << (W - 1)) >> (W - 1);
						x2_temp[i - 1] |= lastBit << (W - 1);
						x2_temp[i] = x2_temp[i] >> 1;
					}

					x2_temp[al] = 0;
					bigCopy(x2, 0, al, x2_temp, 0);

				}
			}

			// if u >= v
			bigSub(sub_temp, &controlBit, u, al, v, nl);
			if(controlBit != 1){

				// u <- u - v
				bigCopy(u, 0, al, sub_temp, 0);


				// x1 <- x1 - x2
				//9 - (-4)
				//TODO: carry bit might overflow
				if(is_negative_x1 == 0 && is_negative_x2 == 0){ // Both Pozitive
					bigSub(sub_temp, &controlBit, x1, nl, x2, nl);
					if(controlBit == 1){
						is_negative_x1 = 1;
						bigSub(sub_temp2, &controlBit, zero_num, nl, sub_temp, nl);
						bigCopy(sub_temp, 0, nl, sub_temp2, 0);
					}
				}
				else if(is_negative_x1 == 1 && is_negative_x2 == 0){
					bigAdd(x1_temp, x1, nl, x2, nl);
					addControl = 1;
				}
				else if(is_negative_x1 == 0 && is_negative_x2 == 1){
					bigAdd(x1_temp, x1, nl, x2, nl);
					addControl = 1;
				}
				else{ // Both negative
					bigSub(sub_temp, &controlBit, x1, nl, x2, nl);
					if(controlBit == 1){
						bigSub(sub_temp2, &controlBit, zero_num, nl, sub_temp, nl);
						bigCopy(sub_temp, 0, nl, sub_temp2, 0);
						is_negative_x1 = 0;
					}

				}

				if(addControl == 1){
					bigCopy(x1, 0, nl, x1_temp, 0);
				}
				else{
					bigCopy(x1, 0, nl, sub_temp, 0);
				}
				addControl = 0;


			}
			else {
				// v <- v - u
				bigSub(sub_temp, &controlBit, v, al, u, nl);
				bigCopy(v, 0, nl, sub_temp, 0);

				//TODO: carry bit might overflow
				// x2 <- x2 - x1
				if(is_negative_x1 == 0 && is_negative_x2 == 0){ // Both Pozitive
					bigSub(sub_temp, &controlBit, x2, nl, x1,  nl);
					if(controlBit == 1){
						is_negative_x2 = 1;
						bigSub(sub_temp2, &controlBit, zero_num, nl, sub_temp, nl);
						bigCopy(sub_temp, 0, nl, sub_temp2, 0);
					}
				}
				else if(is_negative_x2 == 1 && is_negative_x1 == 0){;
					bigAdd(x2_temp, x2, nl, x1, nl);
					addControl = 1;
				}
				else if(is_negative_x2 == 0 && is_negative_x1 == 1){
					bigAdd(x2_temp, x2, nl, x1, nl);
					addControl = 1;
				}
				else{ // Both negative
					bigSub(sub_temp, &controlBit, x2, nl, x1, nl);
					if(controlBit == 1){
						is_negative_x2 = 0;
						bigSub(sub_temp2, &controlBit, zero_num, nl, sub_temp, nl);
						bigCopy(sub_temp, 0, nl, sub_temp2, 0);
					}
				}


				if(addControl == 1){
					bigCopy(x2, 0, nl, x2_temp, 0);
				}
				else{
					bigCopy(x2, 0, nl, sub_temp, 0);
				}
				addControl = 0;

			}

			bigIsEqualUi(&is_one_u, u, al, 1L);
			bigIsEqualUi(&is_one_v, v, nl, 1L);

		} // end of while

		//if u == 1
		if(u[0] == 1){
			ui x1_barret = new ui_t[2 * nl];
			bigCopy(x1_barret, 0, nl, x1, 0);
			for(i = nl; i < 2 * nl; i++)
				x1_barret[i] = 0;
			// -7 mod 73
			if(is_negative_x1 == 1){
				bigSub(z, &controlBit, n_temp, nl, x1, nl);
			}
			else
				barretReduction(z, x1_barret, 2 * nl, n_temp, nl, mu, mul);

		}
		else {
			ui x2_barret = new ui_t[2 * nl];
			bigCopy(x2_barret, 0, nl, x2, 0);
			for(i = nl; i < 2 * nl; i++)
				x2_barret[i] = 0;

			if(is_negative_x2 == 1){
				bigSub(z, &controlBit, n_temp, nl, x2, nl);
			}

			else
				barretReduction(z, x2_barret, 2 * nl, n_temp, nl, mu, mul);
		}
	}
}

void myset(long int* A, long int t, long int len){

	//long int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	//if (t_index < SIZE) {

		long int i;

		A[0] = t;

		for (i = 1; i < len; i++) {
			A[i] = 0;
		}

	//}
}

#define mysar(a) __asm__( \
	"sar $62, %%rax;" \
	: "=a"((a)) \
	: "a"((a)) \
)

void myshift(long int* A, long int len){

	//long int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	//if (t_index < SIZE) {

		long int i;

		mysar(A[0]);

		A[1] += A[0];

		for (i = 0; i < len-1;) {
			A[i] = A[++i];
		}

		A[len-1] = 0;
	//}

}



//void mysar(long int* a) {                       // a    = 11000..1
//
//	//long int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;
//
//	//if (t_index < SIZE) {
//
//		long int sign = a[0] >> (63);     				// sign = 1
//		a[0] >>= (62);
//
//		if(sign == 1) {
//			a[0] |= 0xFFFFFFFC;             				// a = a | 1111...1111
//		}
//	//}
//}

void mycopy(long int *destination, long int *source, long int length) {
	long int i;
	for (i = 0; i < length; ++i) {
		destination[i] = source[i];
	}
}


//Gerekli swap ve hesaplama islemleri
void divsteps2(long int n, long int t, long int* delta, long int f, long int g, long int* uu, long int* vv, long int* qq, long int* rr){

	//long int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	//if (t_index < SIZE) {

		long int d, u, v, q, r, tt, g0, i;

		u=1L; v=0L; q=0L; r=1L; d= *delta;

		for(i = 0L; i < n; i++) {
				if ((d > 0) && ((g >> i) & 1)) {
					d = -d;
					tt = -f;  f = g;  g = tt;
					tt = -u;  u = q;  q = tt;
					tt = -v;  v = r;  r = tt;
				}
				g0 = (g >> i) & 1;
				++d;
				g = g + g0*f;  f = 2*f;
				q = q + g0*u;  u = 2*u;
				r = r + g0*v;  v = 2*v;
			}
			*delta = d; *uu = u; *vv = v; *qq = q; *rr = r;
	//}
}

// uf = u * f, vg = v * g, zn = uf + vg + c

// zn fn ve gn multi-precision
void myMul(long int* zn, long int u, long int* fn, long int v, long int* gn, long int len){

	//long int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	//if (t_index < SIZE) {

		long int i, h, l, f, g, c = 0;

		for (i = 0; i < len; ++i) {
			f = fn[i]; g = gn[i];
			ufvg(u, f, v, g, c, h, l);
			zn[i] = l;
			h = h + (((signed long)h) < 0);
			c = h;
		}
		zn[len] = c;

	//}
}

//void ufvg(long int u, long int f, long int v, long int g, long int c, long int* h, long int* l){
//
//	long int rax, rsi, rdi, rdx, rcx, r10, r11, carryBit, raxBits;
//
//	__mul_lo(rcx, f, u); 					// "imulq %3;"
//	__mul_hi(rsi, f, u);					// "movq %%rax, %%rsi;"
//	rax = v;								// "movq %%rdx, %%rcx;"
//
//	__mul_lo(rdx, g, rax);					// "movq %4, %%rax;"
//	__mul_hi(rax, g, rax);					// "imulq %5;
//
//	__add_cc(rax, rsi, rax);				// "addq %%rsi, %%rax;"
//	__addc_cc(rdx, rcx, rdx);				// "adcq %%rcx, %%rdx;"
//
//	r10 = c;								// "movq %6, %%r10;"
//	carryBit = (r10 >> (W - 1));			// "btq $63, %%r10;"
//
//	r11 = 0;								// "movq $0, %%r11;"
//
//	__sub_cc(r11, r11, carryBit);			// "sbbq $0, %%r11;"
//
//	__add_cc(rax, rax, r10);				// "addq %%r10, %%rax;"
//	__addc_cc(rdx, rdx, r11); 				// "adcq %%r11, %%rdx;"
//
//	rdi = 0x3fffffff;						// "movq $0x3fffffffffffffff, %%rdi;"
//	r11 = 0xc0000000;						// "movq $0xc000000000000000, %%r11;"
//
//	rdx = r10;								// "movq %%rdx, %%r10;"
//	r10 &= r11;								// "andq %%r11, %%r10;"
//
//	raxBits = (rax >> (W - 2));
//	rdx = (rdx << 2);						// "shld $2, %%rax, %%rdx;"
//	rdx |= raxBits;
//
//	rax &= rdi;								// "andq %%rdi, %%rax;
//	rax |= r10;								// "orq %%r10, %%rax;"
//
//	*l = rax;
//	*h = rdx;
//
//}


// f = GCD(f, g), v = ModInv(g, f)
void safegcd(long int* z, long int* f, long int* g, long int* U, long int* V, long int* Q, long int* R, long int* precomp, long int len){

	//long int t_index = (blockDim.x * blockIdx.x) + threadIdx.x;

	//if (t_index < SIZE) {
		long int i, u, v, q, r, delta = 1, iter = 12, t = 30;

		long int *temp_f = (long int *)malloc(sizeof(long int) * len);
		long int *temp_g = (long int *)malloc(sizeof(long int) * len);

		long int *temp_U = (long int *)malloc(sizeof(long int) * len);
		long int *temp_Q = (long int *)malloc(sizeof(long int) * len);
		long int *temp_V = (long int *)malloc(sizeof(long int) * len);
		long int *temp_R = (long int *)malloc(sizeof(long int) * len);

		myset(U, 1, len);
		myset(V, 0, len);
		myset(Q, 0, len);
		myset(R, 1, len);

		for (i = 0; i < iter; i++) {
				divsteps2(t, t, &delta, f[0], g[0], &u, &v, &q, &r);

				myMul(temp_f, u, f, v, g, len);
				myMul(temp_g, q, f, r, g, len);
				mycopy(f, temp_f, len);
				mycopy(g, temp_g, len);

				myMul(temp_U, u, U, v, Q, len);
				myMul(temp_Q, q, U, r, Q, len);
				mycopy(U, temp_U, len);
				mycopy(Q, temp_Q, len);

				myMul(temp_V, u, V, v, R, len);
				myMul(temp_R, q, V, r, R, len);

				mycopy(V, temp_V, len);
				mycopy(R, temp_R, len);

				myshift(f, len);
				myshift(g, len);

			}
	//}
}


void myprint(FILE *file, const char *str, long int* a, long int len) {
	long int i;
	fprintf(file, "%s = ", str);
	for (i = 0; i < len; ++i) {
		fprintf(file, "+(%ld)*(2^62)^%lu", a[i], i);
	}
	fprintf(file, ";\n");
}

void trace() {
	long int len = 8;

	long int *f = (long int *)calloc(sizeof(long int), len);
	long int *g = (long int *)calloc(sizeof(long int), len);
	long int *z = (long int *)calloc(sizeof(long int), len);
	long int *U = (long int *)calloc(sizeof(long int), len);
	long int *V = (long int *)calloc(sizeof(long int), len);
	long int *Q = (long int *)calloc(sizeof(long int), len);
	long int *R = (long int *)calloc(sizeof(long int), len);
	long int *precomp = (long int *)calloc(sizeof(long int), len);

	// MAGMA precomp --> Intseq(25651855896792820957299331057107636096084522772917456926367662764254847716857, 2^62):Hex;
	precomp[0] = 0x1DC1855B1B224DF9UL;
	precomp[1] = 0x329511A76508B241UL;
	precomp[2] = 0x1639D9DB0CCD4719UL;
	precomp[3] = 0x2D9BE62C1DB593D6UL;
	precomp[4] = 0x0000000000000038UL;

	// MAGMA f --> Intseq(2^255-19, 2^62):Hex;
	f[0] = 0x3FFFFFFFFFFFFFEDUL;
	f[1] = 0x3FFFFFFFFFFFFFFFUL;
	f[2] = 0x3FFFFFFFFFFFFFFFUL;
	f[3] = 0x3FFFFFFFFFFFFFFFUL;
	f[4] = 0x000000000000007FUL;

	// MAGMA g --> Intseq(6589194097717408581123106985632999088899568769370928766986301338197408648914,2^62):Hex;
	g[0] = 0x007724FA7FB48ED2UL;
	g[1] = 0x17BCDC94E85E5CA8UL;
	g[2] = 0x253A62FAB2AAAD3BUL;
	g[3] = 0x2456823035855B69UL;
	g[4] = 0x000000000000000EUL;

	/*g[0] = 2878036360140264451;
	g[1] = 2169902097010195153;
	g[2] = 2175766470785808513;
	g[3] = 2076646687522013864;
	g[4] = 90;*/

	safegcd(z, f, g, U, V, Q, R, precomp, len);

	myprint(stdout, "tf", f, len);
	myprint(stdout, "tg", g, len);
	myprint(stdout, "tU", U, len);
	myprint(stdout, "tV", V, len);
	myprint(stdout, "tQ", Q, len);
	myprint(stdout, "tR", R, len);



	free(precomp);
	free(R);
	free(Q);
	free(V);
	free(U);
	free(z);
	free(g);
	free(f);
}

