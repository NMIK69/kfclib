#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "matrix.h"
#include "la_arena.h"


static enum matrix_allocator_type mat = MATRIX_ALLOC_MALLOC;
static struct la_arena *arena = NULL;

void matrix_set_allocator(enum matrix_allocator_type t, void *allocator)
{
	mat = t;
	if(t == MATRIX_ALLOC_LA_ARENA) {
		assert(allocator != NULL);
		arena = (struct la_arena *)allocator;
	}
}

struct matrix *matrix_alloc(size_t H, size_t W)
{
	struct matrix *res = NULL;
	size_t m_size = H * W * sizeof(*(res->m));
	size_t total_size = sizeof(*res) + m_size; 


	switch(mat) {
	case MATRIX_ALLOC_MALLOC:
		res = malloc(total_size);
		break;
	case MATRIX_ALLOC_LA_ARENA:
		assert(arena != NULL);
		res = la_arena_alloc(arena, total_size);	
		break;
	default:
		assert(0);
	}
	
	if(res == NULL)
		return NULL;	
	
	res->H = H;
	res->W = W;
	
	matrix_fill(res, 0.0);

	return res;
}

void matrix_free(struct matrix *M)
{
	if(M == NULL)
		return;

	if(mat == MATRIX_ALLOC_MALLOC)
		free(M);	
}

void matrix_print(struct matrix *M)
{
	for(size_t r = 0; r < M->H; r++) {
		for(size_t c = 0; c < M->W; c++) {
			float val = MATRIX_AT(M, r, c);
			printf("%+9.15lf  ", val);
		}
		printf("\n");
	}
}

void matrix_cpy(struct matrix *dst, struct matrix *src)
{
	assert(dst->H == src->H);
	assert(dst->W == src->W);

	for(size_t i = 0; i < dst->H; i++) {
		for(size_t j = 0; j < dst->W; j++) {
			MATRIX_AT(dst, i, j) = MATRIX_AT(src, i, j);
		}
	}
}

void matrix_fill(struct matrix *M, float val)
{
	for(size_t r = 0; r < M->H; r++) {
		for(size_t c = 0; c < M->W; c++) {
			MATRIX_SET(M, r, c, val);
		}
	}
}

float matrix_det3x3(struct matrix *M)
{
/* ref: https://en.wikipedia.org/wiki/Rule_of_Sarrus */

	assert(M->H == 3);
	assert(M->W == 3);

	float a = MATRIX_AT(M, 0, 0);
	float b = MATRIX_AT(M, 0, 1);
	float c = MATRIX_AT(M, 0, 2);
	float d = MATRIX_AT(M, 1, 0);
	float e = MATRIX_AT(M, 1, 1);
	float f = MATRIX_AT(M, 1, 2);
	float g = MATRIX_AT(M, 2, 0);
	float h = MATRIX_AT(M, 2, 1);
	float i = MATRIX_AT(M, 2, 2);

	return ((a*e*i) + (b*f*g) + (c*d*h) - (g*e*c) - (h*f*a) - (i*d*b));
}

struct matrix *matrix_inv3x3(struct matrix *M)
{

/* ref: https://www.cuemath.com/algebra/inverse-of-3x3-matrix/ 
	https://stackoverflow.com/a/984054 */

	assert(M->H == 3);
	assert(M->W == 3);

	struct matrix *inv = matrix_alloc(3, 3);
	if(inv == NULL)
		return NULL;

	float a = MATRIX_AT(M, 0, 0);
	float b = MATRIX_AT(M, 0, 1);
	float c = MATRIX_AT(M, 0, 2);
	float d = MATRIX_AT(M, 1, 0);
	float e = MATRIX_AT(M, 1, 1);
	float f = MATRIX_AT(M, 1, 2);
	float g = MATRIX_AT(M, 2, 0);
	float h = MATRIX_AT(M, 2, 1);
	float i = MATRIX_AT(M, 2, 2);

	float det_M = ((a*e*i) + (b*f*g) + (c*d*h) - (g*e*c) - (h*f*a) - (i*d*b));

	float det_a = (e*i - f*h) / det_M;
	float det_b = (c*h - b*i) / det_M;
	float det_c = (b*f - c*e) / det_M;
	float det_d = (f*g - d*i) / det_M;
	float det_e = (a*i - c*g) / det_M; 
	float det_f = (c*d - a*f) / det_M; 
	float det_g = (d*h - e*g) / det_M; 
	float det_h = (b*g - a*h) / det_M; 
	float det_i = (a*e - b*d) / det_M; 

	MATRIX_SET(inv, 0, 0, det_a);	
	MATRIX_SET(inv, 0, 1, det_b);	
	MATRIX_SET(inv, 0, 2, det_c);	
	MATRIX_SET(inv, 1, 0, det_d);	
	MATRIX_SET(inv, 1, 1, det_e);	
	MATRIX_SET(inv, 1, 2, det_f);	
	MATRIX_SET(inv, 2, 0, det_g);	
	MATRIX_SET(inv, 2, 1, det_h);	
	MATRIX_SET(inv, 2, 2, det_i);	

	return inv;
}

struct matrix *matrix_mmul(struct matrix *M1, struct matrix *M2)
{
	assert(M1->W == M2->H);
	
	struct matrix *res = matrix_alloc(M1->H, M2->W);
	if(res == NULL)
		return NULL;

	float sum = 0.0f;
	float m1, m2;

	for(size_t i = 0; i < M2->W; i++) {
		for(size_t j = 0; j < M1->H; j++) {
			sum = 0.0f;
			for(size_t k = 0; k < M1->W; k++) {
				m1 = MATRIX_AT(M1, j, k);
				m2 = MATRIX_AT(M2, k, i);
				sum += m1 * m2;
			}
			MATRIX_SET(res, j, i, sum);
		}
	}

	return res;
}

void matrix_fill_diag(struct matrix *M, float val)
{
	assert(M->W == M->H);

	for(size_t i = 0; i < M->H; i++) {
		MATRIX_SET(M, i, i, val);
	}
}

struct matrix *matrix_transpose(struct matrix *M)
{
	struct matrix *res = matrix_alloc(M->W, M->H);
	if(res == NULL)
		return NULL;
	

	float tmp;
	for(size_t i = 0; i < res->H; i++) {
		for(size_t j = 0; j < res->W; j++) {
			tmp = MATRIX_AT(M, j, i);
			MATRIX_SET(res, i, j, tmp);
		}
	}

	return res;
}

void matrix_smul(struct matrix *M, float s)
{
	float tmp;
	for(size_t i = 0; i < M->H; i++) {
		for(size_t j = 0; j < M->W; j++) {
			tmp = MATRIX_AT(M, i, j);
			tmp = tmp * s;
			MATRIX_SET(M, i, j, tmp);
		}
	}
}

struct matrix *matrix_madd(struct matrix *M1, struct matrix *M2)
{
	assert(M1->H == M2->H);
	assert(M1->W == M2->W);

	struct matrix *res = matrix_alloc(M1->H, M1->W);
	if(res == NULL)
		return NULL;


	float sum, m1, m2;
	for(size_t i = 0; i < M1->H; i++) {
		for(size_t j = 0; j < M1->W; j++) {
			m1 = MATRIX_AT(M1, i, j);
			m2 = MATRIX_AT(M2, i, j);
			sum = m1 + m2;
			MATRIX_SET(res, i, j, sum);
		}
	}
	
	return res;
}

struct matrix *matrix_msub(struct matrix *M1, struct matrix *M2)
{
	assert(M1->H == M2->H);
	assert(M1->W == M2->W);

	struct matrix *res = matrix_alloc(M1->H, M1->W);
	if(res == NULL)
		return NULL;


	float sum, m1, m2;
	for(size_t i = 0; i < M1->H; i++) {
		for(size_t j = 0; j < M1->W; j++) {
			m1 = MATRIX_AT(M1, i, j);
			m2 = MATRIX_AT(M2, i, j);
			sum = m1 - m2;
			MATRIX_SET(res, i, j, sum);
		}
	}
	
	return res;
}
