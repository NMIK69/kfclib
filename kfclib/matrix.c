#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "matrix.h"
#include "la_arena.h"


static enum matrix_allocator_type mat = MATRIX_ALLOC_MALLOC;
static struct la_arena *arena = NULL;
static int matrix_alloc_cnt = 0;
static int matrix_free_cnt = 0;

int matrix_get_alloc_free_cnt(void)
{
	return matrix_alloc_cnt - matrix_free_cnt;
}

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

	matrix_alloc_cnt += 1;

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
	
	matrix_free_cnt += 1;
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

struct matrix *matrix_make_I(size_t N)
{
	struct matrix *I = matrix_alloc(N,N);
	if(I == NULL) {
		return NULL;
	}

	matrix_fill(I, 0.0f);

	for(size_t i = 0; i < N; i++) {
		MATRIX_AT(I, i, i) = 1.0f;
	}
	return I;
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

float matrix_det2x2(struct matrix *M)
{
	assert(M->H == 2);
	assert(M->W == 2);

	float a = MATRIX_AT(M, 0, 0);
	float b = MATRIX_AT(M, 0, 1);
	float c = MATRIX_AT(M, 0, 2);
	float d = MATRIX_AT(M, 1, 0);
	
	return ((a*d) - (b*c));
}

/* 3x3 matrix inverse with analytical method */
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


float matrix_det(struct matrix *M)
{
	float det = 0.0f;
	
	assert(M->W == M->H);

	if(M->W == 1)
		return (MATRIX_AT(M, 0, 0));
	else if(M->W == 2)
		return matrix_det2x2(M);
	else if(M->W == 3)
		return matrix_det3x3(M);

	struct matrix *minor = matrix_alloc(M->W-1,  M->H-1);
	assert(minor != NULL);

	
	for(size_t i = 0; i < M->W; i++) {
		matrix_make_minor(M, 0, i, 1, minor);
		assert(minor != NULL);
		float val = (MATRIX_AT(M, 0, i)) * matrix_det(minor);

		det += (i % 2 == 0) ? val : -val;
	}

	matrix_free(minor);
	return det;
}

struct matrix *matrix_make_minor(struct matrix *M, 
		size_t row, size_t col, int in_place,
		struct matrix *M_minor)
{
	int r = 0;
	int c = 0;
	struct matrix *minor;

	if(in_place == 0) {
		minor = matrix_alloc(M->H, M->W);
		if(minor == NULL)
			return NULL;
	}
	else {
		minor = M_minor;	
	}
	

	for(size_t i = 0; i < M->H; i++) {
		c = 0;
		if(i != row) {
			for(size_t j = 0; j < M->W; j++) {
				if(j != col) {
					MATRIX_AT(minor, r, c) = MATRIX_AT(M, i, j);
					c = c + 1; 
				}
			}
			r = r + 1;
		}
	}

	return minor;
}

struct matrix *make_cofactor(struct matrix *M)
{
	struct matrix *co;
	struct matrix *minor;
	float det_minor = 0.0f;

	co = matrix_alloc(M->H, M->W);
	if(co == NULL) {
		return NULL;
	}

	minor = matrix_alloc(M->H-1, M->W-1);
	if(minor == NULL) {
		matrix_free(co);
		return NULL;
	}

	for(size_t i = 0; i < M->H; i++) {
		for(size_t j = 0; j < M->W; j++) {
			matrix_make_minor(M, i, j, 1, minor);			
			assert(minor != NULL);
			det_minor = matrix_det(minor);
			if((i + j) % 2 == 0)
				MATRIX_AT(co, i, j) = det_minor;
			else 
				MATRIX_AT(co, i, j) = -det_minor;
		}
	}

	matrix_free(minor);
	return co;
}

struct matrix *matrix_make_adjunct(struct matrix *M)
{
	struct matrix *co;
	struct matrix *adj;

	co = make_cofactor(M);
	if(co == NULL)
		return NULL;

	adj = matrix_transpose(co);
	if(adj == NULL) {
		matrix_free(co);
		return NULL;
	}

	matrix_free(co);
	return adj;
}

/* general matrix inverse with analytical method */
/* takes up a large amount of memory compared to 
   inverse with gauss jordan elimination. */
struct matrix *matrix_inv_ana(struct matrix *M)
{
	float det;
	struct matrix *adj;
	struct matrix *inv;

	inv = matrix_alloc(M->H, M->W);
	if(inv == NULL) {
		return NULL;
	}

	det = matrix_det(M);
	adj = matrix_make_adjunct(M);
	if(adj == NULL) {
		matrix_free(inv);
		return NULL;
	}

	for(size_t i = 0; i < adj->H; i++) {
		for(size_t j = 0; j < adj->W; j++) {
			MATRIX_AT(inv, i, j) = (MATRIX_AT(adj, i, j)) / det;
		}
	}

	matrix_free(adj);
	return inv;
}

static int find_pivot(struct matrix *M, size_t row, size_t col)
{
	float tmp;
	float max_num = 0.0f;
	int max_idx = -1;

	while(row < M->H) {
		tmp = MATRIX_AT(M, row, col);

		if(fabs(tmp) > fabs(max_num)) {
			max_idx = (int)row;
			max_num = tmp;
		}
		row += 1;
	}

	if(max_idx == -1)
		return 0;
	return max_idx;
}

static void swap_rows(struct matrix *M, size_t r1, size_t r2)
{
	assert(M->H > r1 && M->H > r2);
	
	float tmp;
	for(size_t i = 0; i < M->W; i++) {
		tmp = MATRIX_AT(M, r1, i);	
		MATRIX_AT(M, r1, i) = MATRIX_AT(M, r2, i);
		MATRIX_AT(M, r2, i) = tmp;
	}
}

static void row_norm(struct matrix *M, struct matrix *I, size_t row, size_t col)
{
	float s = 1.0f/MATRIX_AT(M, row, col);
	for(size_t i = 0; i < M->W; i++) {
		MATRIX_AT(M, row, i) = s * MATRIX_AT(M, row, i); 
		MATRIX_AT(I, row, i) = s * MATRIX_AT(I, row, i); 
	}
}

static void row_reduce(struct matrix *M, struct matrix *I,
	size_t src, size_t dst, size_t col)
{
	assert(M->H > src);
	assert(M->H > dst);

	float tmp_M, tmp_I;
	float s_src = MATRIX_AT(M, dst, col);

	/* reduce dst row */
	for(size_t i = 0; i < M->W; i++) {
		tmp_M = MATRIX_AT(M, dst, i) - (s_src * MATRIX_AT(M,src,i));
		tmp_I = MATRIX_AT(I, dst, i) - (s_src * MATRIX_AT(I,src,i));

		MATRIX_AT(M, dst, i) = tmp_M;
		MATRIX_AT(I, dst, i) = tmp_I;
	}
}

/* general matrix inversion with gauss jordan elimination (with partial pivoting) */
struct matrix *matrix_inv_gj(struct matrix *M)
{
	assert(M->W == M->H);
	size_t N = M->H;
	struct matrix *inv;
	struct matrix *I;

	inv = matrix_alloc(N,N);
	if(inv == NULL) return NULL;
	matrix_cpy(inv, M);

	I = matrix_make_I(N);
	if(I == NULL) return NULL;
	
	int pivot_idx;
	for(size_t i = 0; i < N; i++) {
		/* find largest pivot and swap*/
		pivot_idx = find_pivot(inv, i, i);
		assert(pivot_idx != -1);
		swap_rows(inv, pivot_idx, i);
		swap_rows(I, pivot_idx, i);

		/* norm pivot row*/
		row_norm(inv, I, i, i);

		/* row reduce */
		for(size_t j = 0; j < N; j++) {
			if(j != i) {
				row_reduce(inv, I, i, j, i);
			}
		}
	}

	matrix_free(inv);

	return I;
}
