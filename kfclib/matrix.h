#ifndef MATRIX_H
#define MATRIX_H

#define MATRIX_IDX(M, R, C) (((R) * (M)->W) + (C))
#define MATRIX_AT(M, R, C) ((M)->m[MATRIX_IDX((M), (R), (C))])
#define MATRIX_SET(M, R, C, V) ((M)->m[MATRIX_IDX((M), (R), (C))] = (V))


struct matrix
{
	size_t H;
	size_t W;

	float m[];
};

enum matrix_allocator_type
{
	MATRIX_ALLOC_MALLOC = 0,
	MATRIX_ALLOC_LA_ARENA,
};

void matrix_set_allocator(enum matrix_allocator_type t, void *allocator);
struct matrix *matrix_alloc(size_t H, size_t W);
void matrix_free(struct matrix *M);

void matrix_cpy(struct matrix *dst, struct matrix *src);

void matrix_fill(struct matrix *M, float val);
void matrix_fill_diag(struct matrix *M, float val);
void matrix_print(struct matrix *M);

float matrix_det3x3(struct matrix *M);
float matrix_det2x2(struct matrix *M);
float matrix_det(struct matrix *M);
struct matrix *matrix_make_minor(struct matrix *M, 
		size_t row, size_t col, int in_place,
		struct matrix *M_minor);
struct matrix *make_cofactor(struct matrix *M);
struct matrix *matrix_make_adjunct(struct matrix *M);
struct matrix *matrix_inv_ana(struct matrix *M);

int matrix_get_alloc_free_cnt(void);

struct matrix *matrix_inv3x3(struct matrix *M);
struct matrix *matrix_mmul(struct matrix *M1, struct matrix *M2);
struct matrix *matrix_transpose(struct matrix *M);
struct matrix *matrix_madd(struct matrix *M1, struct matrix *M2);
struct matrix *matrix_msub(struct matrix *M1, struct matrix *M2);
void matrix_smul(struct matrix *M, float s);

struct matrix *matrix_make_I(size_t N);
struct matrix *matrix_inv_gj(struct matrix *M);

#endif
