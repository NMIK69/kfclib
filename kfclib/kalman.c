#include <stdlib.h>
#include <stddef.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>

#include "kalman.h"
#include "matrix.h"
#include "la_arena.h"
#include "quaternion.h"


#define KF_NED_AX_REF  0.0f
#define KF_NED_AY_REF  0.0f
#define KF_NED_AZ_REF  1.0f

#define KF_ENU_AX_REF  0.0f
#define KF_ENU_AY_REF  0.0f
#define KF_ENU_AZ_REF -1.0f


struct vec3f
{
	float x, y, z;
};

static struct kalman_filter *kf_alloc(void);
static void kf_norm(struct kalman_filter *kf);
static struct matrix *make_F(float dt, float wx, float wy, float wz);
static struct matrix *make_Q(struct kalman_filter *kf);
static struct matrix *make_W(struct kalman_filter *kf);
static struct matrix *make_H(struct kalman_filter *kf);
static struct vec3f norm_accel(float ax, float ay, float az);

static int kf_update(struct kalman_filter *kf, float ax, float ay, float az);
static int kf_predict(struct kalman_filter *kf, float wx, float wy, float wz);


static struct la_arena *matrix_arena;
static struct matrix *I4;


struct kalman_filter *kf_init(float dt, float var_a, float var_w, float var_P)
{
	struct kalman_filter *kf = kf_alloc();
	if(kf == NULL)
		return NULL;

	kf->dt = dt;
	kf->var_a = var_a;
	kf->var_w = var_w;
	kf->var_P = var_P;

	/* use NED as default reference frame. */
	kf->ax_ref = KF_NED_AX_REF;
	kf->ay_ref = KF_NED_AY_REF; 
	kf->az_ref = KF_NED_AZ_REF; 

	/* use identity quaternion as default */
	MATRIX_SET(kf->x, 0, 0, 1.0);
	MATRIX_SET(kf->x, 1, 0, 0.0);
	MATRIX_SET(kf->x, 2, 0, 0.0);
	MATRIX_SET(kf->x, 3, 0, 0.0);

	matrix_fill_diag(kf->P, kf->var_P);
	matrix_fill_diag(kf->R, kf->var_a);
	matrix_fill_diag(I4, 1.0);

	return kf;
}

/* expecets normed quaternion. */
void kf_set_q(struct kalman_filter *kf, float qw, float qx, float qy, float qz)
{
	MATRIX_SET(kf->x, 0, 0, qw);
	MATRIX_SET(kf->x, 1, 0, qx);
	MATRIX_SET(kf->x, 2, 0, qy);
	MATRIX_SET(kf->x, 3, 0, qz);
}

void kf_set_aref(struct kalman_filter *kf, float ax, float ay, float az)
{
	struct vec3f anorm = norm_accel(ax, ay, az);	

	kf->ax_ref = anorm.x;
	kf->ay_ref = anorm.y; 
	kf->az_ref = anorm.z; 
}


void kf_free(struct kalman_filter *kf)
{
	if(kf == NULL)
		return;

	matrix_set_allocator(MATRIX_ALLOC_MALLOC, NULL);

	matrix_free(kf->x);
	matrix_free(kf->P);
	matrix_free(kf->R);
	matrix_free(kf->K);
	matrix_free(I4);

	la_arena_free(matrix_arena);
}



int kf_filt(struct kalman_filter *kf, 
	   float wx, float wy, float wz,
	   float ax, float ay, float az)
{
	struct vec3f anorm = norm_accel(ax, ay, az);	

	int err = kf_predict(kf, wx, wy, wz);
	if(err != 0)
		return err;
	
	err = kf_update(kf, anorm.x, anorm.y, anorm.z);
	if(err != 0)
		return err;

	/* update state quaternion */
	kf->q.w = MATRIX_AT(kf->x, 0, 0);
	kf->q.x = MATRIX_AT(kf->x, 1, 0);
	kf->q.y = MATRIX_AT(kf->x, 2, 0);
	kf->q.z = MATRIX_AT(kf->x, 3, 0);

	return 0;
}


static int kf_predict(struct kalman_filter *kf, float wx, float wy, float wz)
{
	/* set matrix allocator to linear arena allocator for more efficient and
	 * faster temporary matrix allocations. */
	matrix_set_allocator(MATRIX_ALLOC_LA_ARENA, matrix_arena);

	struct matrix *Q = make_Q(kf);
	if(Q == NULL) 
		goto err_out;
	
	struct matrix *F = make_F(kf->dt, wx, wy, wz);
	if(F == NULL) 
		goto err_out;

	struct matrix *prd_x = matrix_mmul(F, kf->x);
	if(prd_x == NULL) 
		goto err_out;
	
	struct matrix *prd_P = matrix_madd(
			matrix_mmul(
				matrix_mmul(F, kf->P), 
				matrix_transpose(F)
				), 
			Q);

	if(prd_P == NULL)
		goto err_out;

	matrix_cpy(kf->P, prd_P);
	matrix_cpy(kf->x, prd_x);
	kf_norm(kf);

	/* reset arena. */
	la_arena_reset(matrix_arena);
	return 0;

err_out:
	la_arena_reset(matrix_arena);
	return -1;
}

static int kf_update(struct kalman_filter *kf, float ax, float ay, float az)
{
	/* set matrix allocator to linear arena allocator for more efficient and
	 * faster temporary matrix allocations. */
	matrix_set_allocator(MATRIX_ALLOC_LA_ARENA, matrix_arena);

	struct matrix *H = make_H(kf);
	if(H == NULL)
		goto err_out;
	
	struct matrix *z = matrix_alloc(3, 1);
	MATRIX_SET(z, 0, 0, ax);
	MATRIX_SET(z, 1, 0, ay);
	MATRIX_SET(z, 2, 0, az);

	struct matrix *v = matrix_msub(z, matrix_mmul(H, kf->x));
	if(v == NULL)
		goto err_out;

	struct matrix *S = matrix_madd(
				matrix_mmul(
					matrix_mmul(H, kf->P), 
					matrix_transpose(H)),
				kf->R);
	if(S == NULL)
		goto err_out;
	
	struct matrix *K = matrix_mmul( 
				matrix_mmul(kf->P, matrix_transpose(H)), 
				matrix_inv3x3(S)
				);
	if(K == NULL)
		goto err_out;

	struct matrix *est_x = matrix_madd(kf->x, matrix_mmul(K, v));
	if(est_x == NULL)
		goto err_out;

	struct matrix *est_P = matrix_mmul(
					matrix_msub(I4, matrix_mmul(K, H)),
				kf->P);
	if(est_P == NULL)
		goto err_out;
	


	matrix_cpy(kf->x, est_x);
	matrix_cpy(kf->P, est_P);
	matrix_cpy(kf->K, K);

	/* norm state quaternion */
	kf_norm(kf);

	la_arena_reset(matrix_arena);
	return 0;

err_out:
	la_arena_reset(matrix_arena);
	return -1;
}

static void kf_norm(struct kalman_filter *kf)
{
	struct quaternion q;
	q.w = MATRIX_AT(kf->x, 0, 0);
	q.x = MATRIX_AT(kf->x, 1, 0);
	q.y = MATRIX_AT(kf->x, 2, 0);
	q.z = MATRIX_AT(kf->x, 3, 0);

	struct quaternion n = quat_norm(q); 
	MATRIX_SET(kf->x, 0, 0, n.w);
	MATRIX_SET(kf->x, 1, 0, n.x);
	MATRIX_SET(kf->x, 2, 0, n.y);
	MATRIX_SET(kf->x, 3, 0, n.z);
}

static struct kalman_filter *kf_alloc(void)
{
	struct kalman_filter *kf = malloc(sizeof(*kf));
	if(kf == NULL)
		return NULL;

	/* create linear arena allocator that is used for efficient and fast 
	   temporary matrix allocations */
	matrix_arena = la_arena_create(sizeof(float) * 500);
	if(matrix_arena == NULL)
		goto err_out;

	/* set the matrix allocator 'malloc', so that the following matrecies
	 * will be allocated on the heap. */
	matrix_set_allocator(MATRIX_ALLOC_MALLOC, NULL);

	kf->x = matrix_alloc(4, 1);
	kf->P = matrix_alloc(4, 4);
	kf->R = matrix_alloc(3, 3);
	kf->K = matrix_alloc(4, 3);

	/* used as a constant */
	I4 = matrix_alloc(4, 4);

	if(kf->x == NULL || kf->P == NULL || kf->R == NULL || kf->K == NULL || I4 == NULL)
		goto err_out;

	return kf;

err_out:
	kf_free(kf);
	return NULL;
}

static struct matrix *make_H(struct kalman_filter *kf)
{
	/* See README.md for more details. */

	struct matrix *H = matrix_alloc(3, 4);
	if(H == NULL)
		return NULL;

	struct quaternion q;
	q.w = MATRIX_AT(kf->x, 0, 0);
	q.x = MATRIX_AT(kf->x, 1, 0);
	q.y = MATRIX_AT(kf->x, 2, 0);
	q.z = MATRIX_AT(kf->x, 3, 0);

	struct quaternion g;
	g.w = 0.0;
	g.x = kf->ax_ref; 
	g.y = kf->ay_ref; 
	g.z = kf->az_ref;	

	struct quaternion p = quat_mul(q, g);
	
	MATRIX_SET(H, 0, 0, p.x);
	MATRIX_SET(H, 0, 1, -p.w);
	MATRIX_SET(H, 0, 2, p.z);
	MATRIX_SET(H, 0, 3, -p.y);

	MATRIX_SET(H, 1, 0, p.y);
	MATRIX_SET(H, 1, 1, -p.z);
	MATRIX_SET(H, 1, 2, -p.w);
	MATRIX_SET(H, 1, 3, p.x);

	MATRIX_SET(H, 2, 0, p.z);
	MATRIX_SET(H, 2, 1, p.y);
	MATRIX_SET(H, 2, 2, -p.x);
	MATRIX_SET(H, 2, 3, -p.w);

	return H;
}


static struct matrix *make_F(float dt, float wx, float wy, float wz)
{
	/* See README.md for more details. */

	struct matrix *F = matrix_alloc(4, 4);
	if(F == NULL)
		return NULL;

	float s = dt/2.0;

	MATRIX_SET(F, 0, 0, 1.0); 
	MATRIX_SET(F, 0, 1, -s*wx); 
	MATRIX_SET(F, 0, 2, -s*wy); 
	MATRIX_SET(F, 0, 3, -s*wz); 

	MATRIX_SET(F, 1, 0, s*wx); 
	MATRIX_SET(F, 1, 1, 1.0); 
	MATRIX_SET(F, 1, 2, s*wz); 
	MATRIX_SET(F, 1, 3, -s*wy); 

	MATRIX_SET(F, 2, 0, s*wy); 
	MATRIX_SET(F, 2, 1, -s*wz); 
	MATRIX_SET(F, 2, 2, 1.0); 
	MATRIX_SET(F, 2, 3, s*wx); 

	MATRIX_SET(F, 3, 0, s*wz); 
	MATRIX_SET(F, 3, 1, s*wy); 
	MATRIX_SET(F, 3, 2, -s*wx); 
	MATRIX_SET(F, 3, 3, 1.0); 

	return F;
}

static struct matrix *make_Q(struct kalman_filter *kf)
{
	/* See README.md for more details. */

	struct matrix *W = make_W(kf);
	if(W == NULL)
		return NULL;

	struct matrix *Q = matrix_mmul(W, matrix_transpose(W));
	if(Q == NULL)
		return NULL;

	matrix_smul(Q, kf->var_w);

	return Q;
}


static struct matrix *make_W(struct kalman_filter *kf)
{
	/* See README.md for more details. */

	struct matrix *W = matrix_alloc(4, 3);
	if(W == NULL)
		return NULL;

	float qw = MATRIX_AT(kf->x, 0, 0); 
	float qx = MATRIX_AT(kf->x, 1, 0); 
	float qy = MATRIX_AT(kf->x, 2, 0); 
	float qz = MATRIX_AT(kf->x, 3, 0); 
	
	
	MATRIX_SET(W, 0, 0, -qx);	
	MATRIX_SET(W, 0, 1, -qy);	
	MATRIX_SET(W, 0, 2, -qz);	

	MATRIX_SET(W, 1, 0, qw);	
	MATRIX_SET(W, 1, 1, -qz);	
	MATRIX_SET(W, 1, 2, qy);	

	MATRIX_SET(W, 2, 0, qz);	
	MATRIX_SET(W, 2, 1, qw);	
	MATRIX_SET(W, 2, 2, -qx);	

	MATRIX_SET(W, 3, 0, -qy);	
	MATRIX_SET(W, 3, 1, qx);	
	MATRIX_SET(W, 3, 2, qw);	

	float s = kf->dt/2.0;
	matrix_smul(W, s);
	
	return W;
}

static struct vec3f norm_accel(float ax, float ay, float az)
{
	float x2 = ax * ax;
	float y2 = ay * ay;
	float z2 = az * az;

	float n = sqrt((x2 + y2 + z2));

	struct vec3f res = {0};
	res.x = ax / n;
	res.y = ay / n;
	res.z = az / n;
	
	return res;
}
