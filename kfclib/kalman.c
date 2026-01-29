#include <stdlib.h>
#include <stddef.h>
#include <assert.h>
#include <stdio.h>
#include <stdbool.h>
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

static struct kalman_filter *kf_alloc(size_t nx, size_t nz);
static void kf_norm(struct kalman_filter *kf);
static struct matrix *make_F(float dt, float wx, float wy, float wz);
static struct matrix *make_Q(struct kalman_filter *kf);
static struct matrix *make_W(struct kalman_filter *kf);
static struct matrix *make_H(struct kalman_filter *kf);

static int kf_predict(struct kalman_filter *kf, float wx, float wy, float wz);

static int kf_update_6dof(struct kalman_filter *kf, float ax, float ay, float az);

static int kf_update_9dof(struct kalman_filter *kf, 
			float ax, float ay, float az,
			float mx, float my, float mz);

static struct la_arena *matrix_arena;
static struct matrix *I4;


struct kalman_filter *kf_init(float dt, float var_a, float var_w,
				float var_m, float var_P, int type)
{
	struct kalman_filter *kf;

	if(type == KF_6DOF) 
		kf = kf_alloc(4, 3);
	else
		kf = kf_alloc(4, 6);

	if(kf == NULL)
		return NULL;

	kf->is_6dof = type == KF_6DOF ? true : false;
	kf->dt = dt;
	kf->var_a = var_a;
	kf->var_w = var_w;
	kf->var_P = var_P;

	/* use NED as default reference frame. */
	kf->ax_ref = KF_NED_AX_REF;
	kf->ay_ref = KF_NED_AY_REF; 
	kf->az_ref = KF_NED_AZ_REF; 

	kf->mx_ref = 0.0f; 
	kf->my_ref = 0.0f; 
	kf->mz_ref = 0.0f; 

	/* use identity quaternion as default */
	MATRIX_SET(kf->x, 0, 0, 1.0);
	MATRIX_SET(kf->x, 1, 0, 0.0);
	MATRIX_SET(kf->x, 2, 0, 0.0);
	MATRIX_SET(kf->x, 3, 0, 0.0);

	matrix_fill_diag(kf->P, var_P);
	
	MATRIX_SET(kf->R, 0, 0, var_a);
	MATRIX_SET(kf->R, 1, 1, var_a);
	MATRIX_SET(kf->R, 2, 2, var_a);

	if(type == KF_9DOF) {
		MATRIX_SET(kf->R, 3, 3, var_m);
		MATRIX_SET(kf->R, 4, 4, var_m);
		MATRIX_SET(kf->R, 5, 5, var_m);
	}

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
	kf->ax_ref = ax;
	kf->ay_ref = ay; 
	kf->az_ref = az; 
}

void kf_set_mref(struct kalman_filter *kf, float mx, float my, float mz)
{
	kf->mx_ref = mx;
	kf->my_ref = my; 
	kf->mz_ref = mz; 
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

int kf_filt_9dof(struct kalman_filter *kf, 
	   float wx, float wy, float wz,
	   float ax, float ay, float az,
	   float mx, float my, float mz)
{
	int err = kf_predict(kf, wx, wy, wz);
	if(err != 0)
		return err;
	

	err = kf_update_9dof(kf, ax, ay, az, 
				mx, my, mz);
	if(err != 0)
		return err;

	/* update state quaternion */
	kf->q.w = MATRIX_AT(kf->x, 0, 0);
	kf->q.x = MATRIX_AT(kf->x, 1, 0);
	kf->q.y = MATRIX_AT(kf->x, 2, 0);
	kf->q.z = MATRIX_AT(kf->x, 3, 0);

	return 0;
}

int kf_filt_6dof(struct kalman_filter *kf, 
	   float wx, float wy, float wz,
	   float ax, float ay, float az)
{
	int err = kf_predict(kf, wx, wy, wz);
	if(err != 0)
		return err;
	
	err = kf_update_6dof(kf, ax, ay, az);
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

static int kf_update_6dof(struct kalman_filter *kf, float ax, float ay, float az)
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
				matrix_inv_gj(S)
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


static int kf_update_9dof(struct kalman_filter *kf, 
			float ax, float ay, float az,
			float mx, float my, float mz)
{
	/* set matrix allocator to linear arena allocator for more efficient and
	 * faster temporary matrix allocations. */
	matrix_set_allocator(MATRIX_ALLOC_LA_ARENA, matrix_arena);

	struct matrix *H = make_H(kf);
	if(H == NULL)
		goto err_out;
	
	struct matrix *z = matrix_alloc(6, 1);
	MATRIX_SET(z, 0, 0, ax);
	MATRIX_SET(z, 1, 0, ay);
	MATRIX_SET(z, 2, 0, az);
	MATRIX_SET(z, 3, 0, mx);
	MATRIX_SET(z, 4, 0, my);
	MATRIX_SET(z, 5, 0, mz);

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
				matrix_inv_ana(S)
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

static struct kalman_filter *kf_alloc(size_t nx, size_t nz)
{
	struct kalman_filter *kf = malloc(sizeof(*kf));
	if(kf == NULL)
		return NULL;

	/* create linear arena allocator that is used for efficient and fast 
	   temporary matrix allocations */
	matrix_arena = la_arena_create(sizeof(float) * 9000);
	if(matrix_arena == NULL)
		goto err_out;

	/* set the matrix allocator 'malloc', so that the following matrecies
	 * will be allocated on the heap. */
	matrix_set_allocator(MATRIX_ALLOC_MALLOC, NULL);

	kf->x = matrix_alloc(nx, 1);
	kf->P = matrix_alloc(nx, nx);
	kf->R = matrix_alloc(nz, nz);
	kf->K = matrix_alloc(nx, nz);

	/* used as a constant */
	I4 = matrix_alloc(nx, nx);

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

	struct matrix *H;
	if(kf->is_6dof == true) 
		H = matrix_alloc(3, 4);
	else
		H = matrix_alloc(6, 4);
	if(H == NULL)
		return NULL;

	struct quaternion xq;
	xq.w = MATRIX_AT(kf->x, 0, 0);
	xq.x = MATRIX_AT(kf->x, 1, 0);
	xq.y = MATRIX_AT(kf->x, 2, 0);
	xq.z = MATRIX_AT(kf->x, 3, 0);

	struct quaternion aq;
	aq.w = 0.0;
	aq.x = 2.0f * kf->ax_ref; 
	aq.y = 2.0f * kf->ay_ref; 
	aq.z = 2.0f * kf->az_ref;	

	struct quaternion mq;
	mq.w = 0.0;
	mq.x = 2.0f * kf->mx_ref; 
	mq.y = 2.0f * kf->my_ref; 
	mq.z = 2.0f * kf->mz_ref;	

	struct quaternion ap = quat_mul(xq, aq);
	struct quaternion mp = quat_mul(xq, mq);
	
	MATRIX_SET(H, 0, 0, ap.x);
	MATRIX_SET(H, 0, 1, -ap.w);
	MATRIX_SET(H, 0, 2, ap.z);
	MATRIX_SET(H, 0, 3, -ap.y);

	MATRIX_SET(H, 1, 0, ap.y);
	MATRIX_SET(H, 1, 1, -ap.z);
	MATRIX_SET(H, 1, 2, -ap.w);
	MATRIX_SET(H, 1, 3, ap.x);

	MATRIX_SET(H, 2, 0, ap.z);
	MATRIX_SET(H, 2, 1, ap.y);
	MATRIX_SET(H, 2, 2, -ap.x);
	MATRIX_SET(H, 2, 3, -ap.w);

	if(kf->is_6dof == false) {
		MATRIX_SET(H, 3, 0, mp.x);
		MATRIX_SET(H, 3, 1, -mp.w);
		MATRIX_SET(H, 3, 2, mp.z);
		MATRIX_SET(H, 3, 3, -mp.y);

		MATRIX_SET(H, 4, 0, mp.y);
		MATRIX_SET(H, 4, 1, -mp.z);
		MATRIX_SET(H, 4, 2, -mp.w);
		MATRIX_SET(H, 4, 3, mp.x);

		MATRIX_SET(H, 5, 0, mp.z);
		MATRIX_SET(H, 5, 1, mp.y);
		MATRIX_SET(H, 5, 2, -mp.x);
		MATRIX_SET(H, 5, 3, -mp.w);
	}

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

	float s = kf->var_w * (kf->dt/2.0f) * (kf->dt/2.0f);
	matrix_smul(Q, s); 

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

	return W;
}
