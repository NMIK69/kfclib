#ifndef KALMAN_H
#define KALMAN_H

#include "matrix.h"
#include "quaternion.h"

struct kalman_filter
{
	float dt;
	float var_w;
	float var_a;
	float var_P;

	struct matrix *x;	

	struct matrix *P;	
	struct matrix *R;	
	struct matrix *K;

	float ax_ref;
	float ay_ref;
	float az_ref;

	struct quaternion q;
};

struct kalman_filter *kf_init(float dt, float var_a, float var_w, float var_P);
void kf_free(struct kalman_filter *kf);

int kf_filt(struct kalman_filter *kf, 
	   float wx, float wy, float wz,
	   float ax, float ay, float az);

void kf_set_q(struct kalman_filter *kf, float qw, float qx, float qy, float qz);
void kf_set_aref(struct kalman_filter *kf, float ax, float ay, float az);

#endif
