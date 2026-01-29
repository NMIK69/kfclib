#ifndef KALMAN_H
#define KALMAN_H

#include "matrix.h"
#include "quaternion.h"
#include <stdbool.h>

#define KF_6DOF 0x01
#define KF_9DOF 0x02

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

	float mx_ref;
	float my_ref;
	float mz_ref;

	struct quaternion q;

	bool is_6dof;
};

struct kalman_filter *kf_init(float dt, float var_a, float var_w,
				float var_m, float var_P, int type);
void kf_free(struct kalman_filter *kf);

int kf_filt_6dof(struct kalman_filter *kf, 
	   float wx, float wy, float wz,
	   float ax, float ay, float az);

int kf_filt_9dof(struct kalman_filter *kf, 
	   float wx, float wy, float wz,
	   float ax, float ay, float az,
	   float mx, float my, float mz);

void kf_set_q(struct kalman_filter *kf, float qw, float qx, float qy, float qz);
void kf_set_aref(struct kalman_filter *kf, float ax, float ay, float az);
void kf_set_mref(struct kalman_filter *kf, float mx, float my, float mz);

#endif
