#ifndef QUATERNION_H
#define QUATERNION_H

struct quaternion
{
	float w;
	float x;
	float y;
	float z;
};


struct quaternion quat_mul(struct quaternion q1, struct quaternion q2);
struct quaternion quat_norm(struct quaternion q);

#endif
