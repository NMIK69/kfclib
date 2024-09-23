#include <stdlib.h>
#include <math.h>

#include "quaternion.h"

struct quaternion quat_mul(struct quaternion q1, struct quaternion q2)
{
	struct quaternion res = {0};

	res.w = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
	res.x = q1.x*q2.w + q1.w*q2.x - q1.z*q2.y + q1.y*q2.z;
	res.y = q1.y*q2.w + q1.z*q2.x + q1.w*q2.y - q1.x*q2.z;
	res.z = q1.z*q2.w - q1.y*q2.x + q1.x*q2.y + q1.w*q2.z;

	return res;
}

struct quaternion quat_norm(struct quaternion q)
{
	float w2 = q.w * q.w;
	float x2 = q.x * q.x;
	float y2 = q.y * q.y;
	float z2 = q.z * q.z;

	float n = sqrt((w2 + x2 + y2 + z2));

	struct quaternion res = {0};
	res.w = q.w / n;
	res.x = q.x / n;
	res.y = q.y / n;
	res.z = q.z / n;

	return res;
}
