#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <math.h>
#include <time.h>

#include "../matrix.h"
#include "../kalman.h"

#define ARR_SIZE(arr) (sizeof(arr) / sizeof(*arr))
#define MM_PI 3.14159265358979323846f


struct measurement
{
	float wx, wy, wz;
	float ax, ay, az;
	float mx, my, mz;
};

struct vec3f
{
	float x, y, z;
};

static ssize_t get_nlines(const char *fname);
static struct measurement *read_measurements(const char *fname, size_t *len, int);
static struct vec3f get_aref(struct measurement *mea, size_t len);
static struct vec3f get_mref(struct measurement *mea, size_t len);


int main(int argc, char **argv)
{
	int err;
	assert(argc == 3);

	int kf_type = argv[2][0] == '6' ? KF_6DOF : KF_9DOF;

	size_t len = 0;
	struct measurement *mea = read_measurements(argv[1], &len, kf_type);
	assert(mea != NULL);


	/* initilization parameters */
	// {
	struct vec3f aref = get_aref(mea, len);
	struct vec3f mref = get_mref(mea, len);

	float dt = 1.0f/200.0f;
	float var_a = 0.8*0.8;
	float var_w = 0.3*0.3;
	float var_m = 300*300;
	float var_P = 0.001;
	// }

	/* actual filter initilization */
	// {
	struct kalman_filter *kf = kf_init(dt, var_a, var_w, var_m, var_P, kf_type);
	kf_set_aref(kf, aref.x, aref.y, aref.z);
	kf_set_q(kf, 1.0f, 0.0f, 0.0f, 0.0f);
	kf_set_mref(kf, mref.x, mref.y, mref.z);
	// }

	FILE *fout = fopen("out.txt", "w"); 
	assert(fout != NULL);

	/* run filter */
	// {
	for(size_t i = 0; i < len; i++) {
		if(kf_type == KF_9DOF) {
			err = kf_filt_9dof(kf, mea[i].wx, mea[i].wy, mea[i].wz,
				               mea[i].ax, mea[i].ay, mea[i].az,
					       mea[i].mx, mea[i].my, mea[i].mz);
		}
		else {
			err = kf_filt_6dof(kf, mea[i].wx, mea[i].wy, mea[i].wz,
			          	       mea[i].ax, mea[i].ay, mea[i].az);
		}


		assert(err == 0);

		fprintf(fout, "%f,%f,%f,%f\n", kf->q.w, kf->q.x, kf->q.y, kf->q.z);
		assert(ferror(fout) == 0);
	}
	// }

	kf_free(kf);
	free(mea);
	fclose(fout);

	return 0;
}

static ssize_t get_nlines(const char *fname)
{
	char entry[1024];
	FILE *f = fopen(fname, "r");

	if(f == NULL)
		return -1;

	ssize_t count = 0;

	while(feof(f) == 0 && fgets(entry, ARR_SIZE(entry), f) != NULL) {
		count += 1;
	}

	fclose(f);
	return count;
}

static struct measurement *read_measurements(const char *fname, size_t *len, int kf_type)
{
	int ret;
	char entry[1024];
	FILE *f = fopen(fname, "r");
	assert(f != NULL);

	ssize_t nlines = get_nlines(fname);
	assert(nlines != -1);

	struct measurement *mea = malloc(sizeof(*mea) * nlines);
	if(mea == NULL)
		return NULL;
	
	size_t i = 0;
	while(feof(f) == 0 && fgets(entry, ARR_SIZE(entry), f) != NULL) {
		
		if(kf_type == KF_9DOF) {
			ret = sscanf(entry, "%*f,%f,%f,%f,%f,%f,%f,%f,%f,%f", 
					&mea[i].ax, &mea[i].ay, &mea[i].az,
					&mea[i].wx, &mea[i].wy, &mea[i].wz,
					&mea[i].mx, &mea[i].my, &mea[i].mz);
			assert(ret == 9);
		}
		else {
			ret = sscanf(entry, "%*f,%f,%f,%f,%f,%f,%f", 
					&mea[i].ax, &mea[i].ay, &mea[i].az,
					&mea[i].wx, &mea[i].wy, &mea[i].wz);
			assert(ret == 6);
		}
		mea[i].wx = mea[i].wx * MM_PI / 180.0f;
		mea[i].wy = mea[i].wy * MM_PI / 180.0f;
		mea[i].wz = mea[i].wz * MM_PI / 180.0f;

		(*len) += 1;

		i += 1;
	}

	fclose(f);

	return mea;
}

static struct vec3f get_aref(struct measurement *mea, size_t len)
{
	assert(len >= 20);

	struct vec3f a = {0};

	for(size_t i = 0; i < 20; i++) {
		a.x += mea[i].ax;	
		a.y += mea[i].ay;	
		a.z += mea[i].az;	
	}

	a.x = a.x / 20.0f;
	a.y = a.y / 20.0f;
	a.z = a.z / 20.0f;

	return a;
}

static struct vec3f get_mref(struct measurement *mea, size_t len)
{
	assert(len >= 20);

	struct vec3f m = {0};

	for(size_t i = 0; i < 20; i++) {
		m.x += mea[i].mx;	
		m.y += mea[i].my;	
		m.z += mea[i].mz;	
	}

	m.x = m.x / 20.0f;
	m.y = m.y / 20.0f;
	m.z = m.z / 20.0f;

	return m;
}
