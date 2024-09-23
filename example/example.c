#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <math.h>

#include "../kfclib/matrix.h"
#include "../kfclib/kalman.h"

#define ARR_SIZE(arr) (sizeof(arr) / sizeof(*arr))
#define MM_PI 3.14159265358979323846f


struct measurement
{
	float wx, wy, wz;
	float ax, ay, az;
};

struct vec3f
{
	float x, y, z;
};

static ssize_t get_nlines(const char *fname);
static struct measurement *read_measurements(const char *fname, size_t *len);
static struct vec3f get_aref(struct measurement *mea, size_t len);


int main(int argc, char **argv)
{
	assert(argc == 2);

	size_t len = 0;
	struct measurement *mea = read_measurements(argv[1], &len);
	assert(mea != NULL);


	struct vec3f aref = get_aref(mea, len);
	float dt = 1.0f/200.0f;
	float var_a = 0.8*0.8;
	float var_w = 0.1*0.1;
	float var_P = 0.0001;
	struct kalman_filter *kf = kf_init(dt, var_a, var_w, var_P);
	kf_set_aref(kf, aref.x, aref.y, aref.z);


	FILE *fout = fopen("out.txt", "w"); 
	assert(fout != NULL);


	for(size_t i = 0; i < len; i++) {
		int err = kf_filt(kf, mea[i].wx, mea[i].wy, mea[i].wz,
			          mea[i].ax, mea[i].ay, mea[i].az);
		assert(err == 0);

		fprintf(fout, "%f,%f,%f,%f\n", kf->q.w, kf->q.x, kf->q.y, kf->q.z);
		assert(ferror(fout) == 0);
	}

	printf("P:\n");
	matrix_print(kf->P);
	printf("\nK:\n");
	matrix_print(kf->K);

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

static struct measurement *read_measurements(const char *fname, size_t *len)
{
	char entry[1024];
	FILE *f = fopen(fname, "r");
	assert(f != NULL);

	ssize_t nlines = get_nlines(fname);
	assert(nlines != -1);

	struct measurement *mea = malloc(sizeof(*mea) * nlines);
	if(mea == NULL)
		return NULL;
	
	size_t i = 0;
	int ts;
	while(feof(f) == 0 && fgets(entry, ARR_SIZE(entry), f) != NULL) {
		
		int ret = sscanf(entry, "%d,%f,%f,%f,%f,%f,%f", 
				&ts,
				&mea[i].ax, &mea[i].ay, &mea[i].az,
				&mea[i].wx, &mea[i].wy, &mea[i].wz);
		mea[i].wx = mea[i].wx * MM_PI / 180.0f;
		mea[i].wy = mea[i].wy * MM_PI / 180.0f;
		mea[i].wz = mea[i].wz * MM_PI / 180.0f;

		(*len) += 1;
		assert(ret == 7);

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
