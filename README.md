## kfclib

This is a library that provides an extended Kalman filter implementation for
tracking orientation, expressed as a quaternion. It is designed to be used with
a 6DOF IMU. This means that the filter expects angular velocity and acceleration
as input.

In addition to the C implementation, you will also find a Python version that
was used for initial prototyping and design. The C implementation is in
[kfclib](kfclib), and the Python prototype can be found in [kfplib](kfplib).

> [!WARNING] 
> WHILE I HAVE MADE EFFORTS TO ENSURE THAT THE FILTER WORKS AS INTENDED, I CANNOT
> GUARANTEE THE CORRECTNESS OF THE FILTER DESIGN.

**The 9dof version is coming soon.**

## Main Resources and References used
- [https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python](https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python)
- [https://ahrs.readthedocs.io/en/latest/filters/ekf.html#](https://ahrs.readthedocs.io/en/latest/filters/ekf.html#)
- [Basile Graf, arXiv: Dynamical Systems, 2008](https://arxiv.org/pdf/0811.2889)


## Usage 
### Initializing the filter
The filter should be initialized using the ```kf_init``` function:

```c
struct kalman_filter *kf_init(float dt, float var_a, float var_w, float var_P);
```

This function takes the following parameters:
- ```dt```: Time difference between measurements (in seconds).
- ```var_a```: Variance of the accelerometer measurements.
- ```var_w```: Variance of the angular velocity measurements.
- ```var_P```: Variance of the initial state.


#### Setting the initial orientation and acceleration reference vector. 
The filter is designed in a way such that there is no need to know the initial
orientation in relation to a fixed frame of reference (e.g. NED, ENU, etc.).
This leaves the user with two approaches for initializing the state and
expressing the orientation. 

##### Approach 1: Unknown Orientation relative to a fixed reference frame.
In this approach, you might average the first few accelerometer readings
(without inducing any translation or rotation) and use those as the reference
acceleration. The initial orientation, i.e., the initial state vector, should be
the identity quaternion.

##### Approach 2: Known Orientation relative to a fixed reference frame.
If you want to express the orientation relative to a fixed frame of reference
(e.g., NED, ENU, etc.), you should use the gravitational vector of that frame as
the reference acceleration. The initial state vector should be set to the
orientation relative to the fixed frame. This might yield better results depending
on the state of the IMU and how accurately the initial orientation can be
determined

The quaternion that expresses the orientation can be set with:

```c
void kf_set_q(struct kalman_filter *kf, float qw, float qx, float qy, float qz);
```

The reference acceleration can be set with:

```c
void kf_set_aref(struct kalman_filter *kf, float ax, float ay, float az);
```

### Running the filter
Once the filter is initialized, you can use the ```kf_filt``` function to
perform a single filter step (predict + update).

```c
int kf_filt(struct kalman_filter *kf, 
	   float wx, float wy, float wz,
	   float ax, float ay, float az);

```

The function takes a pointer to the ```struct kalman_filter``` that was
obtained from ```kf_init``` along with the angular velocity measurements (```wx,
wy, wz``` in $\frac{\text{rad}}{{s}}$ and the acceleration measurements (```ax,
ay, az```).

After each filter step, the orientation quaternion is updated and can be
accessed via the ```kf->q``` member, which is of type ```struct quaternion```:

```c
struct quaternion
{
	float w, x, y, z;
};
```

### Releasing the filter
To release the filter's allocated memory, use ```kf_free```.

```c
void kf_free(struct kalman_filter *kf);
```

## Example
**Under construction**.


## Kalman Filter Design (6dof-IMU)

### State vector
The state vector is a quaternion that tracks the orientation.

$$
\mathbf{x} = 
\begin{bmatrix}
q_w\\
q_x\\
q_y\\
q_z
\end{bmatrix}
$$

### State Covariance Matrix
The state covariance matrix is initialized as an identity matrix since there
should be no covariance between the quaternion elements.

$$
\mathbf{P} = \sigma^2_p \mathbf{I_4}
$$

### State Transition Matrix
The gyroscope measurements are used as control input to predict the next state (i.e., compute the prior). Therefore the quaternion integration
equation is used to create the state transition matrix:

$$
q_n = q_{n-1} + \Delta t \frac{dq_{n}}{dt}, \quad  
\frac{dq_n}{dt} = \frac{1}{2}q_{n-1} \otimes \Omega_n,
$$

Where:

$$
\Omega_n = {[0, w_x, w_y, w_z]}^T
$$

Here $w_x$, $w_y$ and $w_z$ are the
angular velocity measurements from the gyroscope in $\frac{\text{rad}}{s}$. The
$\otimes$ operator denotes quaternion multiplication. For this kalman filter,
$q_n$ is the current prior, and $q_{n-1}$ is the previous posterior. The state
transition matrix is derived from the quaternion integration equation such that $\mathbf{F}$ is linear and such that:

$$
\mathbf{x_n} = \mathbf{F}\mathbf{x_{n-1}} = \mathbf{x_{n-1}} + \Delta t (\frac{1}{2}\mathbf{x_{n-1}} \otimes \Omega_n)
$$

This leads to:

$$
\mathbf{F} = 
\begin{bmatrix}
1 & \frac{-\Delta t}{2}\omega_x & \frac{-\Delta t}{2}\omega_y & \frac{-\Delta t}{2}\omega_z \\
\frac{\Delta t}{2}\omega_x & 1 & \frac{\Delta t}{2}\omega_z & \frac{-\Delta t}{2}\omega_y \\
\frac{\Delta t}{2}\omega_y & \frac{-\Delta t}{2}\omega_z & 1 & \frac{\Delta t}{2}\omega_x \\
\frac{\Delta t}{2}\omega_z & \frac{\Delta t}{2}\omega_y & \frac{-\Delta t}{2}\omega_x & 1
\end{bmatrix}
$$


### Process Noise Covariance Matrix
In the beginning I just used an identity matrix with the variance of the
gyroscope measurements $\sigma^2_\omega$ on the diagonal. But after using a
design initially found
[here](https://ahrs.readthedocs.io/en/latest/filters/ekf.html#) which is linked
to [this
publication](https://www.researchgate.net/publication/268554081_Gyro_Modeling_and_Estimation_of_Its_Random_Noise_Sources)
I had better results. They designed $\mathbf{Q}$ the following way:

$$
\mathbf{Q_n} = \sigma^2_\omega\mathbf{W_n}\mathbf{W_n^T}
$$


$$
\mathbf{W_n} = \frac{\Delta t}{2}
\begin{bmatrix}
-q_x & -q_y & -q_z\\
q_w & -q_z & q_y\\
q_z & q_w & -q_x\\
-q_y & q_x & q_w\\
\end{bmatrix}
$$

Where $q_w, q_x, q_y, q_z$ are taken from the current state vector.

### Measurement vector
The normalized accelerometer measurements are used for the measurement vector.

$$
\mathbf{z} = 
\begin{bmatrix}
\hat{a_x}\\
\hat{a_y}\\
\hat{a_z}
\end{bmatrix}
$$

### Observation matrix
The state vector is converted into measurement space by rotating the reference acceleration vector using the quaternion of the current state vector. If the prediction was accurate, the rotated $\hat{\mathbf{a_{\text{ref}}}}$ should closely match the measured acceleration $\mathbf{z}$. The observation matrix $\mathbf{H}$ is computed such that $\mathbf{H}$ is linear and such that:

$$
\mathbf{H}\mathbf{x} = \mathbf{x} \otimes 
\begin{bmatrix}
0\\
\hat{\mathbf{a_{\text{ref}}}}
\end{bmatrix}
\otimes \mathbf{x}^*
$$

This results in:

$$
\mathbf{H} = 
\begin{bmatrix}
p_x & -p_w & p_z & -p_y\\
p_y & -p_z & -p_w & p_x\\
p_z & p_y & -p_x & -p_w
\end{bmatrix}
$$

With $\mathbf{p} = {[p_w, p_x, p_y, p_z]}^T$ and:

$$
\mathbf{p} = \mathbf{x} \otimes 
\begin{bmatrix}
0\\
\hat{\mathbf{a_{\text{ref}}}}
\end{bmatrix}
$$

### Measurement Covariance Matrix
The measurement covariance matrix is an identity matrix with the accelerometer
variance $\sigma^2_a$ on the diagonal.

$$
\mathbf{R} = \sigma^2_a \mathbf{I_3},
$$


## Kalman Filter Design (9dof-IMU)
**Coming soon**.
