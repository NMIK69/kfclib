## kfclib

This is a C library that provides an extended Kalman filter implementation for
tracking orientation, expressed as a quaternion. It is designed to be used with
a 6DoF or 9DoF IMU. This means that the filter expects angular velocity, acceleration
and magnetic field measurements as inputs.

In addition to the C implementation, you will also find a Python version that
was used for initial prototyping and design of the 6dof version. However, the Python code will not receive any updates in the future. The C implementation can be found in
[kfclib](kfclib), and the Python prototype can be found in [kfplib](kfplib).

> [!WARNING] 
> WHILE I HAVE MADE EFFORTS TO ENSURE THAT THE FILTER WORKS AS INTENDED, I CANNOT
> GUARANTEE THE CORRECTNESS OF THE FILTER DESIGN.

## Main Resources and References used
- [https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python](https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python)
- [https://ahrs.readthedocs.io/en/latest/filters/ekf.html#](https://ahrs.readthedocs.io/en/latest/filters/ekf.html#)
- [Basile Graf, arXiv: Dynamical Systems, 2008](https://arxiv.org/pdf/0811.2889)
- [Lightweight Extended Kalman Filter for MARG Sensors Attitude
Estimation](https://ieeexplore.ieee.org/document/9402867)
- [An Extended Kalman Filter for Quaternion-Based Orientation Estimation Using MARG Sensors](https://ieeexplore.ieee.org/document/976367)
- [Quaternion-Based Iterative Extended Kalman Filter for Sensor Fusion of Vision Sensor and IMU in 6-DOF Displacement Monitoring](https://ieeexplore.ieee.org/document/9925908)


## Usage 
### Initializing the filter
The filter should be initialized using the ```kf_init``` function:

```c
struct kalman_filter *kf_init(float dt, float var_a, float var_w, float var_m, float var_P, int kf_type);
```

This function takes the following parameters:
- ```dt```: Time difference between measurements (in seconds).
- ```var_a```: Variance of the accelerometer measurements.
- ```var_w```: Variance of the angular velocity measurements.
- ```var_m```: Variance of the magnetometer measurements.
- ```var_P```: Variance of the initial state.
- ```kf_type```: Select 6DoF or 9DoF variant. Can only by ```KF_6DOF``` or ```KF_9DOF```.


#### Setting the initial orientation and reference vectors. 
The filter is designed in a way such that there is no need to know the initial
orientation in relation to a fixed frame of reference (e.g. NED, ENU, etc.).
This leaves the user with two approaches for initializing the state and
expressing the orientation. 

##### Approach 1: Unknown Orientation relative to a fixed reference frame.
In this approach, you might average the first few accelerometer and magnetometer
readings (without inducing any translation or rotation) and use those as the
reference vectors. The initial orientation, i.e., the initial state vector,
should be the identity quaternion.

##### Approach 2: Known Orientation relative to a fixed reference frame.
If you want to express the orientation relative to a fixed frame of reference
(e.g., NED, ENU, etc.), you should use the gravitational and magnetic vector of
that frame as the reference vectors. The initial state vector should be set
to the orientation relative to the fixed frame. This might yield better results
depending on the state of the IMU and how accurately the initial orientation can
be determined

The quaternion that expresses the orientation can be set with:

```c
void kf_set_q(struct kalman_filter *kf, float qw, float qx, float qy, float qz);
```

The reference acceleration and magnetic vector can be set with:

```c
void kf_set_aref(struct kalman_filter *kf, float ax, float ay, float az);
void kf_set_mref(struct kalman_filter *kf, float mx, float my, float mz);
```

### Running the filter
Once the filter is initialized, you can use the ```kf_filt_6dof``` or
```kf_filt_9dof``` function to perform a single filter step (predict + update).

```c
int kf_filt_6dof(struct kalman_filter *kf, 
	   float wx, float wy, float wz,
	   float ax, float ay, float az);

int kf_filt_9dof(struct kalman_filter *kf, 
	   float wx, float wy, float wz,
	   float ax, float ay, float az,
	   float mx, float my, float mz);

```

The function takes a pointer to the ```struct kalman_filter``` that was obtained
from the ```kf_init``` funtion along with the angular velocity measurements
(```wx, wy, wz```) in $\frac{\text{rad}}{{s}}$, acceleration measurements
(```ax, ay, az```) and magnetometer measurements (```mx, my, mz```).

After each filter step, the orientation quaternion is updated and can be
accessed via the ```kf->q``` member, which is of type ```struct quaternion```:

```c
struct quaternion
{
	float w, x, y, z;
};
```

### Releasing memory
To release the filter's allocated memory, use ```kf_free```.

```c
void kf_free(struct kalman_filter *kf);
```

## Example
An example can be found in the [example](example) directory. The code initializes
and runs the filter on some real 6DoF and 9DoF IMU data. The results are then
visualized using `matplotlib`. For the plotting to work, you need `matplotlib`,
`numpy` and possibly `PyQt5`. All requirements can be installed using the
`requirements.txt` file from the example directory. You can run ```make
example_6dof``` or ```make example_9dof``` to run one of the two examples.


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

This leads to the Jacobian:

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
[here](https://ahrs.readthedocs.io/en/latest/filters/ekf.html#), I had better results. They designed $\mathbf{Q}$ the following way:

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
The accelerometer measurements are used for the measurement vector.

$$
\mathbf{z} = 
\begin{bmatrix}
a_x\\
a_y\\
a_z
\end{bmatrix}
$$

### Observation matrix
The state vector is converted into measurement space by rotating the reference
acceleration vector using the quaternion of the current state vector. If the
prediction was accurate, the rotated $\mathbf{a_{\text{ref}}}$ should
closely match the measured acceleration $\mathbf{z}$:

$$
\mathbf{H}\mathbf{x} = \mathbf{x} \otimes 
\begin{bmatrix}
0\\
\mathbf{a_{\text{ref}}}
\end{bmatrix}
\otimes \mathbf{x}^*
$$

This results in the Jacobian:

$$
\mathbf{H} = 2 *
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
\\mathbf{a_{\text{ref}}}
\end{bmatrix}
$$


### Measurement Covariance Matrix
The measurement covariance matrix is an identity matrix with the accelerometer
variance $\sigma^2_a$ on the diagonal.

$$
\mathbf{R} = \sigma^2_a \mathbf{I_3},
$$


## Kalman Filter Design (9dof-IMU)
The only thing that changes in the 9DoF variant is the design of the $H$
and $R$ matrix.

### Observation matrix
$$
\mathbf{H} = 2 *
\begin{bmatrix}
p_x & -p_w &  p_z & -ap_y\\
p_y & -p_z & -p_w &  ap_x\\
p_z &  p_y & -p_x & -ap_w\\
r_x & -r_w &  r_z & -r_y\\
r_y & -r_z & -r_w &  r_x\\
r_z &  r_y & -r_x & -r_w
\end{bmatrix}
$$

With $\mathbf{p} = {[p_w, p_x, p_y, p_z]}^T$ and
$\mathbf{r} = {[r_w, r_x, r_y, r_z]}^T$ and

$$
\mathbf{p} = \mathbf{x} \otimes 
\begin{bmatrix}
0\\
\mathbf{a_{\text{ref}}}
\end{bmatrix}
$$

$$
\mathbf{r} = \mathbf{x} \otimes 
\begin{bmatrix}
0\\
\mathbf{m_{\text{ref}}}
\end{bmatrix}
$$



### Measurement Covariance Matrix
The measurement covariance matrix is an identity matrix with the accelerometer
variance $\sigma^2_a$ on the upper left diagonal and the magnetometer variance
$\sigma^2_m$ on the lower right diagonal:

$$
\mathbf{R} = 
\begin{bmatrix}
\sigma^2_a \mathbf{I_3} & 0\\
0 & \sigma^2_m \mathbf{I_3}
\end{bmatrix}
$$
