## kfclib
> [!WARNING] 
> WHILE I HAVE MADE EFFORTS TO ENSURE THAT THE FILTER WORKS AS INTENDED, I CANNOT
> GUARANTEE THE CORRECTNESS OF THE FILTER DESIGN.

This is a library that provides an extended kalman filter implementation that
tracks orientation expressed as a quaternion. It is intended to be used with a
6dof-IMU. This means, that the filter expects angular velocity and acceleration
as its input. In order to run the filter on a MCU, I tried to make the
implementation as efficient as possible. Thus the use of a linear arena allocator.
Besides the C implementation you can also find a Python implementation, that I
used for the initial prototyping and designing. The C version is in
[kfclib](kfclib) and the Python prototype in [kfplib](kfplib).

**The 9dof version is comming soon.**

## Main Resources and References used
- [https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python](https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python)
- [https://ahrs.readthedocs.io/en/latest/filters/ekf.html#](https://ahrs.readthedocs.io/en/latest/filters/ekf.html#)
- [Basile Graf, arXiv: Dynamical Systems, 2008](https://arxiv.org/pdf/0811.2889)


## Usage 
### Initilizing the filter
The filter should be initilized using the ```kf_init``` function:

```c
struct kalman_filter *kf_init(float dt, float var_a, float var_w, float var_P);
```

The function takes in the time difference between measurements in seconds
(```dt```), the variance of the accelerometer measurements (```var_a```), the
variance of the angular velocity measurements (```var_w```) and the variance of
the initial state (```var_P```).

#### Setting the initial orientation and acceleration reference. 
The filter is designed in a way such that there is no need to know the initial
orientation in relation to a fixed frame of reference (E.g. NED, ENU, etc.).
This leaves the user with two approaches of initilizing the state and
expressing the orientation. 

The quaternion that expresses the orientation can be
set with:
```c
void kf_set_q(struct kalman_filter *kf, float qw, float qx, float qy, float qz);
```

The reference acceleration can be set with:
```c
void kf_set_aref(struct kalman_filter *kf, float ax, float ay, float az);
```

##### Approach 1: Unkonwn Orientation in relation to a fixed reference frame.
When using this apporach, it is probably best to average the first few
accelerometer readings (without inducing any translation or rotation) and using
those as the reference acceleration. The initial orientation aka the initial state
vector should be the identity quaternion.

##### Approach 2: Known Orientation in relation to a fixed reference frame.
If you want express the orientation of in relation to some fixed frame
of reference (E.g. NED, ENU, etc.) you would use it's gravitational vector as the
reference acceleration and the initial orientation (in relation the
the fixed frame) as the initial state vector. This might lead to better results
depending on the state of the IMU and how accurate the initial orientation can
be determined.

### Running the filter
Once you initilized the filter, you can use the ```kf_filt``` function to do one
filter step (predict + update).

```c
int kf_filt(struct kalman_filter *kf, 
	   float wx, float wy, float wz,
	   float ax, float ay, float az);

```

The function takes a poninter to the ```struct kalman_filter``` that was
obtained from ```kf_init``` as well as the angular velocity measurements (```wx,
wy, wz``` in $\frac{\text{rad}}{{s}}$ and the acceleration measuremens (```ax,
ay, az```).

After each filter step, the orientation quaternion is updated and can be
accessed throught ```kf->q```.

### Releasing the filter
The filter's memory can be released with ```kf_free```.

```c
void kf_free(struct kalman_filter *kf);
```

## Example
An example can be found in the [example](example) folder. The example uses
measurements obtained from an mpu6050 IMU. The IMU was turned roughly 90 deg
back and forth around each axis by hand. The IMU was calibrated as well.

### Building and Running the example
- Navigate into the example folder.
- run ```make```
- run ```./example```
This will run the kalman filter and ouput the results in a file called
```out.txt```. In order to visualize the results you can run ```python3 vis.py
out.txt```. This is a Python script that uses matplotlib to visualize the
estimted orientation over time.


## Kalman Filter Design
### 6dof-IMU
#### State vector
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

## State Covariance Matrix
The state covariance matrix is initialized as an identity matrix since there
should be no covariance between the quaternion elements.

$$
\mathbf{P} = \sigma^2_p \mathbf{I_4}
$$

$$
\sigma^2_p = 0.001
$$

## State Transition Matrix
The measurements of the gyroscope are used as control input in order to predict
the next state aka compute the prior. In particular, the quaternion integration
equation is used to create the state transition matrix.

$$
q_n = q_{n-1} + \Delta t \frac{dq_{n}}{dt}, \quad  
\frac{dq_n}{dt} = \frac{1}{2}q_{n-1} \otimes \Omega_n,
$$

With $\Omega_n = {[0, w_x, w_y, w_z]}^T$ where $w_x$, $w_y$ and $w_z$ are the
angular velocity measurements or the gyroscope in $\frac{\text{rad}}{s}$. The
$\otimes$ operator denotes quaternion multiplication. For this kalman filter,
$q_n$ is the current prior and $q_{n-1}$ the previous posterior. The state
transition matrix is derived from the quaternion integration equation such that:

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


## State Transition Covariance Matrix
In the beginning i just used an identity matrix with the variance of the
gyroscope measurements $\sigma^2_\omega$ on the diagonal. But after using a
design initially found
[here](https://ahrs.readthedocs.io/en/latest/filters/ekf.html#) which is linked
to [this
publication](https://www.researchgate.net/publication/268554081_Gyro_Modeling_and_Estimation_of_Its_Random_Noise_Sources)
i had better results. They designed $\mathbf{Q}$ the following way:

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

## Measurement vector
The normed measurements of the accelerometer are used to compute the posterior.

$$
\mathbf{z} = 
\begin{bmatrix}
a_x\\
a_y\\
a_z
\end{bmatrix}
$$

## Measurement model matrix
Here is where the acceleration reference vector $\mathbf{a_{\text{ref}}}$ comes
into play. We somehow need to convert the state vector into measurement space.
This is done by rotating the acceleration reference vector by the quaternion of
the current state vector. If our prediction was accurate, then the rotated
$\mathbf{a_{\text{ref}}}$ should be realy close to the measured acceleration
$\mathbf{z}$. The measurement model matrix $\mathbf{H}$ is computed such that:

$$
\mathbf{H}\mathbf{x} = \mathbf{x} \otimes 
\begin{bmatrix}
0\\
\mathbf{a_{\text{ref}}}
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
\mathbf{a_{\text{ref}}}
\end{bmatrix}
$$

## Measurement Covariance Matrix
The measurement covariance matrix is an identity matrix with the accelerometer
variance $\sigma^2_a$ on the diagonal.

$$
\mathbf{R} = \sigma^2_a \mathbf{I_3},
$$


# Expanding to 9dof
**Coming soon**.
