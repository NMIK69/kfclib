import numpy as np
import sys
import matplotlib.pyplot as plt
import math

np.set_printoptions(suppress=True, linewidth=180)


var_a = 0.8**2
var_g = 0.1**2
var_p = 0.001

#var_a = 0.5**2
#var_g = 0.3**2
#var_p = 0.001

x = np.array([1.0, 0.0, 0.0, 0.0])
R = np.eye(3) * var_a
P = np.eye(4) * var_p
K = []
a_ref = []
dt = 1/200.0
#dt = 1/166.0


def quat_mul(q, p):
    p0, p1, p2, p3 = p
    q0, q1, q2, q3 = q

    w = p0*q0 - q1*p1 - q2*p2 - q3*p3
    x = q1*p0 + q0*p1 + q2*p3 - q3*p2
    y = q2*p0 + q0*p2 + q3*p1 - q1*p3
    z = q3*p0 + q0*p3 + q1*p2 - q2*p1

    return np.array([w,x,y,z])

def norm_quat(q):
    a = q[0]
    b = q[1]
    c = q[2]
    d = q[3]

    n = np.sqrt(a**2 + b**2 + c** 2 + d**2)

    assert(n != 0)	
    return [a / n, b / n, c / n, d / n]

def quat_to_euler(q):
    w, x, y, z = q

    t0 = 2.0 * (w * x + y * z)
    t1 = 1.0 - 2.0 * (x * x + y * y)
    xrot = math.atan2(t0, t1)
    
    t2 = 2.0 * (w * y - z * x)
    t2 = 1.0 if t2 > 1.0 else t2
    t2 = -1.0 if t2 < -1.0 else t2
    yrot = math.asin(t2)
    
    t3 = 2.0 * (w * z + x * y)
    t4 = 1.0 - 2.0 * (y * y + z * z)
    zrot = math.atan2(t3, t4)
    
    return [xrot, yrot, zrot]


def extract_data_from_file(filename):
	fp = open(filename, "r")
	lines = fp.readlines()
	
	accel = []
	gyro = []
	k = 0
	
	while(k < len(lines)):
		measurements = lines[k].split(",")
		accel_xyz = []
		gyro_xyz = []
	
		for i in range(1, 4):
			accel_xyz.append(float(measurements[i]))
	
		for i in range(4, 7):
			gyro_xyz.append(float(measurements[i]) * math.pi / 180)

		accel.append(accel_xyz)
		gyro.append(gyro_xyz)
		k += 1
	
	return accel, gyro


def kalman_predict(wx, wy, wz):
	global x
	global P

	s = dt/2.0

	qw = x[0]
	qx = x[1]
	qy = x[2]
	qz = x[3]


	W = (dt/2.0) * np.array([[-qx, -qy, -qz],
	      		 [qw, -qz, qy],
	      		 [qz, qw, -qx],
	      		 [-qy, qx, qw]])

	Q = var_g * (W@(W.T))

	F = np.array([[1.0, -s*wx, -s*wy, -s*wz],
		      [s*wx, 1.0, s*wz, -s*wy],
		      [s*wy, -s*wz, 1.0, s*wx],
		      [s*wz, s*wy, -s*wx, 1.0]])

	x = F@x

	P = F@P@(F.T) + Q

	x = norm_quat(x)

def kalman_update(ax, ay, az):
	global x
	global P
	global K

	gx = a_ref[0]
	gy = a_ref[1]
	gz = a_ref[2]

	qw = x[0]
	qx = x[1]
	qy = x[2]
	qz = x[3]


	q = np.array([qw, qx, qy, qz])
	g = np.array([0.0, gx, gy, gz])
	p = quat_mul(q, g)
	pw = p[0]
	px = p[1]
	py = p[2]
	pz = p[3]
	
	H = np.array([[px, -pw, pz, -py],
		      [py, -pz, -pw, px],
		      [pz, py, -px, -pw]])


	I4 = np.eye(4)
	z = np.array([ax, ay, az])

	v = z - H@x
	S = H@P@(H.T) + R
	S_inv = np.linalg.inv(S)
	K = P@(H.T)@S_inv
	x = x + K@v
	P = (I4 - K@H)@P

	x = norm_quat(x)


filename = sys.argv[1]
accel, gyro = extract_data_from_file(filename)
time = np.linspace(0, (len(accel) * 1/200.0), num=len(gyro))

accel = np.array(accel)
gyro = np.array(gyro)
a_ref = sum(accel[0:20])/20.0
a_ref = np.array(a_ref)
a_ref = a_ref / np.linalg.norm(a_ref)
print(a_ref)

quats = []

for t in range(len(accel)):
	
	kalman_predict(gyro[t][0], gyro[t][1], gyro[t][2])

	a = np.array([accel[t][0], accel[t][1], accel[t][2]])
	a = a / np.linalg.norm(a)
	kalman_update(a[0], a[1], a[2])

	quats.append([x[0], x[1], x[2], x[3]])


print(P)
print("\n")
print(K)

xr = []
yr = []
zr = []
for q in quats:
	rots = quat_to_euler(q)
	xr.append(rots[0] * 180.0/np.pi)
	yr.append(rots[1] * 180.0/np.pi)
	zr.append(rots[2] * 180.0/np.pi)

plt.plot(time, xr, label="X Rot", color='orange')
plt.plot(time, yr, label="Y Rot", color='blue')
plt.plot(time, zr, label="Z Rot", color='green')
plt.xlabel("time [ $s$ ]")
plt.ylabel("Rotation [ $^{\circ}$ ]")
plt.title("Rotation angles")
plt.legend()
plt.grid()
plt.show()














