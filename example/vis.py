import numpy as np
import sys
import matplotlib.pyplot as plt
import math

np.set_printoptions(suppress=True, linewidth=180)


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
	
	quats = []
	k = 0
	
	while(k < len(lines)):
		measurements = lines[k].split(",")
		q = []
		for i in range(0, 4):
		   q.append(float(measurements[i]))
	
		quats.append(q)
		k += 1
	
	return quats 

filename = sys.argv[1]
quats = extract_data_from_file(filename)
time = np.linspace(0, (len(quats) * 1/200.0), num=len(quats))

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














