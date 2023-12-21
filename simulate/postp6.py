import numpy as np

f0 = open("traj0.txt", 'r')

n0 = 666

fout = open("traj6.txt", 'w')

for i in range(n0):
    line0 = f0.readline()
    data0 = line0.split(' ')
    x = float(data0[1]) + np.random.rand() * 4.0 - 2.0
    y = float(data0[2]) + np.random.rand() * 4.0 - 2.0
    z = float(data0[3]) + np.random.rand() * 4.0 - 2.0
    fout.write(data0[0] + " " + str(x) + " " + str(y) + " " + str(z) + " " + "1.0 0.0 0.0 0.0\n")

fout.close()