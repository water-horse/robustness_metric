import numpy as np

f0 = open("traj0.txt", 'r')

n0 = 666

ts = []
xs = []
ys = []
zs = []

x = np.arange(n0)
xp = np.arange(((n0 - 1)//18)+1) * 18

for i in range(n0):
    line0 = f0.readline()
    data0 = line0.split(' ')
    ts.append(data0[0])
    if i % 18 == 0:
        data0 = line0.split(' ')
        # xs.append(float(data0[1]) + np.random.rand()*1 - 0.5)
        xs.append(float(data0[1]))
        # ys.append(float(data0[2]) + np.random.rand()*1 - 0.5)
        ys.append(float(data0[2]))
        # zs.append(float(data0[3]) + np.random.rand()*1 - 0.5)
        zs.append(float(data0[3]))

f0.close()

xs_interp = np.interp(x, xp, np.array(xs))
ys_interp = np.interp(x, xp, np.array(ys))
zs_interp = np.interp(x, xp, np.array(zs))

fout = open("traj4.txt", 'w')


for i in range(n0):
    fout.write(ts[i] + " " + str(xs_interp[i]) + " " + str(ys_interp[i]) + " " + str(zs_interp[i]) + " 1.0 0.0 0.0 0.0\n")

fout.close()
