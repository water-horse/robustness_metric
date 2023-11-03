f0 = open("traj0.txt", 'r')

n0 = 666

f1 = open("traj_delta2.txt", 'r')

n1 = 330

fout = open("traj2.txt", 'w')

for i in range(n0 - n1):
    fout.write(f0.readline())

for i in range(n1):
    line0 = f0.readline()
    line1 = f1.readline()
    data0 = line0.split(' ')
    data1 = line1.split(' ')
    x = float(data0[1]) + float(data1[1])
    y = float(data0[2]) + float(data1[2])
    z = float(data0[3]) + float(data1[3])
    fout.write(data0[0] + " " + str(x) + " " + str(y) + " " + str(z) + " " + "1.0 0.0 0.0 0.0\n")

f0.close()
f1.close()
fout.close()