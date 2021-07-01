import numpy as np
from sys import argv

tobs = int(argv[1])
p0   = np.zeros(10)
p2   = np.zeros(10)
p1   = np.zeros(10)
Zab  = np.zeros(10)
rate = np.zeros(10)

for i in range(10):
  da = np.loadtxt('tobs%d/reweighted_hist_%d.dat'%(tobs,i))
  p0[i] = np.exp(-da[-2,1])
  p2[i] = np.exp(-da[-1,1])
  p1[i] = np.exp(-da[-3,1])
Zab = p1/(p0+p2)

f = open('tobs%d/path_partition_function_%d.dat'%(tobs,tobs),'w')
for i in range(10):
  f.write('%d %.16f\n'%(i,Zab[i]))

Zab_avg = np.sum(Zab[:])/10.
for i in range(10):
  Zab[i] -= Zab_avg
Zab *= Zab
std_err = np.sqrt(np.sum(Zab[:])/10.)
f.write('%.16f %.16f\n'%(Zab_avg,std_err))
f.close()
