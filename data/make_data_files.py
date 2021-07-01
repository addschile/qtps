import os
from sys import argv

tobs = int(argv[1])

for i in range(6):
  if i != 4:
    for j in range(10):
      if j==0:
        os.system('tail -10000 tobs%d/window%d/hb_list_tobs_%d_task_%d.dat > tobs%d/window%d/window%d.dat'%(tobs,i,tobs,j,tobs,i,i))
      else:
        os.system('tail -10000 tobs%d/window%d/hb_list_tobs_%d_task_%d.dat >> tobs%d/window%d/window%d.dat'%(tobs,i,tobs,j,tobs,i,i))
