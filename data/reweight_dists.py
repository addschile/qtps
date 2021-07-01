import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from scipy.optimize import minimize

def opt_func(cs,data_list,windows):
  err = 0.
  # window 5 and window 0
  err += (data_list[-1][-1] + cs[-1] - data_list[0][0] - cs[0])**2.
  # window 0 and window 1
  err += (data_list[0][-1] + cs[0] - data_list[1][0] - cs[1])**2.
  # window 1 and window 2
  err += (data_list[1][1] + cs[1] - data_list[2][1] - cs[2])**2.
  # window 2 and window 3
  err += (data_list[2][0] + cs[2] - data_list[3][-1] - cs[3])**2.
  # window 3 and window 4
  err += (data_list[3][0] + cs[3] - data_list[4][0] - cs[4])**2.
  return err

tobs = int(argv[1])
windows = [0,1,2,3,4,5]
window_states = [[2,4],[4,5,6,7,8],[3,5],[1,3],[1],[0,2]]

for i in range(10):
  data_list = []
  for window in windows:
    if window != 4:
      da = np.loadtxt('tobs%d/window%d/window%d.dat'%(tobs,window,window))
      #da = np.loadtxt('tobs%d/window%d/window%d.dat'%(tobs,window,window))
      #x = np.bincount(da[:10000,-1].astype(int),minlength=2).astype(float)
      x = np.bincount(da[int(i*10000):int((i+1)*10000),-1].astype(int),minlength=2).astype(float)
      x /= 10000.
      y = -np.log(x[x!=0])
      data_list.append( y.copy() )
    else:
      data_list.append( np.array([0.]) )
     
  cs = np.zeros(len(windows))
  res = minimize(opt_func,cs,args=(data_list,windows))
  for j in range(len(windows)):
    data_list[j][:] += res.x[j]
  
  f = open('tobs%d/reweighted_hist_%d.dat'%(tobs,i),'w')
  #f = open('tobs%d/reweighted_hist.dat'%(tobs),'w')
  for k in range(len(windows)):
    for j in range(len(data_list[k])):
      f.write('%d %.8f\n'%(window_states[k][j],data_list[k][j]))
  f.close()
