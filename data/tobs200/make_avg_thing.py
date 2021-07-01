import numpy as np

Z = np.zeros((10,8))
for i in range(10):
    da = np.loadtxt('reweighted_hist_%d.dat'%(i))
    for j in range(13):
        if da[j,0] < 8:
            Z[i,int(da[j,0])] = np.exp(-da[j,1])
print(Z[0,:])

for i in range(10):
    Z[i,:] /= (Z[i,0]+Z[i,2])
    #for j in range(8):
    #    Z[i,j] /= (Z[i,0]+Z[i,2])
print('')
print(Z[0,:])

Zavg = np.zeros(8)
for i in range(8):
    Zavg[i] = np.sum(Z[:,i])/10.
for i in range(10):
    for j in range(8):
        Z[i,j] -= Zavg[j]
Z *= Z
std_err = np.zeros(8)
for i in range(8):
    std_err[i] = np.sqrt(np.sum(Z[:,i])/10.)

f = open('avg_reweighted_hist.dat','w')
for i in range(8):
    f.write('%d %.16f %.16f\n'%(i,Zavg[i],std_err[i]))
f.close()
