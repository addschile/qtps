import numpy as np
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

da = np.loadtxt('path_ensemble_ratios.dat')
plt.plot(da[:,0],da[:,1]*1.e4,'sb',markersize=15,mec='k')
plt.errorbar(da[:,0],da[:,1]*1.e4,yerr=da[:,2]*1.e4,fmt='none',ecolor='k',capsize=10)
plt.xticks(size=15)
plt.yticks(size=15)
plt.xlabel(r'$t$ / a.u.',size=20)
plt.ylabel(r'$Z_{AB} (t) / Z_A$ / $\times 10^4$',size=20)
plt.tight_layout()
plt.savefig('path_ensemble_ratios.png')
plt.show()
