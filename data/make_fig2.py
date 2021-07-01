import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gspc
from scipy.optimize import curve_fit

def func(t,k):
    return k*t

fs_2_au = 41.3413745758 # au/fs

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

da = np.loadtxt('tobs200/avg_reweighted_hist.dat')
plt.plot(da[:,0],(da[:,1]/200.)*fs_2_au*1000,'sb',markersize=13,mec='k')
#plt.errorbar(da[:,0],da[:,1]/200.,yerr=1.812*da[:,2]/200.,fmt='none',ecolor='k',capsize=5)
plt.ylabel(r'$t^{-1} Z_{Aq} / Z_{A}$ / ps$^{-1}$',size=20)
plt.xlabel(r'$\langle q \rangle$ / $\mathrm{\AA}$',size=20)
plt.yticks(size=15)
plt.xticks(size=15)
plt.yscale('log')

da = np.loadtxt('path_ensemble_ratios.dat')
#print(da[:,1]/(da[:,0]/fs_2_au/1000.))
#print(da[:,2]/(da[:,0]/fs_2_au/1000.))
#k = da[:,1]/(da[:,0]/fs_2_au/1000.)
#print(np.sum(k)/4.)
#kerr = da[:,2]/(da[:,0]/fs_2_au/1000.)
#print(np.sqrt(np.sum(kerr**2.)/4.))
plt.axes([.42,.61,.55,.350])
plt.plot(da[:,0]/fs_2_au/1000.,da[:,1]*1.e4,'sb',markersize=8,mec='k')
plt.errorbar(da[:,0]/fs_2_au/1000.,da[:,1]*1.e4,yerr=da[:,2]*1.e4,fmt='none',ecolor='k',capsize=7)
popt,pocv = curve_fit(func,da[:,0],da[:,1])
plt.plot(da[:,0]/fs_2_au/1000,func(da[:,0],popt[0])*1.e4,'-k',lw=1)
plt.xticks(size=10)
plt.yticks(size=10)
plt.xlabel(r'$t$ / ps',size=12)
plt.ylabel(r'$Z_{AB} (t) / Z_A$ / $\times 10^4$',size=12)
plt.tight_layout()
plt.savefig('rates_from_tps.pdf')
plt.show()

