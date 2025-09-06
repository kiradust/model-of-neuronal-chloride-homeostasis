from plm_singlecomp_withkcc2 import plm,zplm,F
import matplotlib.pyplot as plt
from plotting import clcolor, kcolor, nacolor, wcolor, minithreefig
from matplotlib import gridspec
import numpy as np

pre_gaba = -80.9 

pre_vm_pair_vu = -68.5
pre_vm_unpair = -71.5 

post_gaba_pair_vu = -70.5 
post_gaba_unpair_vu = -69.1 
post_gaba_chabc = -64.5 

post_vm_pair_vu = -69.1 
post_vm_unpair_vu = -72.1 
post_vm_chabc = -56.7 

molinit=plm(gx=1e-8,xt=25,tt=100,two=1,paratwo=True,moldelt=0)
para=zplm(molinit=molinit)

argvm = np.argmin(np.abs(np.array(para[5])-pre_vm_unpair))
argecl = np.argmin(np.abs(np.array(para[3])-pre_gaba))
print(argvm,argecl)
print('evm'+str(para[5][argvm]),'ecl'+str(para[3][argvm]),'ena'+str(para[1][argvm]),'ek'+str(para[2][argvm]))
print('cli'+str(para[8][argvm]),'nai'+str(para[6][argvm]),'ki'+str(para[7][argvm]),'xi'+str(para[9][argvm]))
print(para[5][argecl],para[3][argecl])
print(para[0][argvm],para[0][argecl])

# plotting
gs = gridspec.GridSpec(3, 1, height_ratios=[1.5, 1, 1]) 
plt.subplot(gs[0])
plt.plot(para[0],para[3],color=clcolor,linestyle='-',label='ecl')
plt.plot(para[0],para[2],color=kcolor,linestyle='-',label='ek')
plt.plot(para[0],para[5],'k-',label='vm')
plt.legend()
plt.ylabel('mV')
plt.vlines(para[0][argvm],-100,0,linestyles='--',colors='gray')
plt.subplot(gs[1])
plt.plot(para[0],para[1],color=nacolor,linestyle='-')
plt.ylabel('ena, mV')
plt.vlines(para[0][argvm],0,150,linestyles='--',colors='gray')
plt.subplot(gs[2])
plt.plot(para[0],para[11],color=wcolor,linestyle='-')
plt.vlines(para[0][argvm],1e-12,1e-11,linestyles='--',colors='gray')
plt.ylabel('volume')
plt.xlabel('pump rate')
plt.savefig('default_values.png',dpi=150)
plt.show()

q=10**(para[0][argvm]/1000.0)/F
dg=plm(p=q,tt=1000)
a1,a2,a3=minithreefig([dg[11][4:-1],dg[14][4:-1],dg[13][4:-1],dg[16][4:-1],dg[10][4:-1],dg[24][4:-1],dg[12][4:-1]],'k',yl=[[-100,-70],[1.92e-12,1.98e-12],[0,6e-7]])
plt.savefig('eg_run_z.png',dpi=150)
plt.show()