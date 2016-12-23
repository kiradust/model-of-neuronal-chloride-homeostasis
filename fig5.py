from plm_singlecomp_withkcc2 import plm, zp
from plotting import minifig, minithreefig, minifigtwoaxes,twoaxes,xcolor,clcolor,wcolor,kcolor
import matplotlib.pyplot as plt
from pylab import rcParams
from fig4 import f4e
import numpy as np
rcParams['figure.figsize'] = 8,8

sym=['-b',':r','--g','-.m']

def f5a():
    dez=plm(dz=1e-5,two=1,xt=10)
    minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
    return

def f5b(moldelt=0):
    #XF=[-1.2,-1.1,-1.0,-0.9,-0.85,-0.8,-0.7,-0.6,-0.55]
    DZ=[0e-6,1e-6,2e-6,3e-6,4e-6]
    cl=[]
    vm=[]
    k=[]
    w=[]
    z=[]
    x=[]
    for a in DZ:
        dez=plm(two=1,xt=10,dz=a,Zx=-1.19)
        #minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
        cl.append(dez[7])
        k.append(dez[6])
        vm.append(dez[9])
        w.append(dez[10][-1])
        z.append(dez[22][-1]*100)
        x.append(dez[3])
        dez=plm(two=1,xt=10,dz=a,Zx=-0.51)
        #minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
        cl.append(dez[7])
        k.append(dez[6])
        vm.append(dez[9])
        w.append(dez[10][-1])
        z.append(dez[22][-1]*100)
        x.append(dez[3])
    ax0, ax1, ax2 = f4e(moldelt=moldelt)
    ax0.plot(z,cl,'bo')
    ax0.plot(z,k,'go')
    ax0.plot(z,vm,'ko')
    ax2.plot(z,x,'mo')
    ax1.plot(z,w,'ko')
    plt.show()
    return

def f5c(ratio=0.8,md=1e-13):
    #dez=plm(gx=1e-8,xt=25,tt=200,two=1,xend=0,moldelt=md,ratio=ratio)
    dez=plm(xend=0,two=1,xt=10,xflux=5e-6,moldelt=9e-12)
    minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
    return

def f5d():
    w=[[],[],[],[]]
    z=[[],[],[],[]]
    ZX=[-0.5,-1,-2,-3]
    ZT=[-0.5,-0.55,-0.6,-0.85,-0.9,-0.95,-1,-1.5,-1.75,-1.9,-2.0,-2.5,-2.75,-2.9,-3.0]
    plt.figure()
    for a in ZX:
        for b in ZT:
            dez=plm(xend=0,two=1,xt=10,xflux=5e-6,ztarget=b,Zx=a,moldelt=1,ratio=0.1)
            #minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
            w[abs(int(a))].append(np.log10(dez[10][-1]))
            z[abs(int(a))].append(dez[22][-1])
        plt.plot(z[abs(int(a))],w[abs(int(a))],sym[abs(int(a))]+'o')
    plt.show()
    return w,z