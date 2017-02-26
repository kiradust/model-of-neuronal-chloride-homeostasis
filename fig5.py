from plm_singlecomp_withkcc2 import plm, zp, F
from plotting import minifig, minithreefig, minifigtwoaxes,twoaxes,xcolor,clcolor,wcolor,kcolor
import matplotlib.pyplot as plt
from pylab import rcParams
from fig4 import f4e
import numpy as np
rcParams['figure.figsize'] = 8,8
from matplotlib import gridspec

sym=['-b',':r','--g','-.m']

def f5a():
    dez=plm(dz=1e-6,two=1,xt=30,tt=180)
    minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
    return

def f5b(moldelt=0):
    XF=[-1.2,-1.15,-1.1,-1.05,-1.0,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5]
    XFp=[]
    cl=[]
    vm=[]
    k=[]
    w=[]
    z=[]
    x=[]
    ax0, ax1, ax2, pi, zi, zee, newx = f4e(moldelt=moldelt)
    
    for i in range(len(zi)):
        for b in XF:
            if zi[i]/100.0==b:
                XFp.append(pi[i])
                print pi[i]
    for a in range(len(XF)):
        #dez=plm(two=1,xt=10,dz=3e-6,Zx=a,tt=5000)
        dez=plm(z=XF[a],tt=300,p=10**(XFp[a])/F)
        cl.append(dez[7])
        k.append(dez[6])
        vm.append(dez[9])
        w.append(dez[10][-1])
        z.append(dez[22][-1]*100)
        x.append(dez[3])
    ax0, ax1, ax2, pi, zi, zee, newx = f4e(moldelt=moldelt)
    ax0.plot(z,cl,'bo')
    ax0.plot(z,k,'go')
    ax0.plot(z,vm,'ko')
    ax2.plot(z,x,'mo')
    ax1.plot(z,w,'ko')
    plt.show()
    
    plt.figure()
    gs = gridspec.GridSpec(3,1,height_ratios=[1.5,0.5,0.5])
    plt.subplot(gs[0])
    plt.plot(zi,zee[3],color=clcolor)
    plt.plot(z,cl,'bo')
    plt.plot(z,k,'go')
    plt.plot(z,vm,'ko')
    plt.plot(zi,zee[2],color=kcolor)
    plt.plot(zi,zee[5],'k')
    plt.subplot(gs[1])
    plt.plot(zi,zee[11],color='k') #vol
    plt.plot(z,w,'ko')
    plt.subplot(gs[2])
    plt.plot(zi,newx,color='m') #conc X
    plt.plot(z,x,'mo')
    plt.show()
    return

def f5c(ratio=0.8,md=1e-12):
    dez=plm(gx=1e-8,xt=20,two=1,xend=0,moldelt=md,xflux=0.3*1e-6,ztarget=-0.9,tt=300)
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
            dez=plm(xend=0,two=1,xt=10,xflux=1e-6,ztarget=b,Zx=a,moldelt=1,ratio=0.1)
            #minithreefig([dez[11][1:-1],dez[14][1:-1],dez[13][1:-1],dez[16][1:-1],dez[10][1:-1],dez[22][1:-1]],'k')
            w[abs(int(a))].append(np.log10(dez[10][-1]))
            z[abs(int(a))].append(dez[22][-1])
        plt.plot(z[abs(int(a))],w[abs(int(a))],sym[abs(int(a))]+'o')
    plt.show()
    return w,z