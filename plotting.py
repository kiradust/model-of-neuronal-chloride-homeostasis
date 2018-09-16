# -*- coding: utf-8 -*-
"""
Created on Sun Nov 06 13:34:52 2016

@author: Kira
"""
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pylab import rcParams
rcParams['figure.figsize'] = 8,8

# colours chosen from colorbrewer2.org
clcolor='#1b9e77'
kcolor='#7570b3'
xcolor='#d95f02'
nacolor='#e7298a'
wcolor='k'

ys = [1,1.5,2,2.5,5]
xs = range(len(ys))

def minmax(data):
    return min(data), max(data)
    
def twoaxes(x,y11,y12,y22,x2=None):
    if x2==None:
        x2=x
    fig, ax1 = plt.subplots()
    ax1.plot(x,y11,color='k',linestyle='-')
    ax1.plot(x2,y12,color='k',linestyle='--')
    ax2 = ax1.twinx()
    ax2.plot(x,y22,color=nacolor)
    plt.show()
    return
    
def minifig(delta,x=0, yl=[[-100,40],[1.75e-12,2.75e-12]]):
    plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[1.5, 1]) 
    ax0=plt.subplot(gs[0])
    ax0.plot(delta[0],delta[1],color=clcolor)
    ax0.plot(delta[0],delta[2],color=kcolor)
    ax0.plot(delta[0],delta[5],'k')
    #ax0.plot(delta[0],delta[3],color=nacolor)
    ax0.set_ylim(yl[0])
    if x!=0:
        ax0.axvline(x=x,linestyle='--',color='0.8')
    ax1=plt.subplot(gs[1])
    ax1.plot(delta[0],delta[6],color=wcolor)
    #ax1.set_ylim(yl[1])
    if x!=0:
        ax1.axvline(x=x,linestyle='--',color='0.8')
    return
    
def minifigtwoaxes(delta):
    plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[1.5, 1]) 
    ax0=plt.subplot(gs[0])
    ax0.plot(delta[0],delta[1],color=clcolor,clip_on=False)
    ax0.plot(delta[0],delta[2],color=kcolor,clip_on=False)
    ax0.plot(delta[0],delta[3],'k',clip_on=False)
    ax1=plt.subplot(gs[1])
    #ax0.get_shared_y_axes().join(ax0, ax1)
    ax1.plot(delta[0],delta[4],color=wcolor,clip_on=False)
    ax2=ax1.twinx()
    ax2.plot(delta[0],delta[5],color=xcolor,clip_on=False)
    return ax0, ax1, ax2
    
def minithreefig(delta,colour,x=0,yl=[[-100,-70],[1.0e-13,1.6e-13],[-0.95,-0.8]]):
    plt.figure()
    gs = gridspec.GridSpec(3,1,height_ratios=[1.5,0.5,0.5])
    ax0=plt.subplot(gs[0])
    ax0.plot(delta[0],delta[1],color=clcolor)
    ax0.plot(delta[0],delta[2],color=kcolor)
    ax0.plot(delta[0],delta[3],'k')
    ax0.set_ylim(yl[0])
    if x!=0:
        ax0.axvline(x=x,linestyle='--',color='0.8')
    ax1=plt.subplot(gs[1])
    ax1.plot(delta[0],delta[4],color=wcolor) #volume
    ax1.set_ylim(yl[1])
    if x!=0:
        ax1.axvline(x=x,linestyle='--',color='0.8')
    ax2=plt.subplot(gs[2])
    ax2.plot(delta[0],delta[5],color=colour) #conc X
    if x!=0:
        ax2.axvline(x=x,linestyle='--',color='0.8')
    ax2.set_ylim(yl[2])
    return ax0, ax1, ax2


def fluxplot(timevalues,title=''):
    plt.figure()
    gs = gridspec.GridSpec(5,1,height_ratios=[1, 1, 1, 0.5, 0.5])
    
    ax0=plt.subplot(gs[0])
    ax0.plot(timevalues[11],timevalues[28],color=nacolor)
    ax0.plot(timevalues[11],timevalues[29],color=kcolor)
    ax0.plot(timevalues[11],timevalues[30],color=clcolor)
    ax0.plot(timevalues[11],timevalues[32],color=xcolor)
    ax0.set_ylabel("ion fluxes")
    
    ax1=plt.subplot(gs[1])
    ax1.plot(timevalues[11],timevalues[18],color=nacolor)
    ax1.plot(timevalues[11],timevalues[19],color=kcolor)
    ax1.plot(timevalues[11],timevalues[17],color=clcolor)
    ax1.plot(timevalues[11],timevalues[20],color=xcolor)
    ax1.set_ylabel("concentrations")
    
    ax2=plt.subplot(gs[2])
    ax2.plot(timevalues[11],timevalues[12],color=nacolor)
    ax2.plot(timevalues[11],timevalues[13],color=kcolor)
    ax2.plot(timevalues[11],timevalues[14],color=clcolor)
    ax2.plot(timevalues[11],timevalues[15],color=xcolor)
    ax2.plot(timevalues[11],timevalues[16],color='k')
    ax2.set_ylabel("potentials")
    
    # ws on own
    ax3=plt.subplot(gs[3])
    ax3.plot(timevalues[11],timevalues[31],color=wcolor)
    ax3.set_ylabel("volume flux")
    
    ax4=plt.subplot(gs[4])
    ax4.plot(timevalues[11],timevalues[10],color=wcolor)
    ax4.set_ylabel("volume")
    ax4.set_xlabel("time in seconds")
    if title !='':
        plt.savefig(title)
    plt.show()
    return