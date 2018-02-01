"""
Last modified 26 January 2018

@author: Kira
"""
from pylab import polyfit
import matplotlib.pyplot as plt
from scipy.stats import linregress
import statsmodels.api as sm
import numpy as np

print "\nFigure 3C"

# a - kcc2 change
a=[-0.35,-0.42,-0.37,-0.98,-0.04,-0.17,-0.75,-0.45,-0.20]

# b - chloride change
b=[0.56,0.72,0.59,1.18,0.02,0.37,1.51,0.42,0.20]

# c - DF change (in mV)
c=[-11.20,-12.00,-9.90,-20.90,-8.45,-7.07,-24.60,-9.00]
v=[0.2,2,2,-0.9,-8,1,-1,0]
df = True # set df to True to change y-variable to DF and include Vm change plot, or false to make it for Cl change only

# study names
names=['Tang et al. 2015','Lee et al. 2011','Campbell et al. 2015','Lagostena et al. 2010',
       'MacKenzie et al. 2014','MacKenzie et al. 2014','Coull et al. 2003','Ferrini et al. 2013','Mahadevan et al. 2015']

# weights
w=[31,35,34,36,23,24,18,19,20]

#labels
lab_age_no_path = ['1-3 months', 'less than 1 month','1-3 months','more than 6 months',
           '1-3 months', '1-3 months', '1-3 months', '1-3 months',  
           'more than 6 months']
lab_area_no_path =['spinal cord lamina I','hippocampus','neocortex','hippocampus',
     'hippocampus','hippocampus','spinal cord lamina I','spinal cord lamina I', 'hippocampus']
nkcc1 = [1,1,1,1,0,0,0,0,0,0]
disease = ['stress','excitotoxicity','glioma','anti-NGF','stress',
           'stress','pain','pain','neto2 knockout']

# rough colour output script (matches label chosen in labels to colour ouput; prints legend)
labels = nkcc1
label_use = []
color_use = ['r','b','g','m','k','c']
colors = []

for m in labels:
    true = False
    for i in range(len(label_use)):
        if m==label_use[i]:
            colors.append(color_use[i])
            true = True
    if true == False:
        label_use.append(m)
        colors.append(color_use[label_use.index(m)])

print label_use, color_use
ax = plt.subplot(111)

ylab = "Percentage change in intracellular chloride concentration"

def regplot(x=a,y=b,an=[],bn=[]):
    if len(y) < len(x):
        x.pop()
        if len(names) > len(y):
            names.pop()
            w.pop()
            ylab = "Change in Driving Force (mV)"

    for i in xrange(len(w)):
        for z in xrange(w[i]):
            an.append(x[i])
            bn.append(y[i])

    # linregress (linear regression analysis)
    # Parameters:	
    # x, y : array_like two sets of measurements. Both arrays should have the same length. If only x is given (and y=None), then it must be a two-dimensional array where one dimension has length 2. The two sets of measurements are then found by splitting the array along the length-2 dimension.
    # Returns:	
    # slope : float; slope of the regression line
    # intercept : float; intercept of the regression line
    # rvalue : float; correlation coefficient
    # pvalue : float; two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero
    # stderr : float; Standard error of the estimate
    print linregress(an,bn)

    # scipy's OLS reg analysis
    print reg_m(bn, an).summary()

    # p = polyfit(x,y,n) returns the coefficients for a polynomial p(x) of degree n that is a best fit (in a least-squares sense) for the data in y. The coefficients in p are in descending powers, and the length of p is n+1
    m,c = polyfit(an, bn, 1)

    # weighting depiction in scatter plot and line-of-best-fit plot
    s = [20*(n-10) for n in w]
    plt.scatter(x,y,s=s,c=colors)
    plt.plot(x, np.linalg.linalg.multiply(m,x)+c, color='black')

    # add labels to plot
    #for i, xy in enumerate(zip(x, y)):
        #ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='offset points')
        #ax.annotate(names[i], xy=xy, textcoords='offset points')

    return an,bn

# run scipy's OLS (ordinary least squares) regression analysis
def reg_m(y, x):
    ones = np.ones(len(x))
    X = sm.add_constant(np.column_stack((x, ones)))
    results = sm.OLS(y, X).fit()
    return results