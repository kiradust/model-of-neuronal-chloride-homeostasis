from pylab import polyfit
import matplotlib.pyplot as plt
from scipy.stats import linregress
import statsmodels.api as sm

a = [0.35, 0, 0.42, 0.3720930233, 0.9733333333]
b = [0.5584754086, -0.2071144338, 1.860188312, 0.5891017512, 1.35438295]

print linregress(a,b)

# Parameters:	
# x, y : array_like two sets of measurements. Both arrays should have the same length. If only x is given (and y=None), then it must be a two-dimensional array where one dimension has length 2. The two sets of measurements are then found by splitting the array along the length-2 dimension.
# Returns:	
# slope : float; slope of the regression line
# intercept : float; intercept of the regression line
# rvalue : float; correlation coefficient
# pvalue : float; two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero
# stderr : float; Standard error of the estimate

m,c = polyfit(a, b, 1)

# p = polyfit(x,y,n) returns the coefficients for a polynomial p(x) of degree n that is a best fit (in a least-squares sense) for the data in y. The coefficients in p are in descending powers, and the length of p is n+1

labels = ['25792562','19252497','21532577','25066727','20089897']

plt.plot(a[0], b[0], 'ro',label='Whole cell patch clamp')
plt.plot(a[1], b[1], 'ro')
plt.plot(a[2], b[2], 'bo',label='Gramcidin perforated patch clamp')
plt.plot(a[3], b[3], 'yo',label='Fluorescence imaging')
plt.plot(a[4], b[4], 'go',label='Cell-attached patch clamp')
plt.plot(a, m*a+c, color='black')

plt.legend(loc='lower right')

plt.xlabel("Percentage change in KCC2 expression")
plt.ylabel("Percentage change in intracellular chloride concentration")

plt.show()

def reg_m(y, x):
    ones = np.ones(len(x))
    X = sm.add_constant(np.column_stack((x, ones)))
    results = sm.OLS(y, X).fit()
    return results
    
print reg_m(b, a).summary()