# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 13:48:10 2021

@author: Brown
"""
# Load necessary modules
import sys
import numpy as np
from scipy import interpolate
from scipy import integrate
#-------------------------------------------------------------------

#********************************************************************
# python function to compute the first derivative at the first point
# and the last point of an array. For this computation, the first 4
# points are used for the derivative for the first data point. Similarly
# last 4 (x,y) points are used for the derivative at the last point.

def fd_ends(x,y):
    ''' Parameter:
        x       =       xaxis of the array
        y       =       y data of the array
                (x and y arrays must have atleast 4 elements)
        Returns = first derivative at the first point and the last point
    '''
    if (len(x) < 4 or len(y) < 4):
        print("Error : x and y arrays must have 4 elements")

    subx=np.zeros(4)
    suby=np.zeros(4)

    for i in range(0,4):
        subx[i] = x[i]
        suby[i] = y[i]

    fd1 = ((subx[1]*subx[3]+subx[2]*subx[3]+subx[1]*subx[2]-2*subx[0]*subx[1] \
          -2*subx[0]*subx[3]-2*subx[0]*subx[2]+3*subx[0]**2)/(-subx[3]+subx[0]) \
          /(-subx[1]+subx[0])/(-subx[2]+subx[0])*suby[0]-(-subx[2]+subx[0])    \
          *(-subx[3]+subx[0])/(subx[1]-subx[3])/(-subx[1]+subx[0])/(subx[1]-subx[2])\
          *suby[1]+(-subx[1]+subx[0])*(-subx[3]+subx[0])/(subx[2]-subx[3])          \
          /(subx[1]-subx[2])/(-subx[2]+subx[0])*suby[2]-(-subx[1]+subx[0])          \
          *(-subx[2]+subx[0])/(subx[2]-subx[3])/(subx[1]-subx[3])/(-subx[3]+subx[0])\
          *suby[3] )

    for i in range(0,4):
        subx[i] = x[int(i-4)]
        suby[i] = y[int(i-4)]
#        print (i, int(i-4))

    fdn = ( (subx[1]-subx[3])*(subx[2]-subx[3])/(-subx[3]+subx[0])/(-subx[1] \
           +subx[0])/(-subx[2]+subx[0])*suby[0]-(-subx[3]+subx[0])*(subx[2]   \
           -subx[3])/(subx[1]-subx[3])/(-subx[1]+subx[0])/(subx[1]-subx[2])   \
           *suby[1]+(-subx[3]+subx[0])*(subx[1]-subx[3])/(subx[2]-subx[3])    \
           /(subx[1]-subx[2])/(-subx[2]+subx[0])*suby[2]-(-2*subx[0]*subx[3]  \
           -2*subx[1]*subx[3]-2*subx[2]*subx[3]+subx[0]*subx[1]+subx[0]       \
           *subx[2]+subx[1]*subx[2]+3*subx[3]**2)/(subx[2]-subx[3])/(subx[1]   \
           -subx[3])/(-subx[3]+subx[0])*suby[3]  )

    return(fd1,fdn)

#********************************************************************

# Spline function taken from Numerical Recipes in FORTRAN, page 109 and 110
# Numerical recipes in FORTRAN, Second Edition, Press, Teukolsky, Vetterling, Flannery
# Cambridge University Press, 1992
def spline (x, y, yp1, ypn):

    '''Parameters:
        x       =       1D vector of x-values in increasing order.
        y       =       1D vector of y-values
        n       =       number of elements in xVector (and yVector)
        yp1     =       first derivative of the interpolating function at the first segment
        ypn     =       first derivative of the interpolating function at the last segment
    '''
    nx=len(x)
    ny=len(y)

    if (nx == ny):
        n=nx
    else :
        print("Error : x and y data have different lengths in spline.")
        quit()

    u=np.zeros(n)
    y2=np.zeros(n) # this is the output
    p = 0.0
    sig = 0.0

    if yp1 > 1e30 :     # lower boundar condition 'natural'
        y2[0] = 0.0
        u[0]  = 0.0
    else :              # specified first derivative
        y2[0] = -0.5
        u[0]=(3/(x[1]-x[0])) * ( ((y[1]-y[0])/ (x[1]-x[0])) - yp1)

    for i in range(1,n-1) :     #       Decomposition loop of tridiagonal algorithm. y2 and u are temporary.
        sig = ( (x[i]-x[i-1])/(x[i+1]-x[i-1]) )
        p = (sig * y2[i-1]) + 2.0
        y2[i] = (sig-1.0)/p
        u[i] = ((6.0*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1]) / (x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p)
        # print("first loop:",i)

    if ypn > 1e30 :     # upper boundary condition 'natural'
        qn=0.0
        un=0.0
    else :              # specified first derivative
        qn=0.5
        un = (3.0 /(x[n-1]-x[n-2]))*(ypn-((y[n-1]-y[n-2])/(x[n-1]-x[n-2])) )

    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0)

    for k in range(n-2,-1,-1):           # from second last point to the second point
        y2[k]=y2[k]*y2[k+1]+u[k]        # backsubstitution loop of tridiagonal algorithm
        # print("loop 2 :",k)

    return(y2)


#************************************************************************
# Spline interpolation function taken from Numerical Recipes in FORTRAN, page 109 and 110
# Numerical recipes in FORTRAN, Second Edition, Press, Teukolsky, Vetterling, Flannery
# Cambridge University Press, 1992

def splint(xa, ya, y2a, x):
    ''' Parameters :
        xa      =       original x-axis 1D vector
        ya      =       original y-axis 1D vector
        y2a     =       output of the spline function
        x       =       new x axis, scalar
    '''

    nxa=len(xa)
    nya=len(ya)
    ny2a=len(y2a)

    if (nxa != nya or nxa != ny2a or nya != ny2a):
        print("Error : xa or ya or y2a have incorrect dimension(s).")
        quit()

    n = nxa

    klo = int(0)
    khi = int(n-1)
    k = int(0)
    h = 0.0
    element=0.0

    while ((khi-klo) > 1) :
        k = int((khi+klo)/2)
        element=xa[k]
#        print(element,xa[k],k,x)
        if ( element > x) :
            khi = k
        else:
            klo = k

    h = xa[khi] - xa[klo]

    if h == 0 :
        print("Error : Bad xa input in splint")
        quit()

    a = (xa[khi]-x)/h
    b = (x-xa[klo])/h
    y = a*ya[klo]+b*ya[khi] + (( (a**3)-a)*y2a[klo]+( (b**3)-b)*y2a[khi])*(h**2)/6.0

    return(y) #returns the interpolated value

#************************************************************************


#************************************************************************

# Define the wrappend intergral function which computes the matrix element by evaluating the
#   the numerical integral

# y_parameter should be already interpolated to rwave to that the product may be computed

def compute_int(rwave, psi1, psi2, y_parameter, rMin, rMax):
    p1 = np.multiply(psi1, psi2)
    p2 = np.multiply(p1, rwave)
    p3 = np.multiply(p2, rwave)
    product = np.multiply(p3, y_parameter)

    derivative = fd_ends(rwave, product)
    secarray2 = spline(rwave, product, derivative[0], derivative[1])

    # function defining the integrand which uses the spline coef array to give interpolated values
    def integrand_ME(xpoint):
        '''integrand for the quadrature'''
        result = splint(rwave, product, secarray2, xpoint)
        return result

    res = integrate.quadrature(integrand_ME, rMin, rMax, tol=1.0e-7, vec_func=False, maxiter=1000)
    return res
    # output is a tuple, [0] = integral result, [1]=error

#************************************************************************