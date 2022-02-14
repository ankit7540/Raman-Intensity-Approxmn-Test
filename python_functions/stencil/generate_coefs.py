# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 19:22:00 2021

@author: Brown
"""
import numpy as np
import scipy as sp
import math


###########################################################################
def coef(n,m):
    """
    

    Parameters
    ----------
    n : scalar
        derivative order
    m : scalar
        m-order of the stencil approximation

    Returns
    -------
    coef_matrix : 2D numpy array
        coeffcients for approximating the derivative

    """
    #define offset matrix
    offset_mat= np.zeros((m,m),dtype=int)
    for i in range(m):
        for j in range(m-1,-1,-1):
            offset_mat[i,j]=j-i
    
    
    coef_matrix= np.zeros((m,m),dtype=float)
    a=np.zeros((m,m),dtype=float)
    b=np.zeros(m,dtype=int)
    b[n]=1
    for k in range(m):
        c=np.zeros(m,dtype=float)
        for i in range(m):
            a[:,i]=i
            for j in range(m):
                a[j,i]=offset_mat[k,i]**j  /math.factorial(j)
        c=np.dot(np.linalg.inv(a),b)
        #print(c)
        for s in range(m):
            coef_matrix[k,s]=c[s]
    return coef_matrix

###########################################################################