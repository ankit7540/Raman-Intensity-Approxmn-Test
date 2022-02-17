# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt

from dep import diatomicSE
from dep import transition

import scipy.optimize as opt
import math

######################################################

# load ref data here

N2_ref=np.loadtxt('N2_ref_data/14N15N_ref_J4_v1.txt')

re=2.068

param_opt=np.array([0.71696534, 1.00934466, 2.07636687])

######################################################

def morse_potential(param, xaxis):
    '''
    Morse potential function which takes two parameters
    param : np array of three elements(De, a , re)
    xaxis : 1D array of distance on which the Morese morse_potential
            is to be computed
    '''

    power = -param[1]*(xaxis-param[2])
    return param[0]* (1-np.exp(power))**2


######################################################

def rydberg_potential_3param(param, xaxis):

    # param[0] = De
    # param[1] = ke
    # param[2] = re

    b = math.sqrt(param[1]/param[0])
    power = b*(xaxis-param[2])

    return -param[0]* np.exp(-power)*(1+power)+param[0]



######################################################
def modification_morse_potential(param,xaxis):

    # param[0] = De
    # param[1] = a
    # param[2] = b
    # param[3] = c
    # param[4] = re

    power = -param[1]*(xaxis-param[4])+param[2]*(xaxis-param[4])**2+param[3]*(xaxis-param[4])**3

    return param[0]* (1-np.exp(power))**2


#######################################################
#   optimize function ( optimization of spline coefs)
def run_fit(initial_guess):

    res = opt.minimize( residual_rovibration_3v, initial_guess,\
                       method='Nelder-Mead', \
                              options={'xatol': 1e-8,\
                                       'fatol': 1e-8,\
                                           'maxiter':2500})

    #print (res)
    #opt_params = res.x

    return res



######################################################


def residual_rovibration_3v (parameter):


    xfit = np.arange(0.75, 5.99 , 0.05)
    yfit = morse_potential(parameter,xfit)
    step = 0.010
    Jmax = 4


    e = diatomicSE.get_eigenvalue_J( 7, 14, 7, 15, xfit , yfit, 5, step , 0 )
    dat = np.zeros((e.shape[0], Jmax+1))
    dat[:, 0]=e[:,0]

    for i in range(1, Jmax+1):
        e = diatomicSE.get_eigenvalue_J( 7, 14, 7, 15, xfit , yfit, 5, step , i )

        dat[:, i]=e[:,0]


   #P_v0v0= transition.P_branch(dat, 0, 0, Jmax)
    S_v0v0= transition.S_branch(dat, 0, 0, Jmax)

    O_v0v1= transition.O_branch(dat, 0, 1, Jmax)
    Q_v0v1= transition.Q_branch(dat, 0, 1, Jmax)
    S_v0v1= transition.S_branch(dat, 0, 1, Jmax)


    calc = np.hstack((S_v0v0,O_v0v1,Q_v0v1,S_v0v1))
    #print(calc.shape)

    # subtract with ref
    ref = N2_ref

    #print (calc.shape, ref.shape)
    #print(ref,'\n\n', calc)
    diff = calc - ref
    print(np.sum(diff**2))

    return np.sum(diff**2)


######################################################

def gen_result(xfit,yfit):

    #xfit = np.arange(0.75, 5.99 , 0.05)
    #yfit = rydberg_potential_3param(parameter,xfit)
    step = 0.010
    Jmax = 4

    e = diatomicSE.get_eigenvalue_J( 7, 14, 7, 15, xfit , yfit, 5, step , 0 )
    dat = np.zeros((e.shape[0], Jmax+1))
    dat[:, 0]=e[:,0]

    for i in range(1, Jmax+1):
        e = diatomicSE.get_eigenvalue_J( 7, 14, 7, 15, xfit , yfit, 5, step , i )

        dat[:, i]=e[:,0]


    out = np.zeros((13, 3))

    S_v0v0= transition.S_branch(dat, 0, 0, Jmax)

    O_v0v1= transition.O_branch(dat, 0, 1, Jmax)
    Q_v0v1= transition.Q_branch(dat, 0, 1, Jmax)
    S_v0v1= transition.S_branch(dat, 0, 1, Jmax)





    calc = np.hstack((S_v0v0,O_v0v1,Q_v0v1,S_v0v1))
    #print(calc.shape)

    # subtract with ref
    ref = N2_ref

    #print (calc.shape, ref.shape)

    diff = calc - ref


    for i in range(13):
        out[i,0]=ref[i]
        out[i,1]=calc[i]
        out[i,2]=diff[i]

    return out



######################################################
