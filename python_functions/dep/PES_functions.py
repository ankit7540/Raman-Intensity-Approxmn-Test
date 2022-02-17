
import numpy as np
import math
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