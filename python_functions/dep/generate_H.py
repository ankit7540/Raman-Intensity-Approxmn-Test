# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 22:11:11 2021

@author: Ankit Raj

"""

import periodictable as pt
import numpy as np
from scipy import interpolate
from scipy import constants


from stencil import generate_coefs
import matplotlib.pyplot as plt
#-------------------------------------------------------------------

# temporary code to load files



#-------------------------------------------------------------------

def reduced_mass(zA, iMassA, zB, iMassB ):
    '''
    zA = atomic number for atom A
    iMassA = atomic mass of specific isotope
    zB = atomic number for atom B
    iMassB = atomic mass of specific isotope
    
    reduced_mass( 1, 1, 1, 2) for HD
    reduced_mass( 6, 12, 8, 16 ) for CO
    
    
    '''
    mass = 1822.8884862086923 # amu to  electron mass

    sA=pt.elements[zA][iMassA]
    sB=pt.elements[zB][iMassB]
    
    #print(sA, sB)

    mA=float(sA.mass)
    mB=float(sB.mass)
    
    A=mA*mass
    B=mB*mass

    #reduced mass , atomic unit
    return 1/(1/A + 1/B)

#----------------------------------------------------------

def gen_H_matrix(reduced_mass, rwave , potential, stencil_number, step):
    '''
    Generate the Hamiltonian matrix for the radial
    nuclear equation for diatomic molecule.
    
    Numerical derivatives are used to approximate the 
    partial derivative in the radial nuclear eqn.
    '''
    
    
    print ("reduced mass : ", reduced_mass)
    
    
    # ---- generate full r vector ----------
    rMin = np.amin(rwave)
    rMax = np.amax(rwave)
    
    rFull = np.arange(rMin, rMax, step, dtype=float)
    
    # ---- interpolate potential -----------
    interpolation_data = interpolate.interp1d(rwave, potential)
    
    potentialFull = interpolation_data(rFull) # interpolated potential
    print(potentialFull.shape[0])
    # ----- generate coefs for first and second derivative -----
    fdm = generate_coefs.coef(1, stencil_number)
    sdm = generate_coefs.coef(2, stencil_number)
    
    
    #print(rFull.shape)
    dim = rFull.shape[0]
    # -------------- generate the forward and backward number-------
    maxR = int(stencil_number/2)
    print(maxR)
    # -------------- assign to kinetic matrix ------------------
    k= np.zeros((dim,stencil_number),dtype=float)
    for j in range(maxR):
        for i in range(stencil_number):
        #forward
            k[j,i]=-(fdm[j,i]/(step*rFull[j]*reduced_mass)+sdm[j,i]/(step**2 *reduced_mass*2))
    
    for j in range(maxR,dim-maxR):
        for i in range(stencil_number):
            #central
             k[j,i]=-(fdm[maxR,i]/(step*rFull[j]*reduced_mass )+sdm[maxR,i]/(step**2 *reduced_mass*2) )
    
    for j in range(dim-maxR,dim):
        s=j-dim+stencil_number
        for i in range(stencil_number):
        #backward
            k[j,i]=-(fdm[s,i]/(step*rFull[j]*reduced_mass)+sdm[s,i]/(step**2 *reduced_mass*2))
            
            
    #--------------assign to V matrix ------------
    v = np.zeros((dim, dim), dtype=float)
    for i in range(0, rFull.shape[0]):
        v [i][i]=potentialFull[i]
    #--------------assign to H matrix ------------
    H = np.zeros((dim, dim), dtype=float)
    for i in range(0, rFull.shape[0]):
        H [i][i]=potentialFull[i]
    
    for j in range(maxR+1):
        for i in range(stencil_number):
            H[j,i]=H[j,i]+k[j,i]
    
    for j in range(maxR+1,dim-maxR-1):
        for i in range(stencil_number):
                H[j,i+j-maxR]=H[j,i+j-maxR]+k[j,i]
    
    for j in range(dim-maxR-1,dim):
        for i in range(stencil_number):
            H[j,i+dim-stencil_number]=H[j,i+dim-stencil_number]+k[j,i]
        # assign the diagonal term
    #    J_term = (J*(J+1)) / (2*mass*(rwave[i])**2 )
    #    H[i,i]=J_term + potential[i]

        # assign the derivatives


    return H,k,v, rFull

#-------------------------------------------------------------------


def compute_for_molecule (reducedMass,distance, potential, stencil, step):
    
    out = gen_H_matrix(reducedMass, distance , potential , stencil, step)
    
    H = out [0]
    k = out [1]
    rFull = out[3]
    
    
    w,v = np.linalg.eig(H)

    # sorting eigenvalues and corresponding eigenvectors
    idx = w.argsort()[::-1]   
    w = w[idx]
    v = v[:,idx]
    wReal = w.real
    
    dim=rFull.shape[0]
    
    # print fundamental
    conv_hartree_to_wavenumber = constants.value(u'hartree-inverse meter relationship') / 1e2
    fundamental_in_cm_inv = ((wReal [dim-2] - wReal [dim-1] ) * conv_hartree_to_wavenumber)
    print (wReal [dim-2] , wReal [dim-1])
    print(fundamental_in_cm_inv)
    
    return w,v, rFull
        
    