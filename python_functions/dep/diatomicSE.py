# -*- coding: utf-8 -*-

import numpy as np
from scipy import interpolate
from scipy import constants


from dep import generate_coefs
from dep import generate_H

import scipy.optimize as opt

import matplotlib.pyplot as plt


######################################################

#-------------------------------------------------------------------

def symmetrize_Hmat(inputH, stencil_number):

    maxR = int(stencil_number/2)
    nRows=inputH.shape[0]
    nCols=inputH.shape[1]

    for i in range(maxR):
        inputH=np.delete(inputH, (0), axis=0)
        inputH=np.delete(inputH, (0), axis=1)

        #print( 'test >>', inputH.shape,nRows-2*i-2)

        inputH=np.delete(inputH, (nRows-2*i-2), axis=0)
        inputH=np.delete(inputH, (nRows-2*i-2), axis=1)

        #print('after a loop >> ', inputH.shape)

    # print (' inputH is redimensioned. New shape :', inputH.shape)
    return inputH

#-------------------------------------------------------------------



##############################################################################


#--------------------------------------------------------------------
def get_eigenvalue_J(zA, iMassA, zB, iMassB, rwave , potential, stencil_number, step, J):
    '''


    Parameters
    ----------
    zA : TYPE
        DESCRIPTION.
    iMassA : TYPE
        DESCRIPTION.
    zB : TYPE
        DESCRIPTION.
    iMassB : TYPE
        DESCRIPTION.
    rwave : TYPE
        DESCRIPTION.
    potential : TYPE
        DESCRIPTION.
    stencil_number : TYPE
        DESCRIPTION.
    step : TYPE
        DESCRIPTION.

    Returns
    -------
    rFull : TYPE
        DESCRIPTION.
    w : TYPE
        DESCRIPTION.
    v : TYPE
        DESCRIPTION.

    '''
    #start = time.time()
    #----- calculate reduce mass -----------
    reducedMass = generate_H.reduced_mass(zA, iMassA, zB, iMassB)
    #---- generate full r vector ----------
    rMin = np.amin(rwave)
    rMax = np.amax(rwave)

    rFull = np.arange(rMin, rMax, step, dtype=float)

    # ---- interpolate potential -----------
    interpolation_data = interpolate.interp1d(rwave, potential)

    potentialFull = interpolation_data(rFull) # interpolated potential

    #print(potentialFull.shape[0])
    # ----- generate coefs for first and second derivative -----
    fdm = generate_coefs.coef(1, stencil_number)
    sdm = generate_coefs.coef(2, stencil_number)


    #print(rFull.shape)
    dim = rFull.shape[0]
    # -------------- generate the forward and backward number-------
    maxR = int(stencil_number/2)
    # -------------- assign to kinetic matrix ------------------
    k= np.zeros((dim,stencil_number),dtype=float)
    for j in range(maxR):
        for i in range(stencil_number):
        #forward
            k[j,i]=-(fdm[j,i]/(step*rFull[j]*reducedMass)+sdm[j,i]/(step**2 *reducedMass*2))

    for j in range(maxR,dim-maxR):
        for i in range(stencil_number):
            #central
             k[j,i]=-(fdm[maxR,i]/(step*rFull[j]*reducedMass )+sdm[maxR,i]/(step**2 *reducedMass*2) )

    for j in range(dim-maxR,dim):
        s=j-dim+stencil_number
        for i in range(stencil_number):
        #backward
            k[j,i]=-(fdm[s,i]/(step*rFull[j]*reducedMass)+sdm[s,i]/(step**2 *reducedMass*2))

    #--------------assign to H matrix ----------------------
    H = np.zeros((dim, dim), dtype=float)
    for i in range(0, rFull.shape[0]):
        H [i][i]=potentialFull[i] + (J*(J+1))/(2*reducedMass*rFull[i]**2)

    for j in range(maxR+1):
        for i in range(stencil_number):
            H[j,i]=H[j,i]+k[j,i]

    for j in range(maxR+1,dim-maxR-1):
        for i in range(stencil_number):
                H[j,i+j-maxR]=H[j,i+j-maxR]+k[j,i]

    for j in range(dim-maxR-1,dim):
        for i in range(stencil_number):
            H[j,i+dim-stencil_number]=H[j,i+dim-stencil_number]+k[j,i]

    #--------------eigenvalue in wavenumber------------------------------

    # symmetrize H matrix
    H=symmetrize_Hmat(H, stencil_number)
    Hshape=H.shape[0]
    w,v = np.linalg.eig(H)
    sort_index = np.argsort(w)

    # sorting of eigenvectors is needed
    sorted_energy = np.array(w)[sort_index]
    sorted_wfn = np.zeros((w.shape[0], w.shape[0]))
    for i in range(w.shape[0]):
        sorted_wfn [:, i] = v.real[:, sort_index[i] ]

    # print fundamental
    conv_hartree_to_wavenumber = constants.value(u'hartree-inverse meter relationship') / 1e2
    wcm_inv = np.zeros((Hshape,1),dtype=float)
    for i in range(Hshape):
        wcm_inv[i] = sorted_energy[i] *conv_hartree_to_wavenumber


    #end = time.time()
    #print(end - start)
    #time=6.40555471
    return  wcm_inv


##############################################################################
def get_fullsolution_J(zA, iMassA, zB, iMassB, rwave , potential, stencil_number, step, J):

    #start = time.time()
    #----- calculate reduce mass -----------
    reducedMass = generate_H.reduced_mass(zA, iMassA, zB, iMassB)
    #---- generate full r vector ----------
    rMin = np.amin(rwave)
    rMax = np.amax(rwave)
    
    rFull = np.arange(rMin, rMax, step, dtype=float)
    
    # ---- interpolate potential -----------
    interpolation_data = interpolate.interp1d(rwave, potential)
    
    potentialFull = interpolation_data(rFull) # interpolated potential
    
    #print(potentialFull.shape[0])
    # ----- generate coefs for first and second derivative -----
    fdm = generate_coefs.coef(1, stencil_number)
    sdm = generate_coefs.coef(2, stencil_number)
    
    
    #print(rFull.shape)
    dim = rFull.shape[0]
    # -------------- generate the forward and backward number-------
    maxR = int(stencil_number/2)
    # -------------- assign to kinetic matrix ------------------
    k= np.zeros((dim,stencil_number),dtype=float)
    for j in range(maxR):
        for i in range(stencil_number):
        #forward
            k[j,i]=-(fdm[j,i]/(step*rFull[j]*reducedMass)+sdm[j,i]/(step**2 *reducedMass*2))
    
    for j in range(maxR,dim-maxR):
        for i in range(stencil_number):
            #central
             k[j,i]=-(fdm[maxR,i]/(step*rFull[j]*reducedMass )+sdm[maxR,i]/(step**2 *reducedMass*2) )
    
    for j in range(dim-maxR,dim):
        s=j-dim+stencil_number
        for i in range(stencil_number):
        #backward
            k[j,i]=-(fdm[s,i]/(step*rFull[j]*reducedMass)+sdm[s,i]/(step**2 *reducedMass*2))
            
    #--------------assign to H matrix ----------------------
    H = np.zeros((dim, dim), dtype=float)
    for i in range(0, rFull.shape[0]):
        H [i][i]=potentialFull[i] + (J*(J+1))/(2*reducedMass*rFull[i]**2)
    
    for j in range(maxR+1):
        for i in range(stencil_number):
            H[j,i]=H[j,i]+k[j,i]
    
    for j in range(maxR+1,dim-maxR-1):
        for i in range(stencil_number):
                H[j,i+j-maxR]=H[j,i+j-maxR]+k[j,i]
    
    for j in range(dim-maxR-1,dim):
        for i in range(stencil_number):
            H[j,i+dim-stencil_number]=H[j,i+dim-stencil_number]+k[j,i]
    
    #--------------eigenvalue in wavenumber------------------------------
    
    # symmetrize H matrix
    H=symmetrize_Hmat(H, stencil_number)
    Hshape=H.shape[0]
    w,v = np.linalg.eig(H)
    sort_index = np.argsort(w)
    
    # sorting of eigenvectors is needed 
    sorted_energy = np.array(w)[sort_index]
    sorted_wfn = np.zeros((w.shape[0], w.shape[0]))
    for i in range(w.shape[0]):
        sorted_wfn [:, i] = v.real[:, sort_index[i] ]

    # print fundamental
    conv_hartree_to_wavenumber = constants.value(u'hartree-inverse meter relationship') / 1e2
    wcm_inv = np.zeros((Hshape,1),dtype=float)
    for i in range(Hshape):
        wcm_inv[i] = sorted_energy[i] *conv_hartree_to_wavenumber

    rfull=np.delete(rFull,[0,1,-1,-2],axis=0)
    #end = time.time()
    #print(end - start)
    #time=6.40555471
    return  rfull, wcm_inv, sorted_wfn




#--------------------------------------------------------------------
def get_eigenvalue_harmonic_J(zA, iMassA, zB, iMassB, rwave , potential, stencil_number, step, J):
    '''


    Parameters
    ----------
    zA : TYPE
        DESCRIPTION.
    iMassA : TYPE
        DESCRIPTION.
    zB : TYPE
        DESCRIPTION.
    iMassB : TYPE
        DESCRIPTION.
    rwave : TYPE
        DESCRIPTION.
    potential : TYPE
        DESCRIPTION.
    stencil_number : TYPE
        DESCRIPTION.
    step : TYPE
        DESCRIPTION.

    Returns
    -------
    rFull : TYPE
        DESCRIPTION.
    w : TYPE
        DESCRIPTION.
    v : TYPE
        DESCRIPTION.

    '''
    #start = time.time()
    #----- calculate reduce mass -----------
    reducedMass = generate_H.reduced_mass(zA, iMassA, zB, iMassB)
    #---- generate full r vector ----------
    rMin = np.amin(rwave)
    rMax = np.amax(rwave)

    rFull = np.arange(rMin, rMax, step, dtype=float)

    # ---- interpolate potential -----------
    interpolation_data = interpolate.interp1d(rwave, potential)

    potentialFull = interpolation_data(rFull) # interpolated potential

    #print(potentialFull.shape[0])
    # ----- generate coefs for first and second derivative -----
    fdm = generate_coefs.coef(1, stencil_number)
    sdm = generate_coefs.coef(2, stencil_number)


    #print(rFull.shape)
    dim = rFull.shape[0]
    # -------------- generate the forward and backward number-------
    maxR = int(stencil_number/2)
    # -------------- assign to kinetic matrix ------------------
    k= np.zeros((dim,stencil_number),dtype=float)
    for j in range(maxR):
        for i in range(stencil_number):
        #forward
            k[j,i]=-(fdm[j,i]/(step*rFull[j]*reducedMass)+sdm[j,i]/(step**2 *reducedMass*2))

    for j in range(maxR,dim-maxR):
        for i in range(stencil_number):
            #central
             k[j,i]=-(fdm[maxR,i]/(step*rFull[j]*reducedMass )+sdm[maxR,i]/(step**2 *reducedMass*2) )

    for j in range(dim-maxR,dim):
        s=j-dim+stencil_number
        for i in range(stencil_number):
        #backward
            k[j,i]=-(fdm[s,i]/(step*rFull[j]*reducedMass)+sdm[s,i]/(step**2 *reducedMass*2))

    #--------------assign to H matrix ----------------------
    H = np.zeros((dim, dim), dtype=float)
    for i in range(0, rFull.shape[0]):
        H [i][i]=potentialFull[i] #+ (J*(J+1))/(2*reducedMass*rFull[i]**2)

    for j in range(maxR+1):
        for i in range(stencil_number):
            H[j,i]=H[j,i]+k[j,i]

    for j in range(maxR+1,dim-maxR-1):
        for i in range(stencil_number):
                H[j,i+j-maxR]=H[j,i+j-maxR]+k[j,i]

    for j in range(dim-maxR-1,dim):
        for i in range(stencil_number):
            H[j,i+dim-stencil_number]=H[j,i+dim-stencil_number]+k[j,i]

    #--------------eigenvalue in wavenumber------------------------------

    # symmetrize H matrix
    H=symmetrize_Hmat(H, stencil_number)
    Hshape=H.shape[0]
    w,v = np.linalg.eig(H)
    sort_index = np.argsort(w)

    # sorting of eigenvectors is needed
    sorted_energy = np.array(w)[sort_index]
    sorted_wfn = np.zeros((w.shape[0], w.shape[0]))
    for i in range(w.shape[0]):
        sorted_wfn [:, i] = v.real[:, sort_index[i] ]

    # print fundamental
    conv_hartree_to_wavenumber = constants.value(u'hartree-inverse meter relationship') / 1e2
    wcm_inv = np.zeros((Hshape,1),dtype=float)
    for i in range(Hshape):
        wcm_inv[i] = sorted_energy[i] *conv_hartree_to_wavenumber


    # trimming the rFull after symmetrization procedure
    rFull=np.delete(rFull,[0,1,-1,-2],axis=0)

    #end = time.time()
    #print(end - start)
    #time=6.40555471
    return  wcm_inv, sorted_wfn, rFull


##############################################################################
