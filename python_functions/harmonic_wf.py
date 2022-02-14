# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from dep import diatomicSE
from dep import generate_H


Angs_to_Bohr_conversion = (1/0.52917721067)

######################################################
######################################################


def harmonic_potential(omega, redMass, re, distance):
    return (redMass/2)* (omega**2) * (distance-re)**2

######################################################
# Data from Gilson and coworkers (JRS 1980)

re_Angstrom = 1.09768
re_bohr = re_Angstrom * Angs_to_Bohr_conversion

omega = 2319.0679

print ('\t Data for \n\t\t for 14N 15N : r = ',re_Angstrom, 'A,  (',round(re_bohr,4),'a.u.)' , ' omega = ' , omega)
#######################################################
#######################################################


def gen_harmonic_potential_solve (zA, iMassA, zB, iMassB , omega, re):

    hartree_wavenumber_conversion = 219474.63137020
    omega_Hartree = omega / (hartree_wavenumber_conversion)

    # modify following for each molecule
    xfit = np.arange(0.5, 4.75 , 0.001)

    reducedMass = generate_H.reduced_mass(zA, iMassA, zB, iMassB)
    print('\t Reduced mass : ', reducedMass)

    harmonicPES = harmonic_potential( omega_Hartree , reducedMass, re, xfit)

    step = 0.004

    out = diatomicSE.get_eigenvalue_harmonic_J( zA, iMassA, zB, iMassB\
                                                        , xfit , harmonicPES, 5, step , 0 )

    e=out[0]
    wfn=out[1]
    r_wfn=out[2]

    print('\n\t omega/2 = ', omega/2,'\n')

    # print energies -----------------
    for i in range(8):

        if i==0:
            print( i,'\t', round(e[i][0],4 ) )
        if i>0:
            print(  i,'\t', round(e[i][0],4 ), '\t', round(e[i][0]-e[i-1][0], 4 ))
    #---------------------------------
    #return wavefunctions (harmonic)

    wfn_v0 = wfn[:,0]
    wfn_v1 = wfn[:,1]

    plt.plot(r_wfn, wfn_v0,'b', r_wfn, wfn_v1,'r')

    # export to txt
    np.savetxt('wfns_h/N2_h_v0J0.txt', wfn_v0, fmt='%3.15f')
    np.savetxt('wfns_h/N2_h_v1J0.txt', wfn_v1, fmt='%3.15f')
    np.savetxt('wfns_h/N2_h_r.txt', r_wfn, fmt='%3.15f')

    return wfn, r_wfn


######################################################
