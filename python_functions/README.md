# diatomicSE

**An implementation of collocation method for solution of the 1D Schroedinger equation for diatomic molecules**

Python module based on numpy and scipy, for the solution of 1D radial nuclear equation, for diatomic molecules for a given potential energy curve.

**Table of Contents**


* [diatomicSE](#diatomicse)
     * [About](#about)
         * [Primary function](#primary-function)
     * [Requirements](#requirements)
     * [Tested with](#tested-with)
     * [Usage](#usage)
     * [Standard installation](#standard-installation)
     * [Examples](#examples)
          * [General solution](#general-solution)
          * [Determination of accurate PES using reference experimental data as reference](#determination-of-accurate-pes-using-reference-experimental-data-as-reference)
     * [Authors](#authors)
     * [Bibliography](#bibliography)


---

### About
The solution of the 1D radial nuclear equation gives the energies of the ro-vibrational
states and the corresponding wavefunctions, for diatomic molecules.

![1D radial nuclear equation](https://github.com/ankit7540/Raman-Intensity-Approxmn-Test/blob/main/img/equation.svg)

In the present approach we use the collocation method using numerical derivatives to map the above equation on a matrix. Diagonalization of the matrix gives us the eigenvalues (energies) and the eigenvectors (wavefunctions).


#### Primary function

```
get_eigenvalue_J(zA, iMassA, zB, iMassB, rwave , potential, stencil_number, step, J)

za             : atomic number (atom 1)
iMassA         : atomic mass (atom 1)
zB             : atomic mass (atom 2)
iMassB         : atomic mass (atom 2)
rwave          : 1D vector of internuclear distance (for potential energy curve)
potential      : 1D vector of potential energy curve
stencil_number : order to numerical derivative used to model partial derivatives
step           : step size of the distance vector used to contruct H-matrix
J              : required rotational state (solutions for all vibrational states are found at once)

```


### Requirements

+ `python3`
    + `numpy`
	+ `scipy`
    + `periodictable`


### Tested with

+ `python 3.5+`
    + `numpy 1.15.2+`
	+ `scipy 1.1.0+`
    + `periodictable 1.6.0`


-----


### Usage

Import the module and call the functions defined within the module.

```
import diatomicSE

# for 14N15N molecule
diatomicSE.get_eigenvalue_J( 7, 14, 7, 15, rwave , PES_N14N15, 5, 0.005 , 0 )

```
In the above example, `rwave` is a one-dimensional numpy array representing the internuclear distance. `PES_N14N15` is a one-dimensional numpy array representing the potential energy surface of the molecule, here it is <sup>14</sup>N<sup>15</sup>N.

### Standard installation
Simply clone the repository and use.  Make sure to install the required modules (`numpy`, `scipy` and `periodictable`).

```
git clone https://github.com/ankit7540/Raman-Intensity-Approxmn-Test

```

-----

### Examples

```
# generating the PES curve for CO
distance_vector=np.arange(0.75, 5.25, 0.025)
PES = PES_functions.rydberg_potential_3param(param ,  distance_vector )

step = 0.004 # step (h) in the numerical derivative
J_level = 0  # rotational state

# ---------------------------------------------------------
# calculation for the energies
# Molecule : CO

e = diatomicSE.get_eigenvalue_J( 6, 12, 8, 16, distance_vector, PES, 5, step , J_level  )

# e is a vector of energies (in wavenumbers)
# ---------------------------------------------------------

# ---------------------------------------------------------
# calculation for energies and wavefunctions
# Molecule : CO

solution = diatomicSE.get_fullsolution_J( 6, 12, 8, 16, distance_vector, PES, 5, step , J_level  )

# solution is a tuple.
# solution [0] : distance
# solution [1] : energies (Hartree)
# solution [2] : wavefunctions (across columns), correspond to vibrational levels
# ---------------------------------------------------------

```

**Examples are shown in jupyter notebooks:**

##### General solution for CO molecules [(see file here)](https://github.com/ankit7540/Raman-Intensity-Approxmn-Test/blob/main/python_functions/examples/general_solution.ipynb)

##### Determination of accurate PES using reference experimental data as reference  for N2 molecule [(see file here)](https://github.com/ankit7540/Raman-Intensity-Approxmn-Test/blob/main/python_functions/examples/Obtaining%20PES%20.ipynb)


### Authors

Yen Bang Chao and Ankit Raj

### Bibliography

 - U. M. Ascher and L. R. Petzold, Computer methods for ordinary differential equa-
tions and differential-algebraic equations (Society for Industrial and Applied
Mathematics, 2009).

 - E. Hairer, S. P. Nørsett, and G. Wanner, Solving ordinary differential equations
i: nonstiff problems (Springer, 1993).

- Polarizability tensor invariants of H2, HD, and D2, J. Chem. Phys. 148, 104308 (2018), DOI: [10.1063/1.5011433](https://doi.org/10.1063/1.5011433)

- Vibration–rotation interactions in H2, HD and D2 : centrifugal distortion factors and the derivatives of polarisability invariants, Mol. Phys 2019, DOI:[10.1080/00268976.2019.1632950](https://doi.org/10.1080/00268976.2019.1632950)


---
