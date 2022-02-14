# diatomicSE

**An implementation of collocation method for solution of the 1D Schroedinger equation for diatomic molecules**

Python module based on numpy and scipy, for the solution of 1D radial nuclear equation, for diatomic molecules for a given potential energy curve.

**Table of Contents**

[TOC]

---

### About
The solution of the 1D radial nuclear equation gives the energies of the ro-vibrational
states and the corresponding wavefunctions, for diatomic molecules.

![1D radial nuclear equation](img/equation.svg)

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


### Tested with

+ `python 3.5+`
    + `numpy 1.15.2+`
	+ `scipy 1.1.0+`


-----


### Usage

```

import diatomicSE

# for 14N15N molecule

diatomicSE.get_eigenvalue_J( 7, 14, 7, 15, rwave , PES_N14N15, 5, 0.005 , 0 )

```

### Standard installation

```
git clone https://github.com/ankit7540/Raman-Intensity-Approxmn-Test

```


-----

### Examples



```



```

#### Jupyter notebook

##### General solution
##### Determination of accurate PES using reference experimental data as reference


### Authors

Yen Bang Chao and Ankit Raj

### Bibliography





---
