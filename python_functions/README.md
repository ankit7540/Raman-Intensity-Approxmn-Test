# diatomicSE

**An implementation of collocation method for solution of the 1D Schroedinger equation for diatomic molecules**

Python module based on `numpy` and `scipy`, for solution of the 1D radial nuclear equation for a given potential energy curve.

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
zB             : atomic number (atom 2)
iMassB         : atomic mass (atom 2)
rwave          : 1D vector of internuclear distance (for potential energy curve)
potential      : 1D vector of potential energy curve
stencil_number : order to numerical derivative used to model partial derivatives
step           : step size of the distance vector used to construct H-matrix
J              : required rotational state (solutions for all vibrational states are found at once)

```


and

```
get_fullsolution_J(zA, iMassA, zB, iMassB, rwave , potential, stencil_number, step, J)
# with same argument definitions as above
```


- Other files :
 -- `transition.py`    (in dep directory)
    Contains functions for computing transition frequencies from 1D vector of energies of ro-vibrational states. See [getPES.py](https://github.com/ankit7540/Raman-Intensity-Approxmn-Test/blob/main/python_functions/getPES.py) file and jupyter example for usage.

  -- `PES_functions.py` (in dep directory)
    Contains implementation of the analytical functions for generating potential energy curves of diatomic molecules. See examples below on usage.

  -- `harmonic_wf.py`
     Contains function to generate the harmonic wavefunctions using user provided molecular parameters. See "Constant of Diatomic Molecules" by Herzberg G. and Huber K.P. for tabulation of molecular params.

- Directory structure
```
.
├── dep
│   ├── diatomicSE.py
│   ├── generate_coefs.py
│   ├── generate_H.py
│   ├── __init__.py
│   ├── PES_functions.py
│   ├── quadrature.py
│   └── transition.py
├── examples
│   ├── general_solution.ipynb
│   ├── N2_ref_data
│   │   └── N2_ref_J4_v3.txt
│   └── Obtaining PES .ipynb
├── getPES.py
├── harmonic_wf.py
└── README.md
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
In the above example, `rwave` is a one-dimensional numpy array representing the internuclear distance. `PES_N14N15` is a one-dimensional numpy array representing the potential energy surface of the <sup>14</sup>N<sup>15</sup>N molecule.

### Standard installation
Clone the repository and use.  Make sure to install the required modules (`numpy`, `scipy` and `periodictable`).

```
git clone https://github.com/ankit7540/Raman-Intensity-Approxmn-Test
```

-----

### Examples

Following shows an example for the CO molecule.

```
# generating the PES curve for CO
distance_vector=np.arange(0.75, 5.25, 0.025)
param =np.array([0.37563416, 1.22186557, 2.13081292])  # PES parameters for CO
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

##### General solution for the CO molecule [(see file here)](https://github.com/ankit7540/Raman-Intensity-Approxmn-Test/blob/main/python_functions/examples/general_solution.ipynb)

##### Determination of accurate PES using reference experimental data as reference  for the N<sub>2</sub> molecule [(see file here)](https://github.com/ankit7540/Raman-Intensity-Approxmn-Test/blob/main/python_functions/examples/Obtaining%20PES%20.ipynb)


### Authors

Yen Bang Chao and Ankit Raj

### Bibliography

 - U. M. Ascher and L. R. Petzold, Computer methods for ordinary differential equations and differential-algebraic equations (Society for Industrial and Applied Mathematics, 2009).

 - E. Hairer, S. P. Nørsett, and G. Wanner, Solving ordinary differential equations i: nonstiff problems (Springer, 1993).

- Polarizability tensor invariants of H<sub>2</sub>, HD, and D<sub>2</sub>, J. Chem. Phys. 148, 104308 (2018), DOI: [10.1063/1.5011433](https://doi.org/10.1063/1.5011433)

- Vibration–rotation interactions in H<sub>2</sub>, HD and D<sub>2</sub> : centrifugal distortion factors and the derivatives of polarisability invariants, Mol. Phys 2019, DOI:[10.1080/00268976.2019.1632950](https://doi.org/10.1080/00268976.2019.1632950)


---
