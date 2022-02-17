
# Raman-Intensity-Approximn-Test

Set of python functions and datasets (rovibrational wavefunctions and polarizabilities) used to test the harmonic approximation for computing Raman intensities for the fundamental vibration in selected diatomic molecules.

---

**Table of Contents**

* [Raman-Intensity-Approximn-Test](#raman-intensity-approximn-test)
   * [About](#about)
   * [Bibliography](#bibliography)
   * [Authors](#authors)

---

## About

The double harmonic approximation is widely used in theoretical chemistry to model the frequencies and intensities of vibrational transitions in molecules. With datasets and tool in this repository, we test the harmonic approximation for computed Raman intensities.

**This repository includes :**
 - python module `diatomicSE`  corresponding to the program solving the 1D Schroeginger equation for the diatomic molecule.
 - python functions which use reference datasets on experimental ro-vibrational energies to determine the accurate wavefunctions (based on the `diatomicSE` module).
 - python function which computes harmonic vibrational wavefunctions for given values of molecular parameters.
 - examples to use the diatomicSE module (the mentioned cases)
 - distance-dependent static (wavelength independent) polarizabilities for selected diatomic molecules (H2, HF, HCl, CO, N<sub>2</sub> and F<sub>2</sub>) computed using different *ab initio* techniques.
 - ro-vibrational wavefunctions (anharmonic and those computed under harmonic approximation) for all the above molecules and a few of their isotopologues.



---


## Bibliography

- **NumPy** : Charles R. Harris, K. Jarrod Millman, Stéfan J. van der Walt, Ralf Gommers, Pauli Virtanen, David Cournapeau, Eric Wieser, Julian Taylor, Sebastian Berg, Nathaniel J. Smith, Robert Kern, Matti Picus, Stephan Hoyer, Marten H. van Kerkwijk, Matthew Brett, Allan Haldane, Jaime Fernández del Río, Mark Wiebe, Pearu Peterson, Pierre Gérard-Marchant, Kevin Sheppard, Tyler Reddy, Warren Weckesser, Hameer Abbasi, Christoph Gohlke & Travis E. Oliphant. "Array programming with NumPy", *Nature*, 585, 357–362 (2020) [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)

- **Scipy** : Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) "SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python". *Nature Methods*, 17(3), 261-272. [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)

- **Non-linear optimization in SciPy** :  Travis E. Oliphant. "Python for Scientific Computing, *Computing in Science & Engineering*, 9, 10-20 (2007), DOI:10.1109/MCSE.2007.58


- **Matplotlib**  : J. D. Hunter, "Matplotlib: A 2D Graphics Environment", *Computing in Science & Engineering*, vol. 9, no. 3, pp. 90-95, 2007.


- **Orthogonal Distance Regression as used in SciPy** : (***i***) P. T. Boggs, R. Byrd, R. Schnabel, SIAM J. Sci. Comput. 1987, 8, 1052. (***ii***) P. T. Boggs, J. R. Donaldson, R. h. Byrd, R. B. Schnabel, ACM Trans. Math. Softw. 1989, 15, 348. (***iii***) J. W. Zwolak, P. T. Boggs, L. T. Watson, ACM Trans. Math. Softw. 2007, 33, 27. (***iv***)  P. T. Boggs and J. E. Rogers, “Orthogonal Distance Regression,” in “Statistical analysis of measurement error models and applications: proceedings of the AMS-IMS-SIAM joint summer research conference held June 10-16, 1989,” Contemporary Mathematics, vol. 112, pg. 186, 1990.

## Authors
Yen Bang Chao and Ankit Raj
