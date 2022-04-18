# Details of the provided data

The following directory contains calculated data for the studied molecules. This data includes : (**i**) distance dependent polarizability computed using four methods (Hartree-Fock, DFT, CCSD(T) and CASSCF), and (**ii**) ro-vibrational wavefunctions (exact and harmonic).

-----


**Table of Contents**


[TOC]

-----


## Polarizability data

Polarizability data is placed in the `polarizability` directory. Format of the data is as follows (example from CO_HartreeFock.dat).

| Distance (bohr) | alpha_zz      | alpha_xx      |
|-----------------|---------------|---------------|
| 0.450000000     | 156.143890730 | 100.897769904 |
| 0.550000000     | 185.896547717 | 124.368211396 |
| 0.650000000     | 203.049235366 | 141.128976456 |
| 0.750000000     | 69.656082711  | 73.189204053  |

Polarizability values are given in atomic units.




## Ro-vibrational wavefunctions

Ro-vibrational wavefunctions and the potential energy curve for a molecule are placed in the `wfns` directory. Each data file is a one-dimensional column of data.

Wavefunctions are labelled using the vibrational and rotational quantum number (both for the exact and harmonic type). Distance data (in bohr) is common to both wavefunctions and the PES.

## Potential energy surface

PES corresponds to the anharmonic potential energy curve obtained using analysis scheme developed in this work, where we utilized experimental datasets on transition energies to iteratively obtain the PES via non-linear optimization scheme.

## Studied molecules

H<sub>2</sub>, HD, D<sub>2</sub>, HF,  HCl, DCl, CO, <sup>13</sup>C<sup>16</sup>O   N<sub>2</sub>, <sup>14</sup>N<sup>15</sup>N,  and F<sub>2</sub>
