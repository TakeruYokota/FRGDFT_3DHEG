# Functional-renormalization-group-aided density functional theory (FRG-DFT) for 3D homogeneous electron gas

C++ code for FRG-DFT calculation of the correlation energy of 3D homogeneous electron gas. This returns the results by FRG-DFT and Gell-Mann-Brueckner resummation
for given ranges of the Wigner-Seitz radius and the spin polarization.

For details of the formalism, see [Phys. Rev. Research 3, L012015 (2021)](https://doi.org/10.1103/PhysRevResearch.3.L012015) and [Phys. Rev. B 105, 035105 (2022)](https://doi.org/10.1103/PhysRevB.105.035105).

## Required libraries
GNU Scientific Library (GSL) and an MPI library are required.

## Compilation using cmake
1. Create a directory for compilation (e.g., `build`) and move into it.
2. After running `cmake ..`, the executable `build/prog` is generated.

## Input file for running
To run the calculation, an input file describing the calculation setting is needed. The file `run/input.txt` is an example. The contents are as follows:
```
ls      0
le      0
Nl      0
zs      0
ze      0
Nz      0
out_err 1
```
The first six variables set the range of the Wigner-Seitz radius $r_{\mathrm{s}}$ and the spin polarization $\zeta$. A log mesh is assigned for $r_{\mathrm{s}}$. In other words, an equally spaced mesh for $ l=\log_{10}r_{\mathrm{s}} $ is used. The range is given by `[ls, le]`, and this is divided into `Nl` intervals (thus, the number of the points is `Nl+1`). If `Nl=0`, only the $l=$`ls` case is calculated. For $\zeta$, an equally spaced mesh is assigned. Similarly to $l$, the range is given by `[zs, ze]` and this is divided into `Nz` intervals. If `Nz=0`, only the $\zeta=$`zs` case is calculated. 

GSL functions sometimes give warnings about the convergence of some integral. These warnings are visible if `out_err` is set to 1.

## How to run
As written in run/run.sh, one can run `build/prog` with the following code:
```
NUM_PROCS=1
mpirun -np $NUM_PROCS build/prog input.txt
```
For `NUM_PROCS`$>1$, the calculation is parallelized with respect to the set of $(l,\zeta)$.

## Output file
The results are output in `results.csv`. From the left column, the data mean the following values:
- $r_{\mathrm{s}}$
- $\zeta$
- FRG-DFT result of the correlation energy per particle $E_{\mathrm{corr}}/N$
- Gell-Mann-Bruckner result of $E_{\mathrm{corr}}/N$, 
- FRG-DFT result of the total energy per particle $E_{\mathrm{tot}}/N$
- Gell-Mann-Bruckner result of $E_{\mathrm{tot}}/N$,

## Quickstart
After compilation, one can quickly check if the calculation properly runs through the following procedure:
1. Run `testrun.sh` 
1. Move to `_testrun`
1. Run `run.sh`

## How to cite
Please cite our papers ([Phys. Rev. Research 3, L012015 (2021)](https://doi.org/10.1103/PhysRevResearch.3.L012015) and [Phys. Rev. B 105, 035105 (2022)](https://doi.org/10.1103/PhysRevB.105.035105)) if you use this code.
