# README
This file explains how to replicate the experiments. It is assumed that the system has a working installation of [NLEVP](https://github.com/ftisseur/nlevp).

## Chapter 2
### Examples 2.2/2.3:
  Run RIMrandom.m and RIMbutt.m.
## Chapter 3
### Experiment in Section 3.4.1
For the butterfly problem, run thesisBUTT.m For the hadeler problem, run quadrature_hadeler.m. Please, **make sure to be inside the nlevpChap folder** and to use the the contourSolver.m therein, because it is the only version that approximates the integral for different number of quadrature points.
### Example 3.1
Run poleInfluence.m.
### Figure 3.4
Run filterSurf.m.
### Example 3.2
Run meromorphic2.m.
### Example 3.3
Run recursiveRef.m.
### Example 3.4
Run pdde_symProbing.m.
## Chapter 4
The file do_nlevp_test.m performs all the tests in section 4.4. One can change the accuracy by modifying line 6. In order to test a specific subset of problem, just modify the cell "nep" (see lines 21-31). The file samplePointsPolesPlot.m produces Figure 4.3. The experiments in Section 4.5 are performed by comparenep2rat.m, while the tables are built by the files make_*_table.m.
## Chapter 5
### Example 5.1 and Example 5.5
Run harmonics.m.
### Figure 5.2
Run isolRoots.m.
### Figure 5.3
Run toy1.m.
### Examples 5.7 and 5.8
Run example1.m and example2New.m.
### Example 5.10
Run exponentialPoly.m.
