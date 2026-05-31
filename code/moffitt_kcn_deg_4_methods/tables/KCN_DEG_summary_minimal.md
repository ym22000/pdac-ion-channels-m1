# KCN DEG Summary (Minimal)

This file summarizes the four main DEG frameworks used to prioritize KCN genes.

Minimal columns only:
- `log2FC` when available
- `beta` for the hurdle mixed model
- `FDR`: BH-adjusted significance
- `Lecture`: biological direction (`Activated up`, `Normal up`, `qPSC up`, `myCAF down`, etc.)

Important note:
- `Test 1` and `Test 2`: `log2FC`
- `Test 3`: `beta` from the hurdle mixed model
- `Test 4`: `avg_log2FC` from MAST

## Test 1 - DESeq2 Activated vs Normal

| KCN | log2FC | FDR | Lecture |
|---|---:|---:|---|
| KCNAB1 | -3.029 | 8.53e-20 | Normal up |
| KCNA5 | -3.541 | 3.00e-11 | Normal up |
| KCNMA1 | -2.507 | 5.44e-11 | Normal up |
| KCNK17 | -1.823 | 6.31e-09 | Normal up |
| KCNK6 | 1.208 | 2.25e-06 | Activated up |
| KCNC4 | -1.427 | 1.49e-04 | Normal up |
| KCNQ4 | -1.831 | 7.84e-04 | Normal up |
| KCNT2 | 1.806 | 8.88e-04 | Activated up |
| KCND3 | 1.603 | 1.20e-03 | Activated up |
| KCNN3 | 2.239 | 2.18e-03 | Activated up |
| KCNJ8 | 0.571 | 4.52e-03 | Activated up |
| KCNS3 | -0.462 | 1.40e-02 | Normal up |
| KCNE4 | -0.401 | 3.07e-02 | Normal up |
| KCNB1 | -1.891 | 4.51e-02 | Normal up |

## Test 2 - DESeq2 Subtype vs Rest

| KCN | log2FC | FDR | Lecture |
|---|---:|---:|---|
| KCNJ8 | 1.670 | 9.11e-07 | qPSC up |
| KCNA5 | -2.638 | 3.79e-06 | qPSC down |
| KCNAB1 | -1.894 | 4.27e-05 | qPSC down |
| KCNMB1 | -2.137 | 7.37e-05 | qPSC down |
| KCNK17 | 1.162 | 3.62e-04 | qPSC up |
| KCNE4 | 0.807 | 9.45e-03 | qPSC up |
| KCNK6 | 0.985 | 1.85e-02 | qPSC up |
| KCNK3 | -1.174 | 2.60e-02 | qPSC down |
| KCNAB1 | 3.468 | 1.58e-20 | smPSC up |
| KCNA5 | 3.209 | 2.81e-16 | smPSC up |
| KCNJ8 | -2.080 | 2.48e-14 | smPSC down |
| KCNMA1 | 3.728 | 1.30e-10 | smPSC up |
| KCNK6 | -1.819 | 1.92e-06 | smPSC down |
| KCNMB1 | 1.193 | 3.03e-05 | smPSC up |
| KCND3 | -2.434 | 2.92e-03 | smPSC down |
| KCNC4 | 0.989 | 4.63e-03 | smPSC up |
| KCNQ4 | 1.762 | 2.91e-02 | smPSC up |
| KCNJ8 | -3.164 | 4.04e-11 | myCAF down |
| KCNK17 | -5.620 | 1.50e-07 | myCAF down |
| KCNMB1 | -2.298 | 3.10e-06 | myCAF down |
| KCNK3 | -4.005 | 3.40e-04 | myCAF down |
| KCNA5 | -4.810 | 1.34e-03 | myCAF down |
| KCNMA1 | 2.002 | 3.63e-03 | myCAF up |
| KCND3 | -1.168 | 3.10e-02 | myCAF down |
| KCNK17 | -4.988 | 3.23e-03 | csCAF down |

## Test 3 - Hurdle Mixed Model

| KCN | beta | FDR | Lecture |
|---|---:|---:|---|
| KCND3 | 0.034 | 2.95e-30 | csCAF up |
| KCNMB1 | 0.102 | 1.82e-18 | csCAF up |
| KCNK15 | -0.003 | 1.86e-09 | csCAF down |
| KCNN3 | -0.025 | 2.85e-05 | csCAF down |
| KCNT2 | -0.027 | 7.30e-05 | csCAF down |
| KCND2 | -0.038 | 4.33e-04 | csCAF down |
| KCNK1 | 0.090 | 1.55e-05 | iCAF up |
| KCNC4 | 0.697 | 3.06e-02 | iCAF up |
| KCNMB1 | 0.062 | 9.02e-20 | myCAF up |
| KCND2 | 0.034 | 2.26e-18 | myCAF up |
| KCNK1 | -0.001 | 7.27e-10 | myCAF down |
| KCNMA1 | 0.017 | 2.83e-09 | myCAF up |
| KCNMB4 | 0.017 | 2.75e-06 | myCAF up |
| KCNN3 | 0.056 | 3.59e-05 | myCAF up |
| KCNN4 | 0.048 | 4.37e-05 | myCAF up |
| KCNJ15 | 0.037 | 3.20e-03 | myCAF up |
| KCNJ8 | 0.267 | 0 | qPSC up |
| KCNK17 | 0.084 | 2.41e-162 | qPSC up |
| KCNE4 | 0.073 | 7.60e-103 | qPSC up |
| KCNS3 | 0.041 | 1.01e-22 | qPSC up |
| KCNK6 | 0.014 | 3.03e-21 | qPSC up |
| KCNT2 | 0.053 | 2.32e-03 | qPSC up |
| KCNC4 | -0.013 | 3.93e-02 | qPSC down |
| KCNA7 | -0.012 | 8.72e-02 | qPSC down |
| KCNA5 | 0.036 | 0 | smPSC up |
| KCNAB1 | 0.035 | 2.40e-203 | smPSC up |
| KCNMB1 | -0.063 | 6.37e-113 | smPSC down |
| KCNK3 | 0.001 | 4.80e-41 | smPSC up |
| KCNMA1 | 0.015 | 1.34e-34 | smPSC up |
| KCNK17 | -0.067 | 1.64e-22 | smPSC down |
| KCNC4 | 0.007 | 2.94e-21 | smPSC up |
| KCNS3 | -0.009 | 6.59e-09 | smPSC down |

## Test 4 - MAST

| KCN | log2FC | FDR | Lecture |
|---|---:|---:|---|
| KCNJ8 | 2.288 | 0 | qPSC up |
| KCNK17 | 1.738 | 1.05e-197 | qPSC up |
| KCNE4 | 0.975 | 5.76e-120 | qPSC up |
| KCNMB1 | -1.545 | 2.43e-76 | qPSC down |
| KCNMA1 | -3.133 | 3.90e-46 | qPSC down |
| KCNK6 | 0.773 | 1.04e-24 | qPSC up |
| KCNS3 | 0.629 | 9.71e-24 | qPSC up |
| KCNA5 | -0.950 | 1.16e-21 | qPSC down |
| KCND2 | -1.605 | 1.43e-18 | qPSC down |
| KCNAB1 | -1.107 | 4.59e-16 | qPSC down |
| KCNK1 | -1.428 | 7.22e-15 | qPSC down |
| KCNN3 | -2.165 | 5.87e-12 | qPSC down |
| KCNMB4 | -1.013 | 2.19e-09 | qPSC down |
| KCNK15 | -1.397 | 1.92e-06 | qPSC down |
| KCND3 | -0.334 | 1.71e-04 | qPSC down |
| KCNC4 | 0.320 | 7.34e-04 | qPSC up |
| KCNA7 | 0.661 | 1.13e-03 | qPSC up |
| KCNT2 | 0.665 | 4.10e-02 | qPSC up |
| KCNA5 | 3.659 | 0 | smPSC up |
| KCNAB1 | 3.441 | 4.65e-233 | smPSC up |
| KCNMB1 | 0.578 | 2.15e-129 | smPSC up |
| KCNJ8 | -1.267 | 9.41e-72 | smPSC down |
| KCNK3 | 2.314 | 9.33e-52 | smPSC up |
| KCNK17 | 0.683 | 3.13e-44 | smPSC up |
| KCND2 | -2.741 | 3.03e-34 | smPSC down |
| KCNK6 | -1.330 | 7.71e-28 | smPSC down |
| KCNMA1 | 1.808 | 1.95e-27 | smPSC up |
| KCNC4 | 1.359 | 5.70e-27 | smPSC up |
| KCND3 | -2.084 | 9.64e-25 | smPSC down |
| KCNK1 | -0.852 | 8.59e-11 | smPSC down |
| KCNS3 | 0.372 | 1.05e-08 | smPSC up |
| KCNT2 | -1.332 | 2.07e-07 | smPSC down |
| KCNN3 | -1.560 | 7.25e-06 | smPSC down |
| KCNK15 | -0.374 | 3.50e-03 | smPSC down |
| KCNJ16 | 0.494 | 2.26e-02 | smPSC up |
| KCNK17 | -4.169 | 3.27e-221 | myCAF down |
| KCNJ8 | -1.991 | 3.41e-161 | myCAF down |
| KCNA5 | -5.488 | 2.45e-132 | myCAF down |
| KCNE4 | -0.893 | 1.93e-63 | myCAF down |
| KCNAB1 | -3.495 | 3.54e-62 | myCAF down |
| KCND2 | 1.938 | 1.24e-38 | myCAF up |
| KCNK3 | -3.291 | 8.27e-37 | myCAF down |
| KCNK1 | 1.219 | 2.92e-16 | myCAF up |
| KCNC4 | -1.343 | 5.29e-15 | myCAF down |
| KCNS3 | -0.569 | 1.08e-13 | myCAF down |
| KCNMA1 | 0.604 | 4.93e-13 | myCAF up |
| KCNMB4 | 0.608 | 2.31e-07 | myCAF up |
| KCNN3 | 1.392 | 1.38e-05 | myCAF up |
| KCNN4 | 0.478 | 3.14e-05 | myCAF up |
| KCNA7 | -0.797 | 4.14e-05 | myCAF down |
| KCNAB2 | 0.574 | 3.97e-04 | myCAF up |
| KCNJ16 | -0.612 | 5.11e-04 | myCAF down |
| KCNK17 | -2.339 | 6.47e-114 | csCAF down |
| KCNA5 | -5.500 | 2.05e-111 | csCAF down |
| KCNAB1 | -3.664 | 2.11e-55 | csCAF down |
| KCND3 | 1.571 | 1.27e-35 | csCAF up |
| KCNE4 | -0.681 | 3.07e-26 | csCAF down |
| KCNC4 | -1.790 | 1.66e-23 | csCAF down |
| KCNJ8 | -0.658 | 2.76e-21 | csCAF down |
| KCNMB1 | 0.616 | 5.29e-20 | csCAF up |
| KCNS3 | -0.757 | 1.70e-16 | csCAF down |
| KCNK15 | 1.654 | 2.15e-16 | csCAF up |
| KCNK3 | -1.444 | 1.83e-14 | csCAF down |
| KCNMA1 | -1.681 | 4.64e-13 | csCAF down |
| KCNN3 | 1.152 | 1.06e-10 | csCAF up |
| KCND2 | 0.657 | 9.12e-07 | csCAF up |
| KCNT2 | 0.529 | 2.79e-06 | csCAF up |
| KCNJ15 | -0.336 | 3.44e-06 | csCAF down |
| KCNN4 | -0.861 | 2.38e-03 | csCAF down |
| KCNK1 | 0.353 | 1.17e-02 | csCAF up |
| KCNJ8 | -2.085 | 4.42e-12 | iCAF down |
| KCNS3 | -2.157 | 5.31e-11 | iCAF down |
| KCNA5 | -6.891 | 2.02e-10 | iCAF down |
| KCNK17 | -2.352 | 4.46e-10 | iCAF down |
| KCNMB1 | -1.605 | 6.21e-07 | iCAF down |
| KCNK6 | -2.683 | 4.34e-06 | iCAF down |
| KCNK1 | 1.534 | 1.10e-05 | iCAF up |
| KCNE4 | -0.332 | 2.69e-04 | iCAF down |
| KCNAB1 | -1.988 | 1.39e-02 | iCAF down |
| KCNAB2 | -1.393 | 1.43e-02 | iCAF down |

## Final KCN Union (25 genes)

- `KCNA5`
- `KCNA7`
- `KCNAB1`
- `KCNAB2`
- `KCNB1`
- `KCNC4`
- `KCND2`
- `KCND3`
- `KCNE4`
- `KCNJ15`
- `KCNJ16`
- `KCNJ8`
- `KCNK1`
- `KCNK15`
- `KCNK17`
- `KCNK3`
- `KCNK6`
- `KCNMA1`
- `KCNMB1`
- `KCNMB4`
- `KCNN3`
- `KCNN4`
- `KCNQ4`
- `KCNS3`
- `KCNT2`

Last updated: 2026-05-25
