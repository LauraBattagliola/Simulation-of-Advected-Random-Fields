# Spatio-Temporal Transport Models

'R' code for *Modeling Spatio-Temporal Transport: From Rigid Advection to Realistic Dynamics*  
M.L. Battagliola and S.C. Olhede 

Preprint available [here](https://arxiv.org/abs/2303.02756)

## Files

| File | Model | Figure |
|------|-------|--------|
| `01_FrozenField_Transport.R` | Frozen Field | Fig. 3 |
| `02a_DistributedFF_Cyclone.R` | Distributed FF (cyclone) | Fig. 1 |
| `02b_DistributedFF_TransportDiffusion.R` | Distributed FF (PDE comparison) | Fig. 4 |
| `03_EvolvingFF_Cyclone.R` | Evolving FF | Fig. 2 |
| `04_DampedFF_SpectralFilter.R` | Damped FF | Fig. 7 |

## Requirements

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "MASS", 
                   "circular", "reshape2", "gridExtra", "ReacTran"))
```

## Contact

- laura.battagliola@itam.mx
- sofia.olhede@epfl.ch
