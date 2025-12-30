# Spatio-Temporal Transport Models

'R' code for *Modeling Spatio-Temporal Transport: From Rigid Advection to Realistic Dynamics*  
M.L. Battagliola and S.C. Olhede 

Preprint available [here](https://arxiv.org/abs/2303.02756).

## Mathematical Background

Let $$X_S(s)$$ be a spatial random field with covariance function $$c_XX(h)$$ and spectral density $$S_{XX}(k)$$, where $$s, h \in \mathbb{R}^2$$ are spatial coordinates and $$k \in \mathbb{R}^2$$ is the wavenumber.

The models below define space-time random fields Z(s,t) by advecting X_S in different ways:

### Frozen Field (FF)
A spatial field advected by constant velocity v:

$$Z(s,t) = X_S(s - vt)$$

Covariance: $c_{ZZ}(h,\tau) = c_{XX}(h - v\tau)$  
Spectrum: $S_{ZZ}(k,\omega) = S_{XX}(k) · \delta(\omega + k^\top v)$

### Distributed FF
A weighted superposition of frozen fields with different velocities $\{ v_i\}_i$:

$$Z(s,t) = \sum_{i} p_i X_S(s - v_it)$$

The velocity spread introduces diffusive behavior and breaks perfect temporal correlation.

### Evolving FF
A frozen field with spatially and temporally varying velocity:

$$Z(s,t) = X_S(s - v(s,t)·t)$$

Allows local deformations like rotation or stretching.

### Damped FF
A spectral generalization with tunable temporal decorrelation:

$$S_{ZZ}(k,\omega) = S_{XX}(k) · B(k,\omega; v,\beta)^{-\alpha}$$

where $B(k,\omega; v,\beta)^{-\alpha} = (\omega + k^\top v)^2 + (\beta \omega)^2$, with persistence parameter $\alpha>1/2$ and damping parameter $\beta \geq 0$.

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

