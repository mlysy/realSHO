**(R)obust and (E)fficient c(AL)ibration of (SHO) parameters**


*Bryan Yates, Martin Lysy, Aleksander Labuda*

*June 22, 2017*

---

*Description:*

Provides optimization routines to calibrate the parameters of a simple harmonic oscillator (SHO) model describing the tip displacement of an atomic force microscope (AFM) cantilever.  The most effective of these combines the simplicity of nonlinear least-squares with the statistical efficiency of maximum likelihood.  Also provides a routine for removing most electronic noise in an automated pre-processing step.

*Contents:*

- `./data`: `Q_100_PSD_noise.mat` - Simulated PSD with added sine-noise
- `./functions`: Folder including all relevant functions


*Usage:*

Download all package contents to a common folder and run the examples.


*Examples:*

- `tutorial_basic.m`: PSD generation from time domain & SHO+White Noise (SHOW) fitting
- `tutorial_advanced.m`: SHOW Fitting + noise removal techniques
- `tutorial_CI.m`: Confidence interval calculations

