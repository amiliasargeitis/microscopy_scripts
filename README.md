# Processing of single-cell microscopy data

The provided script is used to calculate the nuclear-to-cytosolic (N/C) ratio of a fluoscently tagged protein in individual budding yeast cells, as described in [Guerra et al. 2021]. 

**Inputs** : 1) a yeast microscopy movie 2) segmentation output from the BudJ plugin of ImageJ

**Outputs**: 1) N/C ratio calculation for each detected cell at each frame 2) interpolated single-cell traces using a fixed number of time points 3) alignment of single-cell traces at specific cell cycle events (bud appearance and kayrokinesis).
