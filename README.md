# Processing of single-cell microscopy data of [Guerra et al. 2021]

**CODE: Single_cells_NC_ratio_and_alignment.ipynb**

The provided script is used to calculate the nuclear-to-cytosolic (N/C) ratio of a fluoscently tagged protein in individual budding yeast cells, it also splits each time series into cell cycle traces and aligns them on a relative time scale as described in [Guerra et al. 2021]. 

**Inputs** : 1) a yeast microscopy movie 2) segmentation output from the BudJ plugin of ImageJ 3) a text file with recorded budding and karyokinesis events in a python dictionary format

**Outputs**: 1) N/C ratio calculation for each detected cell at each frame 2) interpolated single-cell cell cycle traces using a fixed number of time points 3) alignment of single-cell traces at specific cell cycle events (bud appearance and karyokinesis).
