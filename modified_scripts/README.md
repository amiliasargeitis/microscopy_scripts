This folder contains a modification of the Single_cells_NC_ratio_and_alignment.ipynb script presented in [Guerra et al. 2022]. The modified script was used in [Koutsoumpa et al. 2026].
Similarly to the original script, the modified script (General_Sfp1_ratio_aligned_cc.py) calculates the nuclear-to-cytosolic (N/C) ratio of a fluorescently tagged protein (Sfp1) in individual budding yeast cells. The pipeline also splits each time series into cell-cycle traces using annotated cell cycle events and aligns them on a relative time scale.

The original workflow of [Guerra et al. 2022] used a nuclear fluorescent marker (Hta2-mRFP) to facilitate the segmentation of the nucleus. Since no nuclear marker was available in the datasets of [Koutsoumpa et al. 2026], the segmentation step was modified by using an adaptive local thresholding approach. The local threshold offset was iteratively adjusted until the detected nucleus area corresponded to approximately 30% of the total cell area, which was then used as the nuclear mask for downstream N/C ratio calculations.

The script requires the following inputs:
-	A yeast microscopy movie
-	Segmentation output from the BudJ plugin of ImageJ
-	A text file containing recorded birth, budding and cytokinesis events from the Click_cells plugin implemented on ImageJ

The script produces the following outputs:
-	N/C ratio calculation for each detected cell at each frame
-	Interpolated single-cell cell-cycle traces using a fixed number of time points
-	Alignment of single-cell traces at specific cell-cycle events (bud appearance and cytokinesis)

 
