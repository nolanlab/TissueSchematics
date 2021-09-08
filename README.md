# TissueSchematics
Code accompanying the manuscript "Tissue Schematics map the specialization of immune tissue motifs and their appropriation by tumors" by Bhate, Barlow et al.

There are three subfolders:

CellHier:
This contains code for identifying CNs and visualization. This was originally published in Schuerch et al. Cell 2020 and can also be found in the  NeighborhoodCoordination repo.

ipynb:
This contains the ipython notebooks used to perform the analyses in the paper:
- cellularneighborhoods (figure 1)
- spatial_contexts (figure 2)
- assembly_rules (figure 3)
- Figure4 - CRC spatial contexts
- Figure5 - motifs etc in CRC, save MCMC analysis for triangles
- mcmc_generating_samples - code used to call C++ library for generating samples of colorings from the Level 1 null sets
- mcmc_validation_analysis - code used to find higher-order triangles and validate MCMC mixing
- colormc.cpp - C++ library for generating samples of colorings from the Level 1 null sets

htmls:
 - Contains htmls of the crc analysis notebooks to look through, with all the analysis for each figure annotated.
