# Nutrient-switch-analysis
Code for analysis and figures in https://doi.org/10.1101/2024.10.04.616606 <br />
Includes system for categorising and following single-cells and their progeny after a nutrient switch (downstream of CellACDC).
Written in Python using Jupyter notebooks. Details on packages and versions can be found in the [environment.yml](https://github.com/SchmollerLab/Nutrient-switch-analysis/blob/main/environment.yml) file. Microsoft Windows 10 Pro was used as an operating system. 


# General structure
The figures and the corresponding code are grouped into folders based on the underlying experimental datasets used.
Within each folder, the dataprep scripts must be run first, in the sequence assigned in their filenames. 
Once the prepped data is saved, the plotting script for the figures can be run. <br />


While the functions used for dataprep are common between the different replicates of a particular kind of microscopy experiment, unique scripts have been included for every replicate because they carry information about the experimental setup (eg. data location, what position corresponds to what strain etc.). This also allows for saving analysed data for individual replicates before combining them into one file. <br />


Dataprep for figures based on steady-state microscopy [(Fig. 1F)](https://github.com/SchmollerLab/Nutrient-switch-analysis/tree/main/final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy) must be performed before plotting figures based on nutrient-switch microscopy [(Figs. 4, 5 & 7)](https://github.com/SchmollerLab/Nutrient-switch-analysis/tree/main/final_code/Fig4_Fig5_Fig7_FigS4_FigS5). Dataprep for steady-state microscopy uses a custom script instead of functions from Cell-ACDC because this data was analysed before the Cell-ACDC downstream analysis notebook was published. However, the steady-state microscopy data has also been analysed using functions from the Cell-ACDC downstream notebook and both methods of dataprep were found to yield the same results. I recommend using the custom script for dataprepping single-channel microscopy files that were analysed with a very early version of Cell-ACDC (<2021), because the custom script is compatible with their nomenclature. <br />

# Data availability 
The raw microscopy files (.nd2) and the corresponding Cell-ACDC output folders (with manually corrected cell segmentation, cell tracking, cell cycle annotations) are available at ___xyz___. <br />
For each figure, the dates of the experimental replicates (provided in the dataprep code as well as in this README file) can be used to find the underlying dataset. 

# Figure-specific information
  - [Fig1_FigS1_FigS3](https://github.com/SchmollerLab/Nutrient-switch-analysis/tree/main/final_code/Fig1_FigS1_Fig_S3)

  
    - [Fig_1C_S1_smFISH](https://github.com/SchmollerLab/Nutrient-switch-analysis/tree/main/final_code/Fig1_FigS1_Fig_S3/Fig_1C_S1_smFISH) <br />
      - smFISH analysis on datasets 05072021 and 20102022 <br />
      - Two scripts for dataprep of individual replicates [(dataprep_rep1.ipynb](https://github.com/SchmollerLab/Nutrient-switch-analysis/blob/main/final_code/Fig1_FigS1_Fig_S3/Fig_1C_S1_smFISH/dataprep_rep1.ipynb) & [dataprep_rep2.ipynb)](final_code/Fig1_FigS1_Fig_S3/Fig_1C_S1_smFISH/dataprep_rep2.ipynb) and one script for plotting data from both replicates combined [plots_repscombined.ipynb.](final_code/Fig1_FigS1_Fig_S3/Fig_1C_S1_smFISH/plots_repscombined.ipynb)

    - [Fig_1D_smFISH](final_code/Fig1_FigS1_Fig_S3/Fig_1D_smFISH) <br />
      - plotting example images from Position_19 of 05072021 + cell contours and detected spots [smFISH_example_images.ipynb](final_code/Fig1_FigS1_Fig_S3/Fig_1D_smFISH/smFish_example_images.ipynb)

    - [Fig_1F_1G_S3_SteadyStateMicroscopy](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy) <br />
      - Dataprep:  <br />
        - 1 script per medium for steady state microscopy experiments:
          - SCD: datasets 20210120 and 20210121, [1script_SCD.ipynb](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy/dataprep/1script_SCD.ipynb). SCGE: datasets 20210923 and 20230103, [1script_SCGE.ipynb](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy/dataprep/1script_SCGE.ipynb)
        - script2 [(2script_combined_media.ipynb)](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy/dataprep/2script_combined_media.ipynb) and script3 [(3script_combined_media.ipynb)](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy/dataprep/3script_combined_media.ipynb) for combining the media and further downstream analysis  
        - script4 [(4script_additionalCcr4rep.ipynb)](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy/dataprep/4script_additionalCcr4rep.ipynb), script5 [(5script_additionalCcr4rep.ipynb)](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy/dataprep/5script_additionalCcr4rep.ipynb) and script6 [(6script_additionalCcr4rep.ipynb)](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy/dataprep/6script_additionalCr4rep.ipynb) for integrating an extra replicate of SCD microscopy for the Ccr4 mutant (dataset 20231212)
      - 2 scripts for plotting: [plot_Fig1F_and_G.ipynb](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy/plot_Fig1F_and_G.ipynb) and [plot_FigS3.ipynb](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy/plot_FigS3.ipynb)

    - Fig. 1B was plotted from qPCR data in Microsoft Excel. 

 - [Fig2](final_code/Fig2)

   
   - Dataprep: example script for dataprep for one replicate [(concatenating_timept_files_per_rep.ipynb)](final_code/Fig2/dataprep/concatenating_timept_files_per_rep.ipynb). This script concatenates coulter counter data from multiple timepoints and cleans it for use in python. Run this script for datasets 04082022, 05092022, 06092022, 07072022, 09062022 and 15062022. Combine the outputs for all the reps and import resulting df into plotting script.
   - 1 script for plotting Fig. 2B, 2C and 2D [(plot_Fig2.ipynb)](final_code/Fig2/plot_Fig2.ipynb)



 - [Fig4_Fig5_Fig7_FigS4_FigS5](final_code/Fig4_Fig5_Fig7_FigS4_FigS5)

   - Dataprep for nutrient-switch experiments - split into three further categories:
     - [nutrient_switch_original4strains](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_original4strains) - containing WT, *Δwhi5*, *Δbck2* and *Δwhi5Δbck2* <br />
        datasets 09082022, 13102022 and 28032023 - analysed with respective scripts [(script1_rep1_09082022.ipynb](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_original4strains/script1_rep1_09082022.ipynb), [script2_rep2_13102022.ipynb](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_original4strains/script2_rep2_13102022.ipynb), [script3_rep3_28032023.ipynb)](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_original4strains/script3_rep3_28032023.ipynb)<br />
        +1 script to combine the replicates [script4_repscombined_withsteadystate.ipynb](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_original4strains/script4_repscombined_withsteadystate.ipynb)
     - [nutrient_switch_cln3deletions](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_cln3deletions) - containing WT, *Δcln3*, *Δwhi5Δcln3*, *Δwhi5Δcln3Δbck2*, MK40, MK41 and MK44. MK40, MK41 and MK44 were not used in the paper. <br />
        datasets 06042023 and 12042023 - analysed with respective scripts [(script1_rep1_06042023.ipynb](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_cln3deletions/script1_rep1_06042023.ipynb), [script2_rep2_12042023.ipynb](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_cln3deletions/script2_rep2_12042023.ipynb)<br />
        +1 script to combine the replicates [script3_repscombined.ipynb](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_cln3deletions/script3_repscombined.ipynb)
     - [nutrient_switch_cln3stabilisations_with_deletion_mutants](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_cln3stabilisations_with_deletion_mutants) - containing MK40, *Δwhi5_MK40bg*, *Δbck2_MK40bg* and *Δwhi5Δbck2_MK40bg*, MK44, *Δwhi5_MK44bg*, *Δbck2_MK44bg* and *Δwhi5Δbck2_MK44bg* <br /> These strains were not used in the paper. <br />
        datasets 28092023 and 06102023 - analysed with respective scripts [(script1_rep1_28092023.ipynb](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_cln3stabilisations_with_deletion_mutants/script1_rep1_28092023.ipynb), [script2_rep2_06102023.ipynb](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_cln3stabilisations_with_deletion_mutants/script2_rep2_06102023.ipynb)<br />
        +1 script to combine the replicates [nutrient_switch_repscombined_Cln3stabs_withMutants.ipynb](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/dataprep/nutrient_switch_cln3stabilisations_with_deletion_mutants/nutrient_switch_repscombined_Cln3stabs_withMutants.ipynb)

   - Plotting script - [plot_Fig4_Fig5_Fig7_FigS4_FigS5](final_code/Fig4_Fig5_Fig7_FigS4_FigS5/plot_Fig4_Fig5_Fig7_FigS4_FigS5.ipynb). NOTE: Please make sure you have run dataprep for [Fig_1F_1G_S3_SteadyStateMicroscopy](final_code/Fig1_FigS1_Fig_S3/Fig_1F_1G_S3_SteadyStateMicroscopy) before running this script to include steady-state control cells.


- [Fig6_FigS6](final_code/Fig6_FigS6)

  - Dataprep for nutrient switch with cell-cycle reporter strains: 
    - script1 and script2 for dataprep of each replicate (datasets 13102023 and 19102023): [script1_rep1_13102023.ipynb](final_code/Fig6_FigS6/dataprep/script1_rep1_13102023.ipynb) and [script2_rep2_19102023.ipynb](final_code/Fig6_FigS6/dataprep/script2_rep2_19102023.ipynb)
    - script3 [(script3_fluorescenceanalysis.ipynb)](final_code/Fig6_FigS6/dataprep/script3_fluorescenceanalysis.ipynb),script4 [(script4_fluorescenceanalysis_II.ipynb)](final_code/Fig6_FigS6/dataprep/script4_fluorescenceanalysis_II.ipynb) and script5 [(script5.ipynb)](final_code/Fig6_FigS6/dataprep/script5.ipynb) for fluorescence analysis and combining it with nutrient-switch analysis

  - Plotting [Fig6_FigS6_plots.ipynb](final_code/Fig6_FigS6/Fig6_FigS6_plots.ipynb)

- [FigS2](final_code/FigS2)

  - 1 script for dataprep and plotting steady-state SCD coulter counter data for WT, *Δccr4* and *Δwhi5Δccr4* from dataset 04042023 [FigS2_dataprep_and_plotting.ipynb](final_code/FigS2/FigS2_dataprep_and_plotting.ipynb)
