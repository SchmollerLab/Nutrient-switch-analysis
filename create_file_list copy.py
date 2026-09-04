from pathlib import Path
from collections import Counter
import os
import re
import warnings

import pandas as pd

# aligned, background ROI etc should always point source image to the correct source tif
custom_source_image = {
    "_aligned.npz",
    "_bkgrRoiData.npz",
    "_aligned_bkgrRoiData.npz",
    "_tracked.npz"
}

excluded_experiment_folder_names = {"unused position", "test", "unused positions", "recovery", "unused_positions", "unused_position", "unused_positions"}

# Source-image patterns for derived analysis tables. Dated files map to every
# position image from the corresponding acquisition; summary tables do not
# identify one unique source image from their filenames.
analysis_table_source_images = {
    "20231013_20231019_repscombined/final_processed_reps_combined_fluorescence_analysis_dataset.csv": [
        "20231013/Position_*/Images/*_Dia_Ph3_mCherry.tif",
        "20231019/Position_*/Images/*_Dia_Ph3_mCherry.tif",
    ],
    "20231013_20231019_repscombined/reps_combined_fluorescence_analysis_dataset.csv": [
        "20231013/Position_*/Images/*_Dia_Ph3_mCherry.tif",
        "20231019/Position_*/Images/*_Dia_Ph3_mCherry.tif",
    ],
    "bootstrapped_cv_df_2.csv": [
        "20210120/TIFFs/Position_*/Images/*_phase_contr.tif",
        "20210121/MIA_SCD_20210121_6strains.nd2/Position_*/Images/*_phase_contr.tif",
        "20210923/23092021 SCGE/Position_*/Images/*_Dia_Ph3.tif",
        "20220908/Position_*/Images/*_Dia_Ph3.tif",
        "20221013/Position_*/Images/*_Dia_Ph3.tif",
        "20230328/Position_*/Images/*_Dia_Ph3.tif",
        "20230406/Position_*/Images/*_Dia_Ph3.tif",
        "20230412/Position_*/Images/*_Dia_Ph3.tif",
        "20230928/Position_*/Images/*_Dia_Ph3.tif",
        "20231006/Position_*/Images/*_Dia_Ph3.tif",
    ],
    "figure_plotting_df_1.csv": [
        "20210120/TIFFs/Position_*/Images/*_phase_contr.tif",
        "20210121/MIA_SCD_20210121_6strains.nd2/Position_*/Images/*_phase_contr.tif",
        "20210923/23092021 SCGE/Position_*/Images/*_Dia_Ph3.tif",
        "20220908/Position_*/Images/*_Dia_Ph3.tif",
        "20221013/Position_*/Images/*_Dia_Ph3.tif",
        "20230328/Position_*/Images/*_Dia_Ph3.tif",
        "20230406/Position_*/Images/*_Dia_Ph3.tif",
        "20230412/Position_*/Images/*_Dia_Ph3.tif",
        "20230928/Position_*/Images/*_Dia_Ph3.tif",
        "20231006/Position_*/Images/*_Dia_Ph3.tif",
    ],
    "06042023_jupyter_overall_df_with_rel_switch_and_wave_info.csv": "20230406/Position_*/Images/*_Dia_Ph3.tif",
    "06102023_jupyter_overall_df_with_rel_switch_and_wave_info.csv": "20231006/Position_*/Images/*_Dia_Ph3.tif",
    "07032025_jupyter_overall_df_with_rel_switch_and_wave_info.csv": "07032025/Position_*/Images/*_Dia_Ph3_mCherry.tif",
    "09082022_jupyter_overall_df_with_rel_switch_and_wave_info.csv": "20220908/Position_*/Images/*_Dia_Ph3.tif",
    "11122024_jupyter_overall_df_with_rel_and_switch_and_wave_info.csv": "11122024/Position_*/Images/*_Dia_Ph3_mCherry.tif",
    "12042023_jupyter_overall_df_with_rel_switch_and_wave_info.csv": "20230412/Position_*/Images/*_Dia_Ph3.tif",
    "13102022_jupyter_overall_df_with_rel_switch_and_wave_info.csv": "20221013/Position_*/Images/*_Dia_Ph3.tif",
    "13102023_jupyter_overall_df_with_rel_and_switch_and_wave_info.csv": "20231013/Position_*/Images/*_Dia_Ph3_mCherry.tif",
    "19102023_jupyter_overall_df_with_rel_and_switch_and_wave_info.csv": "20231019/Position_*/Images/*_Dia_Ph3_mCherry.tif",
    "27032025_jupyter_overall_df_with_rel_switch_and_wave_info.csv": "27032025/20250327_SCD_2_SCGE_25_MMY_1-5_YCY005_6-10_KSY15_11-15_YCY41_16-20_YCY44_21-25_KSY41_26-30_YCY42_31-35_YCY45_36-40/Position_*/Images/*_Dia_Ph3_mCherry.tif",
    "28032023_jupyter_overall_df_with_rel_switch_and_wave_info.csv": "20230328/Position_*/Images/*_Dia_Ph3.tif",
    "28092023_jupyter_overall_df_with_rel_switch_and_wave_info.csv": "20230928/Position_*/Images/*_Dia_Ph3.tif",
    "ccr4_optimised4MLR_SCD_SCGE_merged_dataset.csv": "Ccr4_data_for_Felix/20231212/Position_*/Images/*_Dia_Ph3.tif",
    "Cln3delsstabs_jupyter_overall_df_with_rel_switch_and_wave_info.csv": [
        "20230406/Position_*/Images/*_Dia_Ph3.tif",
        "20230412/Position_*/Images/*_Dia_Ph3.tif",
    ],
    "Cln3stabs_withMutants_jupyter_overall_df_with_rel_switch_and_wave_info.csv": [
        "20230928/Position_*/Images/*_Dia_Ph3.tif",
        "20231006/Position_*/Images/*_Dia_Ph3.tif",
    ],
    "CombinedCcr4Rep2 SCD LCM analysed.csv": "Ccr4_data_for_Felix/20231212/Position_*/Images/*_Dia_Ph3.tif",
    "CombinedTRs SCD LCM analysed.csv": [
        "20210120/TIFFs/Position_*/Images/*_phase_contr.tif",
        "20210121/MIA_SCD_20210121_6strains.nd2/Position_*/Images/*_phase_contr.tif",
    ],
    "CombinedTRs SCD LCM data.csv": [
        "20210120/TIFFs/Position_*/Images/*_phase_contr.tif",
        "20210121/MIA_SCD_20210121_6strains.nd2/Position_*/Images/*_phase_contr.tif",
    ],
    "CombinedTRs SCGE LCM analysed.csv": [
        "20210923/23092021 SCGE/Position_*/Images/*_Dia_Ph3.tif",
        "20230103/Position_*/Images/*_Dia_Ph3.tif",
    ],
    "CombinedTRs SCGE LCM data.csv": [
        "20210923/23092021 SCGE/Position_*/Images/*_Dia_Ph3.tif",
        "20230103/Position_*/Images/*_Dia_Ph3.tif",
    ],
    "complete_optimised4MLR_SCD_SCGE_merged_dataset.csv": [
        "20210120/TIFFs/Position_*/Images/*_phase_contr.tif",
        "20210121/MIA_SCD_20210121_6strains.nd2/Position_*/Images/*_phase_contr.tif",
        "20210923/23092021 SCGE/Position_*/Images/*_Dia_Ph3.tif",
        "20230103/Position_*/Images/*_Dia_Ph3.tif",
        "Ccr4_data_for_Felix/20231212/Position_*/Images/*_Dia_Ph3.tif",
    ],
    "CoulterCounter_postAnalysis_04042023_sizeRangesandRepsCombined.csv": "Coulter counter measurements",
    "final_10_reps_combined.csv": "Coulter counter measurements",
    "MetprCln311Arepscombined_jupyter_overall_df_with_rel_switch_and_wave_info.csv": [
        "07032025/Position_*/Images/*_Dia_Ph3_mCherry.tif",
        "27032025/20250327_SCD_2_SCGE_25_MMY_1-5_YCY005_6-10_KSY15_11-15_YCY41_16-20_YCY44_21-25_KSY41_26-30_YCY42_31-35_YCY45_36-40/Position_*/Images/*_Dia_Ph3_mCherry.tif",
    ],
    "normalnutrientswitch_jupyter_overall_df_with_rel_switch_and_wave_info.csv": [
        "20220908/Position_*/Images/*_Dia_Ph3.tif",
        "20221013/Position_*/Images/*_Dia_Ph3.tif",
        "20230328/Position_*/Images/*_Dia_Ph3.tif",
    ],
    "optimised4MLR_SCD_SCGE_merged_dataset.csv": [
        "20210120/TIFFs/Position_*/Images/*_phase_contr.tif",
        "20210121/MIA_SCD_20210121_6strains.nd2/Position_*/Images/*_phase_contr.tif",
        "20210923/23092021 SCGE/Position_*/Images/*_Dia_Ph3.tif",
        "20230103/Position_*/Images/*_Dia_Ph3.tif",
    ],
    "smFISH_rep1_postAnalysis.csv": "smFISH analysis",
    "smFISH_rep2_postAnalysis.csv": "smFISH analysis",
    "steadystate_analysed_fornutswitch_jupyter_overall_df_with_rel_switch_and_wave_info.csv": [
        "20210120/TIFFs/Position_*/Images/*_phase_contr.tif",
        "20210121/MIA_SCD_20210121_6strains.nd2/Position_*/Images/*_phase_contr.tif",
        "20210923/23092021 SCGE/Position_*/Images/*_Dia_Ph3.tif",
    ],
}

file_suffix_types = {
    "_acdc_output_all_segm_masks_combined_metrics.csv": "acdc output",
    "_acdc_output_nuc_tracked_cytoplasm.csv": "acdc output",
    "_acdc_output_nuc_tracked.csv": "acdc output",
    "_segm_nuc_tracked_cytoplasm.npz": "nucleus segmentation",
    "_segm_nuc_tracked.npz": "nucleus segmentation",
    "_Dia_Ph3_mCherry_aligned_bkgrRoiData.npz": "background ROI data",
    "_mCitrine_Yagya_aligned_bkgrRoiData.npz": "background ROI data",
    "_mScarlet_Ph3_Yagya_aligned_bkgrRoiData.npz": "background ROI data",
    "_Dia_Ph3_aligned_bkgrRoiData.npz": "background ROI data",
    "_acdc_output.csv": "acdc output",
    "_aligned_bkgrRoiData.npz": "background ROI data",
    "_dataPrepROIs_coords.csv": "ROI coordinates",
    "_dataPrep_bkgrROIs.json": "background ROI settings",
    "_cca_properties_downstream.csv": "downstream CCA properties",
    "_equations_all_segm_masks_combined_metrics.ini": "metric equations",
    "_equations_cell_nuc_combined_metrics.ini": "metric equations",
    "_custom_combine_metrics.ini": "custom metric settings",
    "_custom_annot_params.json": "annotation settings",
    "_segm_hyperparams.ini": "segmentation settings",
    "_tracked_lost_centroids.json": "tracking data",
    "_aligned_bkgrRoiData.npz": "background ROI data",
    "_Dia_Ph3_mCherry.tif": "phase contrast image",
    "_Dia_Ph3.tif": "phase contrast image",
    "_phase_contr.tif": "phase contrast image",
    "_mCitrine_Yagya.tif": "fluorescence image",
    "_mScarlet_Ph3_Yagya.tif": "fluorescence image",
    "_metadataXML.txt": "metadata",
    "_metadata.csv": "metadata",
    "_align_shift.npy": "alignment data",
    "_align_shift.npz": "alignment data",
    "_aligned.npz": "aligned image data",
    "_phc_aligned.npy": "aligned image data",
    "_phc_aligned.npz": "aligned image data",
    "_segm_hyperparams.ini": "segmentation settings",
    "_segmInfo.csv": "segmentation information",
    "_segm_nuc.npz": "nucleus segmentation",
    "_segm.npy": "segmentation",
    "_segm.npz": "segmentation",
    "_cc_stage.csv": "cell-cycle stage",
    "_delROIsInfo.npz": "ROI data",
    "_dataPrep_bkgrValues.csv": "background values",
    "_last_tracked_i.txt": "tracking data",
}

custom_raw_to_annotation_mapper = {
    "060423": "06042023",
    "120423": "12042023",
    "20231212": "Ccr4_data_for_Felix/20231212",
    "20250307": "07032025",
}

file_annotation_dict = {
    "Ccr4_data_for_Felix/20231212": {
        "Date": "12122023",
        "Condition": "SCD",
        "source_image": "_Dia_Ph3.tif",
        "Position": "#Position",
        "position_to_strain": {
            position: "YCY002_1" for position in range(1, 33)
        },
        "experiment": "ccr4_deletion",
    },
    "07032025": {
        "Date" : "07032025",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3_mCherry.tif",
        "Position": "#Position",
        "position_to_strain": {
            1: "MMY",
            2: "MMY",
            3: "MMY",
            4: "MMY",
            5: "MMY",
            6: "YCY005_1_whi5bck2doubledel",
            7: "YCY005_1_whi5bck2doubledel",
            8: "YCY005_1_whi5bck2doubledel",
            9: "YCY005_1_whi5bck2doubledel",
            10: "YCY005_1_whi5bck2doubledel",
            11: "KSY015_5_METpr_CLN",
            12: "KSY015_5_METpr_CLN",
            13: "KSY015_5_METpr_CLN",
            14: "KSY015_5_METpr_CLN",
            15: "KSY015_5_METpr_CLN",
        },
        "switch": 50,
        "experiment": "metpr_cln211a",
    },
    #
    "11122024": {
        # "Figure" : "6",
        "Date" : "11122024",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3_mCherry.tif",
        "Position": "#Position",
        "_special_mapper": {
            "_segm_nuc.npz": "_mScarlet_Ph3_Yagya.tif",
        },
        "position_to_strain": {
            1: "YCY024_1_WTlabelled",
            2: "YCY024_1_WTlabelled",
            3: "YCY024_1_WTlabelled",
            4: "YCY024_1_WTlabelled",
            5: "YCY024_1_WTlabelled",
            6: "YCY024_1_WTlabelled",
            7: "YCY024_1_WTlabelled",
            8: "YCY030_2_Whi5del_labelled",
            9: "YCY030_2_Whi5del_labelled",
            10: "YCY030_2_Whi5del_labelled",
            11: "YCY030_2_Whi5del_labelled",
            12: "YCY030_2_Whi5del_labelled",
            13: "YCY030_2_Whi5del_labelled",
            14: "YCY030_2_Whi5del_labelled",
            15: "YCY037_1_Bck2del_labelled",
            16: "YCY037_1_Bck2del_labelled",
            17: "YCY037_1_Bck2del_labelled",
            18: "YCY037_1_Bck2del_labelled",
            19: "YCY037_1_Bck2del_labelled",
            20: "YCY037_1_Bck2del_labelled",
            21: "YCY037_1_Bck2del_labelled",
            22: "YCY038_1_Whi5Bck2doubledel_labelled",
            23: "YCY038_1_Whi5Bck2doubledel_labelled",
            24: "YCY038_1_Whi5Bck2doubledel_labelled",
            25: "YCY038_1_Whi5Bck2doubledel_labelled",
            26: "YCY038_1_Whi5Bck2doubledel_labelled",
            27: "YCY038_1_Whi5Bck2doubledel_labelled",
            28: "YCY038_1_Whi5Bck2doubledel_labelled",
        },
        "switch": 50,
        "experiment": "cc_marker"
    },
    "20210120": {
        # "Figure" : "1,4,5",
        "Date" : "20012021",
        "Condition" : "SCD",
        "source_image": "_phase_contr.tif",
        "Position": "#Position",
        "_PosFolder": "TIFFs",
        "position_to_strain": {
            **{p: "MMY" for p in range(1, 7)},
            **{p: "BCK2DEL" for p in range(7, 12)},
            **{p: "WHI5DEL" for p in range(12, 20)},
            **{p: "WHI5_BCK2_DD" for p in range(20, 29)},
            **{p: "WHI5_CCR4_DD" for p in range(29, 38)},
        },
        "experiment": "steady_state"
    },
    "20210121": {
        # "Figure" : "1,4,5",
        "Date" : "21012021",
        "Condition" : "SCD",
        "source_image": "_phase_contr.tif",
        "Position": "#Position",
        "_PosFolder": "MIA_SCD_20210121_6strains.nd2",
        "position_to_strain": {
            **{p: "MMY" for p in range(1, 7)},
            **{p: "BCK2DEL" for p in range(7, 12)},
            **{p: "CCR4DEL" for p in range(12, 16)},
            **{p: "WHI5DEL" for p in range(16, 23)},
            **{p: "WHI5_BCK2_DD" for p in range(23, 28)},
            **{p: "WHI5_CCR4_DD" for p in range(28, 34)},
        },
        "experiment": "steady_state"
    },
    "20210923": {
        # "Figure" : "1,4,5",
        "Date" : "23092021",
        "Condition" : "SCGE",
        "source_image": "_Dia_Ph3.tif",
        "Position": "#Position",
        "_PosFolder": "23092021 SCGE",
        "position_to_strain": {
            **{p: "MMY" for p in range(1, 8)},
            **{p: "BCK2DEL" for p in range(8, 15)},
            **{p: "CCR4DEL" for p in range(15, 21)},
            **{p: "WHI5DEL" for p in (22, 23, 24, 26, 27, 28)},
            **{p: "WHI5_BCK2_DD" for p in range(29, 36)},
            **{p: "WHI5_CCR4_DD" for p in (36, 37, 38, 39, 41, 42)},
        },
        "experiment": "steady_state"
        
    },
    "20220908": {
        "Date" : "08092022",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "MMY" for p in range(1, 10)},
            **{p: "BCK2DEL" for p in (10, 11, 13, 14, 15, 16, 17, 18, 19)},
            **{p: "WHI5DEL" for p in (20, 22, 24, 25, 26, 28)},
            **{p: "WHI5_BCK2_DD" for p in range(29, 41)},
        },
        "switch": 50,
        "experiment": "nutrient_switch"
    },
    "20221013": {
        "Date" : "13102022",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "MMY" for p in range(1, 10)},
            **{p: "BCK2DEL" for p in range(10, 23)},
            **{p: "WHI5DEL" for p in (23, 24, 25, 26, 29)},
            **{p: "WHI5_BCK2_DD" for p in (30, 31, 32, 33, 34, 36, 37)},
        },
        "switch": 50,
        "experiment": "nutrient_switch"
    },
    "20230103": { # is equivalent to 3012023
        "Date" : "03012023", 
        "Condition" : "SCGE",
        "source_image": "_Dia_Ph3.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "MMY" for p in range(1, 8)},
            **{p: "BCK2DEL" for p in range(8, 14)},
            **{p: "CCR4DEL" for p in (15, 16, 17, 20, 21)},
            **{p: "WHI5DEL" for p in range(22, 29)},
            **{p: "WHI5_BCK2_DD" for p in range(29, 36)},
        },
        "experiment": "steady_state"
    },
    "20230328": {
        "Date" : "28032023",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "MMY" for p in range(1, 11)},
            **{p: "BCK2DEL" for p in (11, 12, 14, 15, 16, 17, 18, 20)},
            **{p: "WHI5DEL" for p in range(21, 31)},
            **{p: "WHI5_BCK2_DD" for p in range(31, 41)},
        },
        "switch": 50,
        "experiment": "nutrient_switch"
    },
    "20230406": {
        "Date" : "06042023",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "MMY" for p in range(1, 7)},
            **{p: "CLN3DEL" for p in range(7, 13)},
            **{p: "CLN3_WHI5_DD" for p in range(13, 19)},
            **{p: "CLN3_WHI5_BCK2_TD" for p in range(19, 25)},
            **{p: "MK_40" for p in range(25, 31)},
            **{p: "MK_41" for p in range(31, 37)},
            **{p: "MK_44" for p in range(37, 43)},
        },
        "switch": 53,
        "experiment": "cln3_deletion"
    },
    "20230412": {
        "Date" : "12042023",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "MMY" for p in range(1, 7)},
            **{p: "CLN3DEL" for p in range(7, 13)},
            **{p: "CLN3_WHI5_DD" for p in range(13, 19)},
            **{p: "CLN3_WHI5_BCK2_TD" for p in range(19, 23)},
            **{p: "MK_40" for p in range(23, 30)},
            **{p: "MK_41" for p in range(30, 36)},
            **{p: "MK_44" for p in range(36, 42)},
        },
        "switch": 50,
        "experiment": "cln3_deletion"
    },
    "20230928": {
        "Date" : "28092023",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "MK40" for p in range(1, 6)},
            **{p: "YCY027_1_bck2delMK40bg" for p in range(6, 11)},
            **{p: "YCY031_4_whi5delMK40bg" for p in range(11, 16)},
            **{p: "YCY034_1_doubledelMK40bg" for p in range(16, 21)},
            **{p: "MK_44" for p in range(21, 26)},
            **{p: "YCY029_3_bck2delMK44bg" for p in range(26, 31)},
            **{p: "YCY033_2_whi5delMK44bg" for p in range(31, 36)},
            **{p: "YCY036_1_doubledelMK44bg" for p in range(36, 41)},
        },
        "switch": 50,
        "experiment": "cln3_stabilized_deletion"
    },
    "20231006": {
        "Date" : "06102023",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "MK40" for p in range(1, 6)},
            **{p: "YCY027_1_bck2delMK40bg" for p in range(6, 11)},
            **{p: "YCY031_4_whi5delMK40bg" for p in range(11, 16)},
            **{p: "YCY034_1_doubledelMK40bg" for p in range(16, 21)},
            **{p: "MK_44" for p in range(21, 26)},
            **{p: "YCY029_3_bck2delMK44bg" for p in range(26, 31)},
            **{p: "YCY033_2_whi5delMK44bg" for p in range(31, 36)},
            **{p: "YCY036_1_doubledelMK44bg" for p in range(36, 41)},
        },
        "switch": 50,
        "experiment": "cln3_stabilized_deletion"
    },
    "20231013": {
        "Date" : "13102023",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3_mCherry.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "YCY024_1_WTlabelled" for p in range(1, 8)},
            **{p: "YCY030_2_Whi5del_labelled" for p in range(8, 15)},
            **{p: "YCY037_1_Bck2del_labelled" for p in range(15, 22)},
            **{p: "YCY038_1_Whi5Bck2doubledel_labelled" for p in range(22, 29)},
        },
        "switch": 50,
        "experiment": "cc_marker",
        "_special_mapper": {
            "_segm_nuc.npz": "_mScarlet_Ph3_Yagya.tif",
        },
    },
    "20240403": {
        "Date" : "03042024",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3_mCherry.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "JE103" for p in range(2, 7)},
            **{p: "YCY039_3_Whi5del_in_JE103" for p in (12, 13, 14, 15, 19)},
            **{p: "YCY040_1_Bck2del_in_JE103" for p in range(23, 28)},
        },
        "switch": 50,
        "experiment": "JE103_test",
    },
    "27032025": {
        "Date" : "27032025",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3_mCherry.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "MMY" for p in range(1, 6)},
            **{p: "YCY005_1_whi5bck2doubledel" for p in range(6, 11)},
            **{p: "KSY015_5_METpr_CLN3" for p in range(11, 16)},
            **{p: "YCY041_2_METpr_CLN3_bck2del" for p in range(16, 21)},
            **{p: "YCY044_2_METpr_CLN3_whi5bck2doubledel" for p in range(21, 26)},
            **{p: "KSY041_2_METpr_CLN3_11A" for p in range(26, 31)},
            **{p: "YCY042_1_METpr_CLN3_11A_bck2del" for p in range(31, 36)},
            **{p: "YCY045_3_METpr_CLN3_11A_whi5bck2doubledel" for p in range(36, 41)},
        },
        "switch": 50,
        "experiment": "metpr_cln211a",
        "_PosFolder": r"\20250327_SCD_2_SCGE_25_MMY_1-5_YCY005_6-10_KSY15_11-15_YCY41_16-20_YCY44_21-25_KSY41_26-30_YCY42_31-35_YCY45_36-40"
    },
    "20231019": {
        "Date" : "19102023",
        "Condition" : "SCDtoSCGE",
        "source_image": "_Dia_Ph3_mCherry.tif",
        "Position": "#Position",
        "position_to_strain": {
            **{p: "YCY024_1_WTlabelled" for p in range(7, 13)},
            **{p: "YCY030_2_Whi5del_labelled" for p in range(13, 19)},
            **{p: "YCY037_1_Bck2del_labelled" for p in range(19, 25)},
            **{p: "YCY038_1_Whi5Bck2doubledel_labelled" for p in range(25, 31)},
        },
        "switch": 50,
        "experiment": "cc_marker",
        "_special_mapper": {
            "_segm_nuc.npz": "_mScarlet_Ph3_Yagya.tif",
        },
    }
}

file_annotation_dict_by_date = {
    annotation.get("Date"): annotation
    for annotation in file_annotation_dict.values()
    if annotation.get("Date") is not None
}


root_folder = r"G:\Yagya Data\Analysed Live cell microscopy data"
root_raw_microscopy_path = r"G:\Yagya Data\Raw Live cell microscopy data"
root_fish_path = r"G:\Yagya Data\Yagya\FISH"
root_qpcr_path = r"G:\Yagya Data\Yagya\qPCR"

remote_data_folders = {
    "analysed": "analysed_live_cell_microscopy",
    "raw": "raw_microscopy_files",
    "analysis_tables": "analysis_tables",
    "fish": "fish",
    "qpcr": "qpcr",
}


def get_remote_path(data_kind, relative_path):
    """Return a server-side path below the shared New data root."""
    return (Path(remote_data_folders[data_kind]) / relative_path).as_posix()


def get_file_type(file_name):
    """Return the type assigned to *file_name* from its known suffix."""
    for suffix, file_type in sorted(file_suffix_types.items(), key=lambda item: -len(item[0])):
        if file_name.endswith(suffix):
            return file_type
    return "unclassified"


def get_recursive_source_image(file_path, position_root, annotation):
    """Return a best-effort source image and reason for a nested file."""
    if file_path.suffix.lower() in {".tif", ".tiff", ".czi"}:
        return file_path, ""

    position_path = next(
        (
            parent
            for parent in file_path.parents
            if parent.parent == position_root
            and re.fullmatch(r"Position_(\d+)", parent.name)
        ),
        None,
    )
    if position_path is None:
        return None, "no Position_<number> ancestor"

    images_path = position_path / "Images"
    if not images_path.is_dir():
        return None, "position has no Images folder"

    source_candidates = sorted(
        path
        for path in images_path.iterdir()
        if path.is_file() and path.name.endswith(annotation["source_image"])
    )
    if len(source_candidates) == 1:
        return source_candidates[0], ""
    if not source_candidates:
        return None, "no matching default source image"
    return None, "multiple matching default source images"


def create_file_table(
    data_root=root_folder, 
    raw_microscopy_path=root_raw_microscopy_path,
    file_suffix_counts=None, 
    unclassified_files=None,
    unresolved_annotations=None):
    """Return annotated files under *data_root* with validated source images."""
    data_root = Path(data_root)
    rows = []
    custom_suffixes = tuple(sorted(custom_source_image, key=len, reverse=True))

    for experiment_folder, annotation in file_annotation_dict.items():
        experiment_path = data_root / experiment_folder
        if not experiment_path.is_dir():
            warnings.warn(f"Experiment folder does not exist: {experiment_path}")
            continue

        included_files = set()
        position_root = experiment_path
        if "_PosFolder" in annotation:
            position_root /= annotation["_PosFolder"].lstrip("\\\\/")

        if not position_root.is_dir():
            warnings.warn(f"Position folder does not exist: {position_root}")
            continue

        for position_path in sorted(position_root.glob("Position_*")):
            if not position_path.is_dir():
                continue

            position_match = re.fullmatch(r"Position_(\d+)", position_path.name)
            if position_match is None:
                warnings.warn(f"Skipping folder with an unrecognised position: {position_path}")
                continue

            position = int(position_match.group(1))
            strain = annotation["position_to_strain"].get(position, "")
            if not strain:
                warnings.warn(
                    f"Position {position} is not annotated for experiment "
                    f"{experiment_folder}; files will be marked unresolved."
                )

            images_path = position_path / "Images"
            if not images_path.is_dir():
                warnings.warn(f"Images folder does not exist: {images_path}")
                continue

            image_files = sorted(path for path in images_path.iterdir() if path.is_file())
            shared_part = os.path.commonprefix(
                [path.name for path in image_files]
            ).rstrip("_")
            if file_suffix_counts is not None:
                file_suffix_counts.update(
                    path.name.removeprefix(shared_part) for path in image_files
                )
            if unclassified_files is not None:
                unclassified_files.update(
                    path.relative_to(data_root).as_posix()
                    for path in image_files
                    if get_file_type(path.name) == "unclassified"
                )
            default_sources = [
                path
                for path in image_files
                if path.name.endswith(annotation["source_image"])
            ]
            if len(default_sources) != 1:
                raise FileNotFoundError(
                    f"Expected exactly one source image ending in "
                    f"{annotation['source_image']!r} in {images_path}, found "
                    f"{len(default_sources)}."
                )
            default_source = default_sources[0]

            for file_path in image_files:
                if not file_path.is_file():
                    continue

                source_path = default_source
                for mapped_suffix, mapped_source in annotation.get("_special_mapper", {}).items():
                    if file_path.name.endswith(mapped_suffix):
                        source_path = images_path / (
                            file_path.name.removesuffix(mapped_suffix) + mapped_source
                        )
                        break
                else:
                    for suffix in custom_suffixes:
                        if file_path.name.endswith(suffix):
                            source_stem = file_path.name.removesuffix(suffix)
                            source_candidates = [
                                path
                                for path in image_files
                                if path.stem.lstrip(".") == source_stem.lstrip(".")
                            ]
                            if not source_candidates and "_phc" in source_stem:
                                warnings.warn(
                                    f"_phc found in source stem {source_stem!r}, trying phase contrast alternative."
                                )
                                phase_contrast_stem = source_stem.replace(
                                    "_phc", "_phase_contr"
                                )
                                source_candidates = [
                                    path
                                    for path in image_files
                                    if path.stem.lstrip(".")
                                    == phase_contrast_stem.lstrip(".")
                                ]
                            if not source_candidates and "_phc" in source_stem:
                                source_path = default_source
                                break
                            if len(source_candidates) != 1:
                                raise FileNotFoundError(
                                    f"Expected exactly one custom source for {file_path}, "
                                    f"found {len(source_candidates)} files with stem "
                                    f"{source_stem!r}."
                                )
                            source_path = source_candidates[0]
                            break

                if not source_path.is_file():
                    raise FileNotFoundError(
                        f"Source image for {file_path} does not exist: {source_path}"
                    )

                rows.append(
                    {
                        "Files": get_remote_path(
                            "analysed", file_path.relative_to(data_root)
                        ),
                        "source_image": get_remote_path(
                            "analysed", source_path.relative_to(data_root)
                        ),
                        "Date": annotation["Date"],
                        "condition": annotation["Condition"],
                        "experiment": annotation["experiment"],
                        "strain": strain,
                        "Position": position,
                        "file_type": get_file_type(file_path.name),
                        "annotation_status": "annotated" if strain else "unresolved",
                        "annotation_reason": (
                            "" if strain else f"position {position} is not annotated"
                        ),
                        "_add_aspera_to_file_table": True,
                        "switch": annotation.get("switch", ""),
                        "_aspera_path": file_path.relative_to(data_root).as_posix()
                    }
                )
                included_files.add(file_path)

        for file_path in sorted(experiment_path.rglob("*")):
            relative_parts = file_path.relative_to(experiment_path).parts
            if (
                not file_path.is_file()
                or file_path in included_files
                or any(
                    part.casefold() in excluded_experiment_folder_names
                    for part in relative_parts
                )
            ):
                continue

            position_match = next(
                (
                    re.fullmatch(r"Position_(\d+)", parent.name)
                    for parent in file_path.parents
                    if parent.parent == position_root
                ),
                None,
            )
            position = int(position_match.group(1)) if position_match else ""
            strain = annotation["position_to_strain"].get(position, "")
            source_path, reason = get_recursive_source_image(
                file_path, position_root, annotation
            )
            if position_match is None:
                reason = "no Position_<number> ancestor"
            elif not strain:
                reason = reason or f"position {position} is not annotated"

            status = "annotated" if source_path is not None and not reason else "unresolved"
            relative_path = file_path.relative_to(data_root).as_posix()
            if unresolved_annotations is not None and status == "unresolved":
                unresolved_annotations[f"{reason}: {relative_path}"] += 1

            rows.append(
                {
                    "Files": get_remote_path("analysed", relative_path),
                    "source_image": (
                        get_remote_path(
                            "analysed", source_path.relative_to(data_root)
                        )
                        if source_path is not None
                        else ""
                    ),
                    "Date": annotation["Date"],
                    "condition": annotation["Condition"],
                    "experiment": annotation["experiment"],
                    "strain": strain,
                    "Position": position,
                    "file_type": (
                        "Raw microscopy file"
                        if file_path.suffix.lower() == ".czi"
                        else get_file_type(file_path.name)
                    ),
                    "annotation_status": status,
                    "annotation_reason": reason,
                    "_add_aspera_to_file_table": True,
                    "switch": annotation.get("switch", ""),
                    "_aspera_path": relative_path,
                }
            )
       
    df1 = pd.DataFrame(
        rows,
        columns=[
            "Files",
            "source_image",
            "Date",
            "condition",
            "experiment",
            "strain",
            "Position",
            "file_type",
            "annotation_status",
            "annotation_reason",
            "switch",
            "_add_aspera_to_file_table",
            "_aspera_path",
        ],
    )
    rows_2 = []
    raw_microscopy_path = Path(raw_microscopy_path)
    for file_path in raw_microscopy_path.glob("*"):
        if not file_path.is_file():
            continue

        numeric_parts = re.findall(r"\d+", file_path.name)
        date = max(numeric_parts, key=len, default=None)
        date = custom_raw_to_annotation_mapper.get(date, date)
        annotation = file_annotation_dict.get(date, {})
        if not annotation:
            annotation = file_annotation_dict_by_date.get(date, {})
        if not annotation:
            warnings.warn(f"No annotation found for date {date} in file {file_path}")

        relative_path = file_path.relative_to(raw_microscopy_path).as_posix()
        rows_2.append(
            {
                "Files": get_remote_path("raw", relative_path),
                "Date": annotation.get("Date", ""),
                "condition": annotation.get("Condition", ""),
                "experiment": annotation.get("experiment", ""),
                "strain": annotation.get("position_to_strain", ""),
                "file_type": "Raw microscopy file",
                "annotation_status": "annotated" if annotation else "unresolved",
                "annotation_reason": "" if annotation else "no date annotation",
                "switch": annotation.get("switch", ""),
                "_add_aspera_to_file_table": True,
                "_aspera_path": relative_path,
            }
        )

    df2 = pd.DataFrame(
        rows_2,
        columns=[
            "Files",
            "Date",
            "condition",
            "experiment",
            "strain",
            "file_type",
            "annotation_status",
            "annotation_reason",
            "switch",
            "_add_aspera_to_file_table",
            "_aspera_path",
        ],
    )
    return df1, df2


def write_aspera_file_list(file_table, output_path):
    """Write the relative source paths that Aspera should upload."""
    aspera_list = file_table[
        file_table["_add_aspera_to_file_table"]
        ]["_aspera_path"].drop_duplicates()
    aspera_list.to_csv(
        output_path,
        index=False,
        header=False,
        lineterminator="\n",
    )
    return aspera_list


def create_analysis_table_file_table(analysis_folder="test_dfs_for_final_code"):
    """Return a manifest of all analysis tables and their source images."""
    rows = []
    for table_path, source_images in analysis_table_source_images.items():
        source_relative_path = f"{analysis_folder}/{table_path}"
        if source_images is None:
            source_images = []
        elif isinstance(source_images, str):
            source_images = [source_images]

        rows.append(
            {
                "Files": get_remote_path("analysis_tables", source_relative_path),
                "source_images": "; ".join(
                    get_remote_path("analysed", source_image)
                    for source_image in source_images
                ),
                "file_type": "derived analysis table",
                "_add_aspera_to_file_table": True,
                "_aspera_path": source_relative_path,
            }
        )

    return pd.DataFrame(
        rows,
        columns=[
            "Files",
            "source_images",
            "file_type",
            "_add_aspera_to_file_table",
            "_aspera_path",
        ],
    )


def create_fish_file_table(fish_root=root_fish_path):
    """Return all FISH files with strain and source-image annotations."""
    fish_root = Path(fish_root)
    source_images_by_position = {}
    strain_names = {
        "MMY": "MMY",
        "YCY001_1": "bck2d",
        "YCY001-1": "bck2d",
        "YCY002_1": "ccr4d",
        "YCY002-1": "ccr4d",
        "YCY004_5": "whi5d",
        "YCY004-5": "whi5d",
    }
    rows = []

    for file_path in fish_root.rglob("*"):
        if not file_path.is_file():
            continue
        relative_path = file_path.relative_to(fish_root).as_posix()

        # skip the BCK2 nutrient switch timecourse FISH fodler
        if "BCK2_nutrient_switch" in relative_path:
            continue

        tiffs_path = next(
            (parent for parent in file_path.parents if parent.name.lower() == "tiffs"),
            None,
        )
        spotmax_path = next(
            (
                parent
                for parent in file_path.parents
                if "spotmax" in parent.name.lower()
            ),
            None,
        )
        strain_folder = next(
            (
                parent
                for parent in file_path.parents
                if parent.name.lower() in {name.lower() for name in strain_names}
            ),
            None,
        )
        strain = strain_names.get(str(strain_folder), "")
        condition = "SCGE" if "SCGE" in relative_path.upper() else "SCD"

        position_path = None
        for parent in file_path.parents:
            images_path = parent / "Images"
            if images_path.is_dir():
                position_path = parent
                break

        if file_path.suffix.lower() == ".czi":
            source_images = []
            file_type = "Raw microscopy file"
        elif spotmax_path is not None:
            spotmax_tiffs_path = spotmax_path.parent / "TIFFs"
            source_images = [
                spotmax_tiffs_path.relative_to(fish_root).as_posix()
            ] if spotmax_tiffs_path.is_dir() else []
            file_type = "spotMAX output"
        elif position_path is not None and position_path not in source_images_by_position:
            images_path = position_path / "Images"
            source_images_by_position[position_path] = sorted(
                image.relative_to(fish_root).as_posix()
                for image in images_path.iterdir()
                if image.is_file() and image.suffix.lower() in {".tif", ".tiff", ".czi"}
            )
        
        if (
            spotmax_path is None
            and position_path is not None
            and file_path.suffix.lower() != ".czi"
        ):
            position_source_images = source_images_by_position[position_path]
            if file_path.suffix.lower() in {".tif", ".tiff", ".czi"}:
                source_images = [relative_path]
                file_type = get_file_type(file_path.name)
            else:
                brightfield_sources = [
                    image
                    for image in position_source_images
                    if image.lower().endswith(("_bf.tif", "_phase_contr.tif", "_dia_ph3.tif"))
                ]
                source_images = brightfield_sources or position_source_images
                file_type = get_file_type(file_path.name)
        elif spotmax_path is None and tiffs_path is not None:
            source_images = [tiffs_path.relative_to(fish_root).as_posix()]
            file_type = get_file_type(file_path.name)
        elif spotmax_path is None:
            summary_path = next(
                (parent.parent for parent in file_path.parents if parent.name == "results_combined"),
                None,
            )
            source_images = [summary_path.relative_to(fish_root).as_posix()] if summary_path else []
            file_type = get_file_type(file_path.name)

        if file_path.suffix.lower() == ".czi":
            source_images = []
            file_type = "Raw microscopy file"

        rows.append(
            {
                "Files": get_remote_path("fish", relative_path),
                "source_images": "; ".join(
                    get_remote_path("fish", source_image)
                    for source_image in source_images
                ),
                "strain": strain,
                "condition": condition,
                "file_type": file_type,
                "_add_aspera_to_file_table": True,
                "_aspera_path": relative_path,
            }
        )

    return pd.DataFrame(
        rows,
        columns=[
            "Files",
            "source_images",
            "strain",
            "condition",
            "file_type",
            "_add_aspera_to_file_table",
            "_aspera_path",
        ],
    )


def create_qpcr_file_table(qpcr_root=root_qpcr_path, unresolved_annotations=None):
    """Return a manifest of qPCR files with their experiment and run folders."""
    qpcr_root = Path(qpcr_root)
    rows = []

    for file_path in sorted(qpcr_root.rglob("*")):
        if not file_path.is_file():
            continue

        relative_path = file_path.relative_to(qpcr_root)
        qpcr_experiment = relative_path.parts[0] if relative_path.parts else ""
        qpcr_run = relative_path.parts[1] if len(relative_path.parts) > 2 else ""
        is_final_paper_summary = (
            qpcr_experiment.casefold() == "final paper data"
            and len(relative_path.parts) == 2
            and file_path.suffix.casefold() == ".xlsx"
        )
        file_type = "summary" if is_final_paper_summary or (
            qpcr_experiment.casefold() != "final paper data"
            and file_path.suffix.casefold() == ".xlsx"
        ) else "qPCR output"
        relative_path_string = relative_path.as_posix()
        if not qpcr_experiment:
            annotation_reason = "no qPCR experiment folder"
        elif len(relative_path.parts) == 1:
            annotation_reason = "no qPCR experiment folder"
        elif not qpcr_run and not is_final_paper_summary:
            annotation_reason = "no qPCR run folder"
        else:
            annotation_reason = ""
        annotation_status = "annotated" if not annotation_reason else "unresolved"
        if unresolved_annotations is not None and annotation_reason:
            unresolved_annotations[
                f"{annotation_reason}: {relative_path_string}"
            ] += 1

        rows.append(
            {
                "Files": get_remote_path("qpcr", relative_path_string),
                "Type": file_type,
                "qPCR_experiment": qpcr_experiment,
                "qPCR_run": qpcr_run,
                "annotation_status": annotation_status,
                "annotation_reason": annotation_reason,
                "_add_aspera_to_file_table": True,
                "_aspera_path": relative_path_string,
            }
        )

    return pd.DataFrame(
        rows,
        columns=[
            "Files",
            "Type",
            "qPCR_experiment",
            "qPCR_run",
            "annotation_status",
            "annotation_reason",
            "_add_aspera_to_file_table",
            "_aspera_path",
        ],
    )


def print_fish_orphaned_files(file_table):
    """Print FISH records without strain, source, or a known file type."""
    orphaned_files = file_table[
        file_table["strain"].eq("")
        & file_table["source_images"].eq("")
        & file_table["file_type"].eq("unclassified")
    ]
    print(f"\nOrphaned FISH files: {len(orphaned_files)}")
    for file_path in orphaned_files["Files"]:
        print(f"  {file_path}")

if __name__ == "__main__":
    file_suffix_counts = Counter()
    unclassified_files = Counter()
    unresolved_annotations = Counter()
    df1, df2 = create_file_table(
        file_suffix_counts=file_suffix_counts,
        unclassified_files=unclassified_files,
        unresolved_annotations=unresolved_annotations,
    )
    output_path_1 = Path(__file__).with_name("file_table.tsv")
    aspera_list_path_1 = Path(__file__).with_name("aspera_file_list.txt")
    aspera_list_1 = write_aspera_file_list(df1, aspera_list_path_1)
    df1 = df1.drop(columns=["_add_aspera_to_file_table", "_aspera_path"])
    df1.to_csv(output_path_1, sep="\t", index=False)
    print(f"Saved {len(df1)} rows to {output_path_1}")
    print(f"Saved {len(df1['Files'].drop_duplicates())} paths to {aspera_list_path_1}")
    
    output_path_2 = Path(__file__).with_name("file_table_raw.tsv")
    aspera_list_path_2 = Path(__file__).with_name("aspera_file_list_raw.txt")
    aspera_list_2 = write_aspera_file_list(df2, aspera_list_path_2)
    df2 = df2.drop(columns=["_add_aspera_to_file_table", "_aspera_path"])
    df2.to_csv(output_path_2, sep="\t", index=False)
    print(f"Saved {len(df2)} rows to {output_path_2}")
    print(f"Saved {len(df2['Files'].drop_duplicates())} paths to {aspera_list_path_2}")

    df3 = create_analysis_table_file_table()
    output_path_3 = Path(__file__).with_name("file_table_dfs.tsv")
    aspera_list_path_3 = Path(__file__).with_name("aspera_file_list_dfs.txt")
    aspera_list_3 = write_aspera_file_list(df3, aspera_list_path_3)
    df3 = df3.drop(columns=["_add_aspera_to_file_table", "_aspera_path"])
    df3.to_csv(output_path_3, sep="\t", index=False)
    print(f"Saved {len(df3)} rows to {output_path_3}")
    print(f"Saved {len(aspera_list_3)} paths to {aspera_list_path_3}")

    df4 = create_fish_file_table()
    print_fish_orphaned_files(df4)
    df4_raw = df4[df4["file_type"] == "Raw microscopy file"]
    df4_data = df4[df4["file_type"] != "Raw microscopy file"]

    output_path_4 = Path(__file__).with_name("file_table_fish.tsv")
    aspera_list_path_4 = Path(__file__).with_name("aspera_file_list_fish.txt")
    aspera_list_4 = write_aspera_file_list(df4_data, aspera_list_path_4)
    df4_data.drop(columns=["_add_aspera_to_file_table", "_aspera_path"]).to_csv(
        output_path_4, sep="\t", index=False
    )
    print(f"Saved {len(df4_data)} rows to {output_path_4}")
    print(f"Saved {len(aspera_list_4)} paths to {aspera_list_path_4}")

    output_path_4_raw = Path(__file__).with_name("file_table_fish_raw.tsv")
    aspera_list_path_4_raw = Path(__file__).with_name("aspera_file_list_fish_raw.txt")
    aspera_list_4_raw = write_aspera_file_list(df4_raw, aspera_list_path_4_raw)
    df4_raw.drop(columns=["_add_aspera_to_file_table", "_aspera_path"]).to_csv(
        output_path_4_raw, sep="\t", index=False
    )
    print(f"Saved {len(df4_raw)} rows to {output_path_4_raw}")
    print(f"Saved {len(aspera_list_4_raw)} paths to {aspera_list_path_4_raw}")

    unresolved_qpcr_annotations = Counter()
    df5 = create_qpcr_file_table(
        unresolved_annotations=unresolved_qpcr_annotations
    )
    output_path_5 = Path(__file__).with_name("file_table_qpcr.tsv")
    aspera_list_path_5 = Path(__file__).with_name("aspera_file_list_qpcr.txt")
    aspera_list_5 = write_aspera_file_list(df5, aspera_list_path_5)
    df5.drop(columns=["_add_aspera_to_file_table", "_aspera_path"]).to_csv(
        output_path_5, sep="\t", index=False
    )
    print(f"Saved {len(df5)} rows to {output_path_5}")
    print(f"Saved {len(aspera_list_5)} paths to {aspera_list_path_5}")
    print("\nUsed file suffixes:")
    
    for suffix, count in sorted(file_suffix_counts.items()):
        print(f"{count:>6}  {suffix}")
    print("\nUnclassified files:")
    if unclassified_files:
        for file_path, count in sorted(unclassified_files.items()):
            print(f"{count:>6}  {file_path}")
    else:
        print("None")
    print("\nUnresolved annotations:")
    if unresolved_annotations:
        for file_path, count in sorted(unresolved_annotations.items()):
            print(f"{count:>6}  {file_path}")
    else:
        print("None")
    print("\nUnresolved qPCR annotations:")
    if unresolved_qpcr_annotations:
        for file_path, count in sorted(unresolved_qpcr_annotations.items()):
            print(f"{count:>6}  {file_path}")
    else:
        print("None")
    print("Length annotation dictionary:", len(file_annotation_dict))
    print("Len Aspera list:", len(aspera_list_1))
    print("Len Aspera list for raw:", len(aspera_list_2))
    print("Len Aspera list for analysis tables:", len(aspera_list_3))
    print("Len Aspera list for fish:", len(aspera_list_4))
    print("Len Aspera list for fish raw:", len(aspera_list_4_raw))
    print("Len Aspera list for qPCR:", len(aspera_list_5))

