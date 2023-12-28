## Step 1: Prepare for the ChIP-seq data
Pre_download_data.R, Pre_TFs_and_md5_rmv_data.R
## Step 2: Prepare for the WGBS
Make_WGBS.R, including chrseparate.py
## Step 3: Call methylation levels of the peaks
TF_methyl_aligment.R
## Step 4: Calculate Overlap ratio and consine, identify the significant interactions
Call_Significant_overlapratio_consine.R
## Step 5: Call USFs and the features of USFs
Stripe_TF.R
## Step 6: Call co-binding TF pairs and modules
Paired_TF.R
## Step 7: Get the DNA sequencece and call raw motifs by Homer
Separete_peaksToHomer.R
## Step 8: Identified the TF bindingsites and rebuild the motif and plot the 5-letter DNA motifs
Rebulid_motif.R, Rebulid_logo.R
