## Step 1:Prepare for the ChIP-seq data
using Pre_download_data.R, Pre_TFs_and_md5_rmv_data.R
## Step 2:Prepare for the WGBS
using Make_WGBS.R, including chrseparate.py
## Step 3:call methylation levels of the peaks
TF_methyl_aligment.R
## Step 4:calculate Overlap ratio and consine, identify the significant interactions
Call_Significant_overlapratio_consine.R
