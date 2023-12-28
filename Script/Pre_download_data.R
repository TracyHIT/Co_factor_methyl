#first we use xargs -L 1 curl -O -J -L < files.txt to run the first line to get the metadata.tsv
make_download_meta <- function(Cell="GM12878")
{
  metadata <- read.delim(paste0("data/TF/meta/",Cell,"_metadata.tsv"), stringsAsFactors=FALSE)
  Experiment.accession <- unique(metadata[,"Experiment.accession"])
  TF_all <- unique(metadata[,"Experiment.target"])###here all the TFs

  metadata_GRCh38 <- metadata[metadata[,"File.assembly"]=="GRCh38",]
  # GRCh38_IDR thresholded peaks part1
  metadata_GRCh38_IDR <-metadata_GRCh38[metadata_GRCh38[,"File.format"] =="bed narrowPeak" &
                                                   metadata_GRCh38[,"Output.type"] =="IDR thresholded peaks",]
  GRCh38_TF <- unique(metadata_GRCh38_IDR[,"Experiment.target"])

 # GRCh38_conservative_IDR thresholded peaks part2
  metadata_GRCh38_conservative_IDR <-metadata_GRCh38[metadata_GRCh38[,"File.format"] =="bed narrowPeak" &
                                                       !metadata_GRCh38[,"Experiment.target"] %in% GRCh38_TF &
                                                       metadata_GRCh38[,"Output.type"] =="conservative IDR thresholded peaks",]
  GRCh38_TF <- unique(c(GRCh38_TF,metadata_GRCh38_conservative_IDR[,"Experiment.target"]))

  # GRCh38_optimal_IDR thresholded peaks part3
  metadata_GRCh38_optimal_IDR <-metadata_GRCh38[metadata_GRCh38[,"File.format"] =="bed narrowPeak" &
                                                       !metadata_GRCh38[,"Experiment.target"] %in% GRCh38_TF &
                                                       metadata_GRCh38[,"Output.type"] =="optimal IDR thresholded peaks",]
  GRCh38_TF <- unique(c(GRCh38_TF,metadata_GRCh38_optimal_IDR[,"Experiment.target"]))


  # hg38 less no replicate and add hg19
  # hg19_optimal_IDR thresholded peaks part4
  metadata_hg19 <- metadata[metadata[,"File.assembly"]=="hg19"&
                          !metadata[,"Experiment.target"] %in% GRCh38_TF, ]
  Experiment.accession_hg19 <- unique(metadata_hg19[,"Experiment.target"])
  metadata_hg19_IDR <-metadata_hg19[metadata_hg19[,"File.format"] =="bed narrowPeak" &
                                                     metadata_hg19[,"Output.type"] =="optimal IDR thresholded peaks",]

  #summary omit TF
  include_TF <- c(unique(metadata_GRCh38_IDR[,"Experiment.target"]),
                  unique(metadata_GRCh38_conservative_IDR[,"Experiment.target"]),
                  unique(metadata_GRCh38_optimal_IDR[,"Experiment.target"]),
                 unique(metadata_hg19_IDR[,"Experiment.target"]))
  unused_Exp_TF <- metadata[!metadata[,"Experiment.target"] %in% include_TF,]

  write.csv(rbind(metadata_GRCh38_IDR,metadata_GRCh38_conservative_IDR,metadata_GRCh38_optimal_IDR,metadata_hg19_IDR),file = paste0("data/TF/meta/",Cell,"_metadata_peak.csv"))
  write.csv( unused_Exp_TF, file = paste0("data/TF/meta/",Cell,"_unused_Exp_TF.csv"))

  hg19_download <- paste0("https://www.encodeproject.org/files/", metadata_hg19_IDR[,"File.accession"],"/@@download/",metadata_hg19_IDR[,"File.accession"],".bed.gz")
  GRCH38_download <- c(paste0("https://www.encodeproject.org/files/", metadata_GRCh38_IDR[,"File.accession"],"/@@download/",metadata_GRCh38_IDR[,"File.accession"],".bed.gz"),
                       paste0("https://www.encodeproject.org/files/", metadata_GRCh38_conservative_IDR[,"File.accession"],"/@@download/",metadata_GRCh38_conservative_IDR[,"File.accession"],".bed.gz"),
                       paste0("https://www.encodeproject.org/files/", metadata_GRCh38_optimal_IDR[,"File.accession"],"/@@download/",metadata_GRCh38_optimal_IDR[,"File.accession"],".bed.gz"))

  write.table(hg19_download , file = paste0("data/TF/meta/",Cell,"_hg19_download.files"),col.names = F,row.names = F,quote = F)
  write.table(GRCH38_download , file = paste0("data/TF/meta/",Cell,"_GRCH38_download.files"),col.names = F,row.names = F,quote = F)

  system(paste0("mkdir data/TF/",Cell))
  system(paste0("mkdir data/TF/",Cell))
  system(paste0("mkdir data/TF/",Cell,"/Raw_hg19"))
  system(paste0("mkdir data/TF/",Cell,"/Raw_hg38"))

}
####call the function#####
celllist = c("HepG2","K562","GM12878","MCF-7","HEK293","H1","A549","HeLa-S3","SK-N-SH","WTC11")
for(cell in celllist)
{
  print(paste("Start", cell,"!!!!!!!!!!!!!!!!!!!!!"))
  make_download_meta(cell)
  system(paste0("xargs -L 1 curl -O -J -L < ",paste0("data/TF/meta/",cell,"_hg19_download.files")))
  system(paste0("xargs -L 1 curl -O -J -L < ",paste0("data/TF/meta/",cell,"_GRCH38_download.files")))
  print(paste("End", cell,"!!!!!!!!!!!!!!!!!!!!!"))
}



