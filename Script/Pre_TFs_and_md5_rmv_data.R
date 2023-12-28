#it is a temp function to move the files, it may be included in call_check md5
temp_mv_files <- function( cell = "GM12878")
{
  cell_metadata <- read.csv(paste0("data/TF/meta/",cell,"_metadata_peak.csv"), stringsAsFactors=FALSE)[,c(2,3,4,5,6,7,8,12,24,45,47)]
  for(i in 1:nrow(cell_metadata))
  {
    switch(cell_metadata[i,"File.assembly"],
             "GRCh38"=system(paste0("mv ",cell_metadata[i,"File.accession"],".bed.gz data/TF/", cell,"/Raw_hg38/",cell_metadata[i,"File.accession"],".bed.gz")),
             "hg19"=system(paste0("mv ",cell_metadata[i,"File.accession"],".bed.gz data/TF/", cell,"/Raw_hg19/",cell_metadata[i,"File.accession"],".bed.gz"))
    )
  }
}

############check md5###################
#run this to confirm the file is downloaded with no error, please run it iterationly by hand.
check_md5 <- function( cell = "GM12878")
{
  cell_metadata <- read.csv(paste0("data/TF/meta/",cell,"_metadata_peak.csv"), stringsAsFactors=FALSE)[,c(2,3,4,5,6,7,8,12,24,45,47)]
  if(length(grep("POLR2A", cell_metadata[,"Experiment.target"])) >0)
  {
    cell_metadata <- cell_metadata[-grep("POLR2A", cell_metadata[,"Experiment.target"]),]
  }
  Check_md5sum <- rep("",nrow(cell_metadata))
  for(i in 1: nrow(cell_metadata))
  {
    if(file.exists(paste0("data/TF/", cell,"/Raw_hg38/",cell_metadata[i,"File.accession"],".bed"))|
       file.exists(paste0("data/TF/", cell,"/Raw_hg19/",cell_metadata[i,"File.accession"],".bed")))
    {
      Check_md5sum[i] <- cell_metadata[i,"md5sum"]
      next()
    }
    if(file.exists(paste0("data/TF/", cell,"/Raw_hg38/",cell_metadata[i,"File.accession"],".bed.gz"))|
       file.exists(paste0("data/TF/", cell,"/Raw_hg19/",cell_metadata[i,"File.accession"],".bed.gz")))
    {
      switch(cell_metadata[i,"File.assembly"],
           "GRCh38"=system(paste0("md5sum data/TF/", cell,"/Raw_hg38/",cell_metadata[i,"File.accession"],".bed.gz > temp.txt")),
           "hg19"=system(paste0("md5sum data/TF/", cell,"/Raw_hg19/",cell_metadata[i,"File.accession"],".bed.gz > temp.txt")))
      temp <- read.table("temp.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)
      system("rm temp.txt")
      if(temp[1,1] == cell_metadata[i,"md5sum"])
      {
        Check_md5sum[i] <-temp[1,1]
      }else
      {
        Check_md5sum[i] <- "Wrong md5"
      }
    }else
    {
      Check_md5sum[i] <- "No_file"
    }
  }
  return(cbind(cell_metadata,Check_md5sum))
}
call_check_md5 <- function()
{
  all_lost_file <- NULL
  celllist = c("HepG2","K562","GM12878","MCF-7","HEK293","H1","A549","HeLa-S3","SK-N-SH","WTC11")
  for(cell in celllist )
  {
    #temp_mv_files(cell)
    check_list <- check_md5(cell)
    temp_lost <- check_list[which(check_list[,12]=="No_file" |check_list[,12]=="Wrong md5"),]
    all_lost_file <-rbind(all_lost_file,temp_lost)
  }
   write.csv( all_lost_file,file = "data/TF/meta/Download_error_files.csv")
}
#call_check_md5()
###########end check md5 ##############

##########overlift hg19Tohg38 and make annotation#############
overLift_annotation <- function()
{
  annotation <- NULL
  celllist = c("GM12878","HepG2","K562","MCF-7","HEK293","H1","A549","HeLa-S3","SK-N-SH","WTC11")
  for(cell in celllist)
  {
    print(cell)
  check_list <- read.csv(paste0("data/TF/meta/",cell,"_metadata_peak.csv"), stringsAsFactors=FALSE)[,c(2,3,4,5,6,7,8,12,24,45,47)]
  if(length(grep("POLR2A",check_list[,"Experiment.target"])) >0)
  {
    check_list <- check_list[-grep("POLR2A", check_list[,"Experiment.target"]),]
  }

  check_list[,"Biosample.term.name"] <- cell
  temp_hg19 <- check_list[check_list[,"File.assembly"] == "hg19",]
  temp_GRCh38 <- check_list[check_list[,"File.assembly"] == "GRCh38",]
  system(paste0("gzip -d data/TF/",cell,"/Raw_hg19/*"))
  system(paste0("gzip -d data/TF/",cell,"/Raw_hg38/*"))
  peaks_count_raw <- NULL
  peaks_count_run <- NULL
  for(i in 1:nrow(temp_GRCh38))
  {
    TF_files <- read.delim(paste0("data/TF/",cell,"/Raw_hg38/",temp_GRCh38[i,"File.accession"],".bed"), header=FALSE)
    peaks_count_raw <- c(peaks_count_raw,nrow(TF_files))
    peaks_count_run <- c(peaks_count_run,nrow(TF_files))
  }
  if(nrow(temp_hg19) > 0)
  {
    system(paste0("mkdir data/TF/",cell,"/Raw_hg19Tohg38"))
    for(j in 1:nrow(temp_hg19))
    {
      #system(paste0("liftOver/liftOver data/TF/",cell,"/Raw_hg19/",temp_hg19[j,"File.accession"],".bed liftOver/hg19ToHg38.over.chain data/TF/",cell,"/Raw_hg19Tohg38/",temp_hg19[j,"File.accession"],".bed unmaped.bed -bedPlus=4"))
      TF_files_raw <- read.delim(paste0("data/TF/",cell,"/Raw_hg19/",temp_hg19[j,"File.accession"],".bed"), header=FALSE)
      #TF_files_run <- read.delim(paste0("data/TF/",cell,"/Raw_hg19Tohg38/",temp_hg19[j,"File.accession"],".bed"), header=FALSE)
      peaks_count_raw <- c(peaks_count_raw,nrow(TF_files_raw))
      peaks_count_run <- c(peaks_count_run,0)#nrow(TF_files_run))
    }

  }
  check_list <- cbind(check_list,peaks_count_raw,peaks_count_run)
  annotation <- rbind(annotation,check_list)
}
  write.csv(annotation,file = "data/annotation/annotation.raw.csv",quote = F,row.names = F)
  annotation_fit <- annotation[which(annotation[,"peaks_count_run"] >= 500),]
  write.csv(annotation_fit,file = "data/annotation/annotation.fit.csv",quote = F,row.names = F)
}
overLift_annotation()
#########end overlift hg19Tohg38 and make annotation##########
