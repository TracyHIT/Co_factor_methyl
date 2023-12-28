#get the DNA sequence and generate the shell script to call raw motif
generate.high.fasta_bed.file <- function()
{
  if(!("package:Biostrings" %in%  search()))
  {
    library(Biostrings);
  }
  files_list <- list.files("Result/TFBS.Peak.methyl/")
  low_methyl_peaks_count <- c()
  high_methyl_peaks_count <- c()
  system_lines <-c()
  for(i in 1:length(files_list))
  {
    print(i)
    tfbs.seq.filename <- files_list[i]
    Peaks.seq <- read.csv(paste0("Result/TFBS.Peak.methyl/",tfbs.seq.filename),stringsAsFactors = F);
    Peaks.index <- paste0(Peaks.seq[,1],"_", Peaks.seq[,2],"_",Peaks.seq[,3])
    tfbs.seq <- read.csv(paste0("Result/TFBS.Seq/",tfbs.seq.filename),stringsAsFactors = F);
    tfbs.seq.index <- paste0(tfbs.seq[,1],"_", tfbs.seq[,2],"_",tfbs.seq[,3])
    Peaks.seq <- cbind(Peaks.seq,tfbs.seq[match(Peaks.index, tfbs.seq.index),"Seq"])
    rm(tfbs.seq.index)

    high.methyl.Peak <- Peaks.seq[Peaks.seq[,"avg.na.methyl.ratio"] >= 0.6 & !is.na(Peaks.seq[,"avg.na.methyl.ratio"]), ];
    low.methyl.Peak <- Peaks.seq[Peaks.seq[,"avg.na.methyl.ratio"] < 0.6 | is.na(Peaks.seq[,"avg.na.methyl.ratio"]), ];
    low_methyl_peaks_count <- c(low_methyl_peaks_count,nrow(low.methyl.Peak))
    high_methyl_peaks_count <- c(high_methyl_peaks_count,nrow(high.methyl.Peak))
    
    name_index <- substr(tfbs.seq.filename,1,nchar(tfbs.seq.filename)-4)
    system(paste0("mkdir Result/Homer/High_0_6/",name_index))
    write.table(high.methyl.Peak,paste0("Result/Homer/High_0_6/",name_index,"/",name_index,".bed"),quote = F,sep = "\t",row.names = F,col.names = F)
    
    if(nrow(high.methyl.Peak)> 100)
    {
      high.methyl.Peak <-  high.methyl.Peak[order(as.numeric(as.character(high.methyl.Peak[,"Signalvalue"])),decreasing = T),]
      Seq <- high.methyl.Peak[1:min(500,nrow(high.methyl.Peak)),22]
      dna = DNAStringSet(Seq);
      writeXStringSet(dna, filepath=paste0("Result/Homer/High_0_6/",name_index,"/",name_index,".fasta"));
      system_lines <-c( system_lines,paste0("findMotifs.pl Result/Homer/High_0_6/",name_index,"/",name_index,".fasta fasta Result/Homer/High_0_6/",name_index," -noknown -S 2")) 
    }
    if(length(system_lines) %% 100 == 0)
    {
      write.table(system_lines[(length(system_lines)-99):length(system_lines)],file = paste0("Result/Homer/Homer_command_lines_",length(system_lines),".sh"),quote = F,row.names = F,col.names = F)
    }
    
    system(paste0("mkdir Result/Homer/Low_0_6/",name_index))
    write.table(low.methyl.Peak,paste0("Result/Homer/Low_0_6/",name_index,"/",name_index,".bed"),quote = F,sep = "\t",row.names = F,col.names = F)
    #system(paste0("findMotifsGenome.pl Result/Homer/High_0_6/",name_index,"/",name_index,".bed hg38 Result/Homer/High_0_6/",name_index,"/ -size 200")) 
    if(nrow(low.methyl.Peak) > 100)
    {
      low.methyl.Peak <-  low.methyl.Peak[order(as.numeric(as.character(low.methyl.Peak[,"Signalvalue"])),decreasing = T),]
      Seq <- low.methyl.Peak[1:min(500,nrow(low.methyl.Peak)),22]
      dna = DNAStringSet(Seq);
      writeXStringSet(dna, filepath=paste0("Result/Homer/Low_0_6/",name_index,"/",name_index,".fasta"));
      system_lines <-c( system_lines,paste0("findMotifs.pl Result/Homer/Low_0_6/",name_index,"/",name_index,".fasta fasta Result/Homer/Low_0_6/",name_index," -noknown -S 2"))
    }
    if(length(system_lines) %% 100 == 0)
    {
      write.table(system_lines[(length(system_lines)-99):length(system_lines)],file = paste0("Result/Homer/Homer_command_lines_",length(system_lines),".sh"),quote = F,row.names = F,col.names = F)
    }
  }
  summary_mehtyl_peaks <- cbind(files_list,low_methyl_peaks_count,high_methyl_peaks_count)
  write.csv(summary_mehtyl_peaks,file = "Result/Homer/summary_mehtyl_peaks.csv",quote = F,row.names = F)
  write.table(system_lines,file = "Result/Homer/Homer_command_lines.sh",quote = F,row.names = F,col.names = F)
}
# generate.high.fasta_bed.file()

#select some scripts just for some TFs motif calling
get_involved_Homer_shell <- function(summary_file="Result/Summary/Strip_TF_all_involved_fileID.csv")
{
  TF_involved <- read.csv(file = summary_file)
  data <- NULL
  for(i in 1:43)
  {
    all_shell_file <- paste0("Result/Homer/Homer_command_lines_",i,"00.sh")
    temp <- read.delim(all_shell_file,header = F)
    data <- rbind(data,temp)
  }
  
  matching_rows <- NULL
  for(i in 1:nrow(TF_involved))
  {
    keyword <- TF_involved[i,"File.accession"]
    matching_rows_temp <- data[grep(keyword, data$V1), ]
    matching_rows <- c(matching_rows,matching_rows_temp)
  }
  write.table(matching_rows, file = "Result/Homer/Homer_command_lines_involved.sh", sep = "\t", row.names = FALSE,quote = F)
}
# summary_file="Result/Summary/Strip_TF_all_involved_fileID.csv"
# get_involved_Homer_shell(summary_file)

#get the DNA sequence overalap peaks and generate the shell script to call raw motif
generate.highOrLow.cobind.fasta_bed.file <- function(TF_1,TF_2,cell,lowORhigh_index="High")
{
  #tfbs.seq.filename_1 <- "GM12878_EGR1-human_ENCFF637UJN_GRCh38.csv"
  #tfbs.seq.filename_2 <- "GM12878_IKZF1-human_ENCFF819VMH_GRCh38.csv"
  if(!("package:Biostrings" %in%  search()))
  {
    library(Biostrings);
  }
  library(GenomicRanges)
  files_list <- list.files("Result/TFBS.Peak.methyl/")
  Peaks.seq_1 <- NULL
  for(i in files_list[grepl(TF_1,files_list)])
  {
    Peaks.seq <- read.csv(paste0("Result/TFBS.Peak.methyl/",i),stringsAsFactors = F);
    Peaks.index <- paste0(Peaks.seq[,1],"_", Peaks.seq[,2],"_",Peaks.seq[,3])
    tfbs.seq <- read.csv(paste0("Result/TFBS.Seq/",i),stringsAsFactors = F);
    tfbs.seq.index <- paste0(tfbs.seq[,1],"_", tfbs.seq[,2],"_",tfbs.seq[,3])
    Peaks.seq_1 <- rbind(Peaks.seq_1,cbind(Peaks.seq,tfbs.seq[match(Peaks.index, tfbs.seq.index),"Seq"]))
  }
  rm(tfbs.seq.index)
    
  Peaks.seq_2<- NULL
  for(i in files_list[grepl(TF_2,files_list)])
  {
    Peaks.seq <- read.csv(paste0("Result/TFBS.Peak.methyl/",i),stringsAsFactors = F);
    Peaks.index <- paste0(Peaks.seq[,1],"_", Peaks.seq[,2],"_",Peaks.seq[,3])
    tfbs.seq <- read.csv(paste0("Result/TFBS.Seq/",i),stringsAsFactors = F);
    tfbs.seq.index <- paste0(tfbs.seq[,1],"_", tfbs.seq[,2],"_",tfbs.seq[,3])
    Peaks.seq_2 <- rbind(Peaks.seq_2,cbind(Peaks.seq,tfbs.seq[match(Peaks.index, tfbs.seq.index),"Seq"]))
  }
  rm(tfbs.seq.index)
    

  if(lowORhigh_index=="High")
  {
  high.methyl.Peak_1 <- Peaks.seq_1[Peaks.seq_1[,"avg.na.methyl.ratio"] >= 0.6 & !is.na(Peaks.seq_1[,"avg.na.methyl.ratio"]), ];
  high.methyl.Peak_2 <- Peaks.seq_2[Peaks.seq_2[,"avg.na.methyl.ratio"] >= 0.6 & !is.na(Peaks.seq_2[,"avg.na.methyl.ratio"]), ];
  }else
  {
    high.methyl.Peak_1 <- Peaks.seq_1[Peaks.seq_1[,"avg.na.methyl.ratio"] < 0.6 | is.na(Peaks.seq_1[,"avg.na.methyl.ratio"]), ];
    high.methyl.Peak_2 <- Peaks.seq_2[Peaks.seq_2[,"avg.na.methyl.ratio"] < 0.6 | is.na(Peaks.seq_2[,"avg.na.methyl.ratio"]), ];
  }
  high.methyl.Peak_1 <- high.methyl.Peak_1[order(as.numeric(high.methyl.Peak_1[,"Signalvalue"]),decreasing = T),]
  high.methyl.Peak_2 <- high.methyl.Peak_2[order(as.numeric(high.methyl.Peak_2[,"Signalvalue"]),decreasing = T),]
  regions_1 <- GRanges(seqnames =  high.methyl.Peak_1[,1],IRanges(start =  high.methyl.Peak_1[,2],end = high.methyl.Peak_1[,3]))
  regions_2 <- GRanges(seqnames =  high.methyl.Peak_2[,1],IRanges(start =  high.methyl.Peak_2[,2],end = high.methyl.Peak_2[,3]))
  Overlap_index <- findOverlaps(regions_1,regions_2)
  temp <- rbind(high.methyl.Peak_1[unique(Overlap_index@from)[1:(min(length(unique(Overlap_index@from)),500))],],
                high.methyl.Peak_2[unique(Overlap_index@to)[1:(min(length(unique(Overlap_index@to)),500))],])
  colnames(temp)[22]<-"Seq"
  Seq<-as.character(temp[,"Seq"])
  dna = DNAStringSet(Seq);
  
  name_index <- paste0(cell,"_",TF_1,"_",TF_2)
  print(name_index)
  system(paste0("mkdir Result/Homer/Co_bind/",name_index,"/"))
  write.table(temp,paste0("Result/Homer/Co_bind/",name_index,"/",name_index,".bed"),quote = F,sep = "\t",row.names = F,col.names = F)
  writeXStringSet(dna, filepath=paste0("Result/Homer/Co_bind/",name_index,"/",name_index,".fasta"));
  system_lines<-paste0("findMotifs.pl Result/Homer/Co_bind/",name_index,"/",name_index,".fasta fasta Result/Homer/Co_bind/",name_index," -noknown -S 2")
}

#call function and make the shell script.
call_generate.highOrLow.cobind.fasta_bed.file <- function()
{
  High_Low <- read.csv(file = "Result/Pair_group/ASummary_High_Low.csv",stringsAsFactors = F)
  High_Low <- High_Low[High_Low[,"cell_count"]>=2,]
  sys_lines <- NULL
  for(i in 1:nrow(High_Low ))
  {
    temp <- High_Low[i,]
    TF_1 <-  unlist(strsplit( temp[,"index"]," _ "))[1]
    TF_2 <-  unlist(strsplit( temp[,"index"]," _ "))[2]
    for(cell in unlist(strsplit( temp[,"cell_list"],"_")))
    {
      sys_line <- generate.highOrLow.cobind.fasta_bed.file(TF_1,TF_2,cell,lowORhigh_index="High")
      sys_lines <- c(sys_lines,sys_line)
    }
  }
  write.table(sys_lines,file = "Result/Homer/Homer_command_cobind_lines_High_Low.sh",quote = F,row.names = F,col.names = F)
  
  High_Low <- read.csv(file = "Result/Pair_group/ASummary_Low_High.csv",stringsAsFactors = F)
  High_Low <- High_Low[High_Low[,"cell_count"]>=2,]
  sys_lines <- NULL
  for(i in 1:nrow(High_Low ))
  {
    temp <- High_Low[i,]
    TF_1 <-  unlist(strsplit( temp[,"index"]," _ "))[1]
    TF_2 <-  unlist(strsplit( temp[,"index"]," _ "))[2]
    for(cell in unlist(strsplit( temp[,"cell_list"],"_")))
    {
      sys_line <- generate.high.cobind.fasta_bed.file(TF_1,TF_2,cell,lowORhigh_index="Low")
      sys_lines <- c(sys_lines,sys_line)
    }
  }
  write.table(sys_lines,file = "Result/Homer/Homer_command_cobind_lines_Low_High.sh",quote = F,row.names = F,col.names = F)
}
