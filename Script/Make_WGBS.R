#########################
# merge wgbs: before merge we separete the file based on chr i by using the 'python chrseparate.py '
###########################
combined.wgbs.file <- function(pathway,filename,pathwayw,filenamew)
{
  WGBS.rep1.filename <- paste0(pathway, filename[1]);
  WGBS.rep2.filename <- paste0(pathway, filename[2]);
  wgbs1 <- read.delim(WGBS.rep1.filename, header=FALSE)
  wgbs2 <- read.delim(WGBS.rep2.filename, header=FALSE)


  wgbs1.methy.read <- round(wgbs1[,10]*wgbs1[,11]/100);
  wgbs2.methy.read <- round(wgbs2[,10]*wgbs2[,11]/100);
  wgbs1 <- cbind(wgbs1, wgbs1.methy.read);
  wgbs2 <- cbind(wgbs2, wgbs2.methy.read);

  tempchr <- as.character(unique(wgbs1[,1]));#get chr 1:22,X Y

  position <- unique(c(as.character(wgbs1[,2]), as.character(wgbs2[,2])));
  strand <- array("", dim=length(position));
  strand[match(wgbs1[,2], position)] <- as.character(wgbs1[,6]);
  strand[match(wgbs2[,2], position)] <- as.character(wgbs2[,6]);
  all.read <- array(0, dim=length(position));
  methy.read <- array(0, dim=length(position));
  all.read[match(wgbs1[,2], position)] <- all.read[match(wgbs1[,2], position)] + as.numeric(wgbs1[,10]);
  all.read[match(wgbs2[,2], position)] <- all.read[match(wgbs2[,2], position)] + as.numeric(wgbs2[,10]);
  methy.read[match(wgbs1[,2], position)] <- methy.read[match(wgbs1[,2], position)] + as.numeric(wgbs1[,12]);
  methy.read[match(wgbs2[,2], position)] <- methy.read[match(wgbs2[,2], position)] + as.numeric(wgbs2[,12]);
  chr<- array(tempchr, dim=length(position));
  combined.wgbs <- cbind(chr, position, strand, all.read, methy.read);


  methyl.pos <- combined.wgbs [combined.wgbs [,3]=="+", ];
  methyl.neg <- combined.wgbs [combined.wgbs [,3]=="-", ];
  methyl.neg <- data.frame(methyl.neg[,1], as.numeric(methyl.neg[,2])-1, methyl.neg[,3:5]);
  chr.position <- unique(c(as.character(methyl.pos[,2]), as.character(methyl.neg[,2])));
  all.read <- array(0, dim=length(chr.position));
  methy.read <- array(0, dim=length(chr.position));
  all.read[match(methyl.pos[,2], chr.position)] <- all.read[match(methyl.pos[,2], chr.position)] + as.numeric(as.character(methyl.pos[,4]));
  all.read[match(methyl.neg[,2], chr.position)] <- all.read[match(methyl.neg[,2], chr.position)] + as.numeric(as.character(methyl.neg[,4]));
  methy.read[match(methyl.pos[,2], chr.position)] <- methy.read[match(methyl.pos[,2], chr.position)] + as.numeric(as.character(methyl.pos[,5]));
  methy.read[match(methyl.neg[,2], chr.position)] <- methy.read[match(methyl.neg[,2], chr.position)] + as.numeric(as.character(methyl.neg[,5]));
  chr<- array(tempchr, dim=length(chr.position));
  combined.wgbs<- cbind(chr, chr.position, all.read, methy.read);
  cat(paste0("End ",as.character(filenamew),"\n"));
  output.filename <- paste0(pathwayw,filenamew);
  write.csv(combined.wgbs, file=output.filename, row.names=FALSE, quote=FALSE)
}
plot.WGBS.methylation.levels <- function(pathway,filename,pathwayw)
{
  # plot all sites figure
  wgbs.filename <- paste0(pathway,filename);
  png.filename <- paste0(pathwayw, filename,".png");
  wgbs <- read.csv(wgbs.filename);
  methyl.per <- wgbs[,"methy.read"]/wgbs[,"all.read"];
  png(file=paste0(pathwayw,filename,".png"), width = 960, height = 480);
  par(mfrow=c(1,2))
  hist(methyl.per, nclass=20,
       main=paste0("WGBS distribution within ",  filename),
       xlab="methylation ratio of all sites")
  hist(methyl.per[methyl.per>=0.1], nclass=20,
       main=paste0("WGBS distribution within ",  filename),
       xlab="methylation ratio>=0.1 of all sites")
  dev.off()


  methyl.per <- methyl.per[which(wgbs[,"all.read"]>5)]
  png(file=paste0(pathwayw,filename,"(>5 reads).png"),width = 960, height = 480);
  par(mfrow=c(1,2))
  hist(methyl.per, nclass=20,
       main=paste0("WGBS distribution within ",  filename),
       xlab="methylation ratio of >5 reads sites")
  hist(methyl.per[methyl.per>=0.1], nclass=20,
       main=paste0("WGBS distribution within ",  filename),
       xlab="methylation ratio>=0.1 of >5 reads sites")
  dev.off()


}
calling_merge_WGBS <- function(cell,fileID1,fileID2)
{
  system(paste0("mkdir data/WGBS/",cell,"_Merge_GRCh38/"))
  ALLWGBS <- NULL
  for(i in c(1:22,"X","Y"))
  {
    print(paste0("Start to merge chr",i,".out.bed"))
    pathway <-paste0("data/WGBS/",cell,"_Separete_GRCh38/")
    filename <- c(paste0(fileID1,"chr",i,"out.bed"),paste0(fileID2,"chr",i,"out.bed"))
    pathwayw <- paste0("data/WGBS/",cell,"_Merge_GRCh38/")
    filenamew <- c(paste0("chr",i,".bed"))
    combined.wgbs.file(pathway,filename,pathwayw,filenamew)
    plot.WGBS.methylation.levels(pathwayw,filenamew,pathwayw)
    temp <- read.csv(paste0("data/WGBS/",cell,"_Merge_GRCh38/chr",i,".bed"))
    ALLWGBS  <-rbind(ALLWGBS ,temp)
  }
  colnames(ALLWGBS) <- c("chr","position","all.read","methy.read")
  write.csv(ALLWGBS,file = paste0("data/WGBS/",cell,".GRCh38.WGBS.csv"),quote = F,row.names = F)
}
#calling_merge_WGBS("SK-N-SH","ENCFF179VKR","ENCFF940XWW")
#calling_merge_WGBS("A549","ENCFF005TID","ENCFF003JVR")

temp_merge_MCF7 <- function(M_file,cove_file)
{
  MCF7xten_cov <- read.delim("data/WGBS/GSM3336908_MCF7xten.cov.bedGraph", header=FALSE)
  MCF7xten_cov[which(MCF7xten_cov[,4]==-1),4]<-0
  MCF7xten <- read.delim("data/WGBS/GSM3336908_MCF7xten.bedGraph", header=FALSE)
  MCF7xten[which(MCF7xten[,4]==-1),4]<-0
  MCF7xten_cov_index <- paste0(MCF7xten_cov[,1],"_",MCF7xten_cov[,2])
  MCF7xten_index <- paste0(MCF7xten[,1],"_",MCF7xten[,2])
  match_index <- match(MCF7xten_cov_index,MCF7xten_index)
  ALLWGBS <- cbind(MCF7xten_cov[,c(1,2,4)],MCF7xten[match_index,4])
  ALLWGBS[is.na(ALLWGBS[,4]),4] <- 0
  ALLWGBS[,4] <- round(ALLWGBS[,3]*ALLWGBS[,4])
  colnames(ALLWGBS) <- c("chr","position","all.read","methy.read")
  write.csv(ALLWGBS,file = paste0("data/WGBS/MCF7.GRCh38.WGBS.csv"),quote = F,row.names = F)
}

temp_make_bed_to_Hek293 <- function(M_file,cove_file)
{
  Hek293_cov <- read.delim("data/WGBS/GSE207255_WGBS_CpG_hek293.csv", sep = ",",header=TRUE,stringsAsFactors = FALSE)
  Hek293_cov <- Hek293_cov[,c(2,3,4,7,8)]#chr start   end hek293_wt_nTotal hek293_wt_nMeth
  write.table(Hek293_cov ,quote = F,sep = "\t",file ="data/WGBS/GSE207255_WGBS_CpG_hek293.bed",col.names = F,row.names = F)
  #!!!!!!!!!!!then litfOver##########
  Hek293_cov <- read.delim("data/WGBS/GSE207255_WGBS_CpG_hek293_hg38.bed",header=TRUE,stringsAsFactors = FALSE)
  
  Hek293_cov[,2] <- Hek293_cov[,2]-1
  Hek293_cov <- Hek293_cov[,c(1,2,4,5)]
  colnames( Hek293_cov ) <- c("chr","position","all.read","methy.read")
  write.csv( Hek293_cov ,file = paste0("data/WGBS/HEK293.GRCh38.WGBS.csv"),quote = F,row.names = F)
}