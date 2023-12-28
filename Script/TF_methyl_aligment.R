# get the DNA sequence in the peaks calling from ChIP_seq data
TFBS.SequenceW <- function(TFBS.filename=NULL, output.filename=NULL,dic="/Users/lxmpro/Desktop/hg38/")
{
  library("Biostrings")
  chrindex<-c(as.character(c(1:22)),"X","Y")
  TFBS <- read.table(TFBS.filename, sep="\t")
  TFBS <- TFBS[order(TFBS[,2]), ]
  TFBS <- TFBS[order(TFBS[,1]), ]
  TFBS<-TFBS[which((as.character(TFBS[,1])) %in% paste0("chr",chrindex)),]

  tmp.chr <-""
  Seq <- array("", dim=nrow(TFBS))
  for(i in 1:nrow(TFBS))
  {
    if(tmp.chr != as.character(TFBS[i,1]))
    {
      tmp.chr <- as.character(TFBS[i,1])
      s <- readDNAStringSet(paste0(dic, tmp.chr, '.fa'))
    }
    tryCatch({
      Seq[i] <- toString(subseq(s, start=TFBS[i,2], end=TFBS[i,3]));
    },warning= function(w){
    },error=function(e){
      Seq[i] <- ""
      cat (e, "\n")
    },finally={}
    )
    Seq[i] <- toString(subseq(s, start=TFBS[i,2], end=TFBS[i,3]))
    if (floor(i/1000) == i/1000)
    {
      cat ("Has retrieved ",i, " peaks sequences! \n")
    }
  }
  TFBS <- data.frame(TFBS, Seq)
  colnames(TFBS) <- c("Chrom", "Start", "End", "Name", "Score", "Strand", "Signalvalue", "Pvalue", "Qvalue", "Peak", "Seq")
  TFBS<-TFBS[which(TFBS[,"Seq"]!=""),]
  write.csv(TFBS, file=output.filename, row.names=FALSE)
}
# get the CpG in the peaks 
TFBS.CpG.site <- function(TFBS.filename=NULL, WGBS=NULL, output.filename=NULL)
{
  TFBS <- read.csv(TFBS.filename)
  chrname <- as.character(unique(TFBS[,1]))
  TFBS.CpG.methyl <- NULL
  for (i in 1:length(chrname))
  {
    methyl.chr <-  WGBS[WGBS[,1]==chrname[i], ]
    TFBS.chr <- TFBS[TFBS[,1]==chrname[i], ]
    # retrieve methylation of individual CpG site
    CpG.pos.chr <- NULL;
    TFBS.rep.chr <- NULL;
    for(j in 1:nrow(TFBS.chr))
    {
      CpG.pos <- unlist(gregexpr("CG", as.character(TFBS.chr$Seq[j])))#有多少CG
      if(CpG.pos[1]!=-1)
      {
        CpG.pos <- TFBS.chr$Start[j]+CpG.pos-1
        CpG.pos.chr <- c(CpG.pos.chr, CpG.pos)#竖列位置
        TFBS.rep <- rep(j, each=length(CpG.pos))
        TFBS.rep.chr <- c(TFBS.rep.chr, TFBS.rep)#TFBS的第几个
        if (floor(j/1000) == j/1000)
        {
          cat (j, "\n")
          cat ("Has retrieved methylation sites of ",i, " peaks! \n")
        }
      }
    }
    CpG.all.read <- methyl.chr[match(CpG.pos.chr, methyl.chr[,2]+1), "all.read"]
    CpG.methyl.read <- methyl.chr[match(CpG.pos.chr, methyl.chr[,2]+1), "methy.read"]
    TFBS.CpG.methyl.chr <- cbind(TFBS.chr[TFBS.rep.chr,1:10], CpG.pos.chr, CpG.all.read, CpG.methyl.read)
    TFBS.CpG.methyl <- rbind(TFBS.CpG.methyl, TFBS.CpG.methyl.chr)
    cat(chrname[i], "\tDone!\n");
  }
  colnames(TFBS.CpG.methyl) <- c("Chrom", "Start", "End", "Name", "Score", "Strand", "Signalvalue", "Pvalue", "Qvalue", "Peak", "CpG.pos", "CpG.all.read", "CpG.methyl.read");
  write.csv(TFBS.CpG.methyl, file=output.filename, row.names=FALSE);
}
# get the peaks' methylation 
TFBS.Peak.methylation <- function(TFBS.filename=NULL, TFBS.CpG.filename=NULL, output.filename=NULL)
{
  TFBS <- read.csv(TFBS.filename);
  TFBS.CpG <- read.csv(TFBS.CpG.filename);
  methyl.ratio <- TFBS.CpG$CpG.methyl.read/TFBS.CpG$CpG.all.read;
  TFBS.CpG <- cbind(TFBS.CpG, methyl.ratio);
  CpG.num <- array("", dim=nrow(TFBS));
  methyl.CpG.num <- array("", dim=nrow(TFBS));
  na.CpG.num <- array("", dim=nrow(TFBS));
  max.read.num <- array("", dim=nrow(TFBS));
  max.methyl.num <- array("", dim=nrow(TFBS));
  min.read.num <- array("", dim=nrow(TFBS));
  min.methyl.num <- array("", dim=nrow(TFBS));
  all.read.num <- array("", dim=nrow(TFBS));
  all.methyl.num <- array("", dim=nrow(TFBS));
  avg.nona.methyl.ratio <- array("", dim=nrow(TFBS));
  avg.na.methyl.ratio <- array("", dim=nrow(TFBS));

  tmp.chr <- "";
  for(i in 1:nrow(TFBS))
  {
    if(tmp.chr != as.character(TFBS[i,1]))
    {
      tmp.chr <- as.character(TFBS[i,1]);
      chr.TFBS.CpG <- TFBS.CpG[TFBS.CpG$Chrom==tmp.chr, ];
      temp.TFBS<-TFBS[TFBS$Chrom==tmp.chr,]
      cat(i, "\t", tmp.chr, "\n");
    }
    #chr.TFBS.CpG：取出染色质的所有数据  取出在methyl中有数据的。
    idv.TFBS.CpG.na <- chr.TFBS.CpG[chr.TFBS.CpG$Start==TFBS$Start[i]&chr.TFBS.CpG$End==TFBS$End[i], ];
    CpG.pos <- unlist(gregexpr("CG", as.character(TFBS$Seq[i])));

    if(CpG.pos[1]!=-1)
    {
      CpG.num[i] <- length(CpG.pos);
      methyl.CpG.num[i] <- sum(!is.na(idv.TFBS.CpG.na$CpG.all.read));
      na.CpG.num[i] <- sum(is.na(idv.TFBS.CpG.na$CpG.all.read));
      idv.TFBS.CpG <- idv.TFBS.CpG.na[!is.na(idv.TFBS.CpG.na$CpG.all.read), ];
      if (nrow(idv.TFBS.CpG)>0)
      {
        max.TFBS.CpG <- idv.TFBS.CpG[idv.TFBS.CpG$methyl.ratio==max(idv.TFBS.CpG$methyl.ratio), ];
        max.read.num[i] <- as.numeric(max.TFBS.CpG$CpG.all.read[1]);
        max.methyl.num[i] <- as.numeric(max.TFBS.CpG$CpG.methyl.read[1]);
        min.TFBS.CpG <- idv.TFBS.CpG[idv.TFBS.CpG$methyl.ratio==min(idv.TFBS.CpG$methyl.ratio), ];
        min.read.num[i] <- as.numeric(min.TFBS.CpG$CpG.all.read[1]);
        min.methyl.num[i] <- as.numeric(min.TFBS.CpG$CpG.methyl.read[1]);
        all.read.num[i] <- sum(idv.TFBS.CpG$CpG.all.read);
        all.methyl.num[i] <- sum(idv.TFBS.CpG$CpG.methyl.read);
        avg.nona.methyl.ratio[i] <- mean(idv.TFBS.CpG$methyl.ratio);
        avg.na.methyl.ratio[i] <- sum(idv.TFBS.CpG$methyl.ratio)/nrow(idv.TFBS.CpG.na);
      }
    }
    cat("Has calculated the methylation levels of peaks on ", as.character(TFBS$Chrom[i]),"! \n");
  }
  TFBS.Peak.methyl <- data.frame(TFBS[,1:10], CpG.num, methyl.CpG.num, na.CpG.num, max.read.num, max.methyl.num, min.read.num, min.methyl.num, all.read.num, all.methyl.num, avg.nona.methyl.ratio, avg.na.methyl.ratio);
  write.csv(TFBS.Peak.methyl, file=output.filename, row.names=FALSE);
}


plot.TFBS.CpG.site.methylation.figure <- function(pathway,tfbs.rrbs.filename,startrow,endrow)
{
  TF_meta <-read.csv(tfbs.rrbs.filename);
  if(endrow==(-1))
  {
    endrow=nrow(TF_meta)
  }
  for(i in startrow:endrow)
  {
    TFBS.Peak.filename <- paste0(pathway,"TFBS.Peak.methyl/", as.character(TF_meta[i,"Biosample.term.name"]), ".",as.character(TF_meta[i,"Experiment.target"]),
                                 ".", as.character(TF_meta[i,1]),".csv");
    png.filename <- paste0(pathway,"TFBS.Peak.methyl/","Peak.methyl.figure/",as.character(TF_meta[i,"Biosample.term.name"]),
                           ".",as.character(TF_meta[i,"Experiment.target"]),
                           ".", as.character(TF_meta[i,1]), ".png");
    TFBS.Peak <- read.csv(TFBS.Peak.filename);
    Peak.methyl.ratio <- TFBS.Peak$avg.nona.methyl.ratio;
    Peak.methyl.ratio <- Peak.methyl.ratio[!is.na(Peak.methyl.ratio)];
    png(file=png.filename);
    h <- hist(Peak.methyl.ratio, breaks = 40, plot=FALSE);
    h$counts=h$counts/sum(h$counts);
    main.text <- paste0(as.character(TF_meta[i,"Experiment.target"]),　" peaks in ",　as.character(TF_meta[i,"Biosample.term.name"]), " cell line");
    plot(h, main=main.text, xlab="Average methylation levels of each peak", col="lightseagreen", font.axis=2, font.lab=2, cex.axis=1.5, cex.lab=1.5, lwd=2, cex.main=1.5);
    dev.off(); # cat(i, "\t", as.character(tfbs.rrbs.file$Cell[i]), "\t", as.character(TF_meta[i,"Experiment target"]), "\n");
  }
}

RRBS.TFBS.methylation.levels.distribution <- function(pathway,tfbs.rrbs.filename)
{
  TF_meta  <-read.csv(tfbs.rrbs.filename);
  all.methyl.num <- NULL;
  for(i in 1:nrow(TF_meta))
  {
    TFBS.Peak.filename <- paste0(pathway,"TFBS.Peak.methyl/", as.character(TF_meta[i,"Biosample.term.name"]), ".",as.character(TF_meta[i,"Experiment.target"]),
                                 ".", as.character(TF_meta[i,1]),".csv");
    TFBS.Peak <- read.csv(TFBS.Peak.filename);
    Peak.methyl.ratio <- TFBS.Peak$avg.nona.methyl.ratio;
    Peak.methyl.ratio <- Peak.methyl.ratio[!is.na(Peak.methyl.ratio)];
    methyl.range <- c(0:50)*0.02;
    methyl.num <- array(0, dim=(length(methyl.range)-1));
    methyl.num[1] <- sum(Peak.methyl.ratio>=methyl.range[1] & Peak.methyl.ratio<=methyl.range[2]);
    for(j in 2:length(methyl.num))
    {
      methyl.num[j] <- sum(Peak.methyl.ratio>methyl.range[j] & Peak.methyl.ratio<=methyl.range[j+1]);
    }
    all.methyl.num <- rbind(all.methyl.num, methyl.num);
    # cat(i, "\t", as.character(tfbs.rrbs.file$Cell[i]), "\t", as.character(tfbs.rrbs.file$Factor[i]), "\n");
  }
  tfbs.rrbs.methyl.num <- cbind(TF_meta[,c(1:12,17:25,29,36:42,45:47)], all.methyl.num);
  write.csv(tfbs.rrbs.methyl.num, paste0(pathway,"TFBS.RRBS.methylation.levels.distribution.csv"), row.names=FALSE);
}

###########
# main to call the function
###########
##############################step1##########################
call_TFBS.SequenceW_froPar <- function(i)
{
  print(paste0("strat " ,i))
  annotaion <- read.csv("data/annotation/annotation.fit.csv",stringsAsFactors =F )
  raw_dir<- switch (annotaion[i,"File.assembly"],
                    "GRCh38" = "Raw_hg38",
                    "hg19" = "Raw_hg19Tohg38")
  pathway <- paste0("data/TF/",annotaion[i,"Biosample.term.name"],"/",raw_dir,"/")
  TFBS.filename=paste0(pathway,annotaion[i,1],".bed")
  output.filename=paste0("Result/TFBS.Seq/",annotaion[i,"Biosample.term.name"],"_",annotaion[i,"Experiment.target"],"_",annotaion[i,1],"_","GRCh38",".csv")
  TFBS.SequenceW(TFBS.filename, output.filename,dic="/ref_hg38/hg38/")
}
library(foreach)
library(doParallel)
detectCores(logical=F)
cl <- makeCluster(4)
registerDoParallel(cl)
foreach(i = 1:2337) %dopar% call_TFBS.SequenceW_froPar(i)
stopCluster(cl)
##############################step2##########################
system("mkdir Result/TFBS.CpG.methyl/")
call_TFBS.CpG.site_froPar <- function(cell="HepG2")
{
  print(paste0("strat " ,i))
  annotaion <- read.csv("data/annotation/annotation.fit.csv",stringsAsFactors =F )
  annotaion <-  annotaion[annotaion[,"Biosample.term.name"]==cell,]
  WGBS <- read.csv(paste0("data/WGBS/",cell,".GRCh38.WGBS.csv"))
  for(i in 1:nrow(annotaion ))
  {
    TFBS.filename=paste0("Result/TFBS.Seq/",annotaion[i,"Biosample.term.name"],"_",annotaion[i,"Experiment.target"],"_",annotaion[i,1],"_","GRCh38",".csv")
    output.filename <- paste0("Result/TFBS.CpG.methyl/",annotaion[i,"Biosample.term.name"],"_",annotaion[i,"Experiment.target"],"_",annotaion[i,1],"_","GRCh38",".csv")
    TFBS.CpG.site(TFBS.filename, WGBS, output.filename)
  }
}
cl <- makeCluster(5)
registerDoParallel(cl)
foreach(i=c("GM12878","H1","HeLa-S3","K562","HepG2","A549","HEK293","MCF-7","SK-N-SH")) %dopar% call_TFBS.CpG.site_froPar(i)
stopCluster(cl)

##############################step3##########################
system("mkdir Result/TFBS.Peak.methyl/")
all_files <- list.files("Result/TFBS.CpG.methyl/")
call_TFBS.Peak.methylation_froPar <- function(i)
{
  TFBS.Peak.methylation(TFBS.filename=paste0("Result/TFBS.Seq/",all_files[i]),
                        TFBS.CpG.filename=paste0("Result/TFBS.CpG.methyl/",all_files[i]),
                        output.filename=paste0("Result/TFBS.Peak.methyl/",all_files[i]))
}
cl <- makeCluster(10)
registerDoParallel(cl)
foreach(i=c(1:length(all_files))) %dopar% call_TFBS.Peak.methylation_froPar(i)
stopCluster(cl)

