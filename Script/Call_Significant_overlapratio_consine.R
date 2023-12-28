make_annotation_add_uniqueTFnames <- function()
{
  annotaion <- read.csv("Co_factor_methyl/data/annotation/annotation.fit.csv")
  #to map the raw files names
  raw_names<-paste0(annotaion[,"Biosample.term.name"],"_",annotaion[,"Experiment.target"],"_",annotaion[,1],"_","GRCh38",".csv")
  Unique_TF_name <- rep("",length(raw_names))
  for(cell in c("GM12878","H1","HeLa-S3","K562","HepG2","A549","HEK293","SK-N-SH","MCF-7"))
  {
    High_1<- read.csv(file = paste0("Result/Summary/",cell,"_high_overlap_count_to.csv"),stringsAsFactors = F,row.names = 1)
    TF_name_list <- unlist(lapply(row.names(High_1),function(x)
    {
      temp <- unlist(strsplit(x,"-human"))[1]
      temp <- unlist(strsplit(temp,"_"))[2]
      return(temp)
    }))
    duplicated<- tapply(TF_name_list,as.factor(TF_name_list),function(x){length(x)})
    base_TF_name_list <- TF_name_list
    if(max(duplicated)>1)
    {
      for(i in 2:max(duplicated))
      {
        TF_name_list[which(duplicated(TF_name_list))] <- paste0(base_TF_name_list[which(duplicated(TF_name_list))],"_",i)
      }
    }
   match_index <-match(raw_names,row.names(High_1))
   Unique_TF_name[which(!is.na(match_index))] <- TF_name_list[match_index[which(!is.na(match_index))]]
  }
  annotaion <- cbind(annotaion,Unique_TF_name)
  write.csv(annotaion,file = "data/annotation/annotation.fit_add_names.csv")
}
#count the numbers of overlap peaks
Summary.high.low.count_to <- function(cell)
{
  library(GenomicRanges)
  if(!("package:Biostrings" %in%  search()))
  {
    library(Biostrings);
  }
  files_list <- list.files("Result/TFBS.Peak.methyl/")
  files_list <- files_list[grep(cell,files_list)]
  print(paste0("Here have ",length(files_list)," ChIP-seq dataset!"))
  low_methyl_peaks_count <- c()
  high_methyl_peaks_count  <- c()
  system_lines <-c()
  overlap_cout_matrix_high <- matrix(data=NA, nrow = length(files_list),ncol = length(files_list)+1)
  overlap_cout_matrix_low <- matrix(data=NA, nrow = length(files_list),ncol = length(files_list)+1)
  for(i in 1:length(files_list))
  {
    print(i)
    tfbs.seq.filename <- files_list[i]
    Peaks.seq <- read.csv(paste0("Result/TFBS.Peak.methyl/",tfbs.seq.filename),stringsAsFactors = F);
    Peaks.seq <-Peaks.seq[Peaks.seq[,3]-Peaks.seq[,2] < 800,]
    high.methyl.Peak <- Peaks.seq[Peaks.seq[,"avg.na.methyl.ratio"] >= 0.6 & !is.na(Peaks.seq[,"avg.na.methyl.ratio"]), ];
    low.methyl.Peak <- Peaks.seq[Peaks.seq[,"avg.na.methyl.ratio"] < 0.6 | is.na(Peaks.seq[,"avg.na.methyl.ratio"]), ];

    overlap_cout_matrix_low[i,1]<- nrow(low.methyl.Peak)
    overlap_cout_matrix_high[i,1]<- nrow(high.methyl.Peak)
    
    if(nrow(high.methyl.Peak)!=0)
    {
      high_regions <- GRanges(seqnames = high.methyl.Peak[,1],IRanges(start = high.methyl.Peak[,2],end = high.methyl.Peak[,3]))
      low_regions <- GRanges(seqnames = low.methyl.Peak[,1],IRanges(start = low.methyl.Peak[,2],end = low.methyl.Peak[,3]))
      for(j in 1:length(files_list))
      {
        print(j)
        tfbs.seq.filename_j <- files_list[j]
        Peaks.seq_j <- read.csv(paste0("Result/TFBS.Peak.methyl/",tfbs.seq.filename_j),stringsAsFactors = F);
        Peaks.seq_j <-Peaks.seq_j [Peaks.seq_j [,3]-Peaks.seq_j [,2] < 800,]
        Peaks.seq_j_regions <- GRanges(seqnames =  Peaks.seq_j [,1],IRanges(start =  Peaks.seq_j [,2],end =  Peaks.seq_j [,3]))
        high_overlap_regions <- findOverlaps(Peaks.seq_j_regions,high_regions)
        low_overlap_regions <- findOverlaps(Peaks.seq_j_regions,low_regions)
        high_overlap_count <- length(unique(high_overlap_regions@to))
        low_overlap_count <- length(unique(low_overlap_regions@to))
        overlap_cout_matrix_high[i,j+1] <-high_overlap_count
        overlap_cout_matrix_low[i,j+1] <- low_overlap_count
      }
    }else
    {
      low_regions <- GRanges(seqnames = low.methyl.Peak[,1],IRanges(start = low.methyl.Peak[,2],end = low.methyl.Peak[,3]))
      for(j in 1:length(files_list))
      {
        print(j)
        tfbs.seq.filename_j <- files_list[j]
        Peaks.seq_j <- read.csv(paste0("Result/TFBS.Peak.methyl/",tfbs.seq.filename_j),stringsAsFactors = F);
        Peaks.seq_j <-Peaks.seq_j [Peaks.seq_j [,3]-Peaks.seq_j [,2] < 800,]
        Peaks.seq_j_regions <- GRanges(seqnames =  Peaks.seq_j [,1],IRanges(start =  Peaks.seq_j [,2],end =  Peaks.seq_j [,3]))
        low_overlap_regions <- findOverlaps(Peaks.seq_j_regions,low_regions)
        low_overlap_count <- length(unique(low_overlap_regions@to))
        overlap_cout_matrix_low[i,j+1] <- low_overlap_count
      }
    }
  }
  row.names(overlap_cout_matrix_low) <- files_list 
  colnames(overlap_cout_matrix_low)<-c("total_low_peak",paste0(files_list,"_overlap") )
  
  row.names(overlap_cout_matrix_high) <- files_list 
  colnames(overlap_cout_matrix_high)<-c("total_high_peak",paste0(files_list,"_overlap") )
  
  write.csv(overlap_cout_matrix_high,file = paste0("Result/Summary/",cell,"_high_overlap_count_to.csv"),quote = F)
  write.csv(overlap_cout_matrix_low,file = paste0("Result/Summary/",cell,"_low_overlap_count_to.csv"),quote = F)
}
#get the significant overlap ratio and consion ratio
Draw_heatmap<- function(cell="GM12878")#ratio and consine
{
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  library(fitdistrplus)
  High<- read.csv(file = paste0("Result/Summary/",cell,"_high_overlap_count_to.csv"),stringsAsFactors = F,row.names = 1)
  Low<- read.csv(file = paste0("Result/Summary/",cell,"_low_overlap_count_to.csv"),stringsAsFactors = F,row.names = 1)
  TF_name_list <- unlist(lapply(row.names(High),function(x)
    {
     temp <- unlist(strsplit(x,"-human"))[1]
     temp <- unlist(strsplit(temp,"_"))[2]
     return(temp)
  }))

  duplicated<- tapply(TF_name_list,as.factor(TF_name_list),function(x){length(x)})
  base_TF_name_list <- TF_name_list
  if(max(duplicated)>1)
  {
    for(i in 2:max(duplicated))
      {
      TF_name_list[which(duplicated(TF_name_list))] <- paste0(base_TF_name_list[which(duplicated(TF_name_list))],"_",i)
      }
  }
    
  row.names(High) <- TF_name_list 
  colnames(High)<-c("total",TF_name_list )
  row.names(Low) <- TF_name_list 
  colnames(Low)<-c("total",TF_name_list )
  #start make overlap ratio and cosine matrix for High and low
  KeepID <- which(High[,1]>= 500)
  High <- High[KeepID,c(1, KeepID+1)]
  over_High <- High[,-1]
  over_High_cosine <- High[,-1]
  for(j in 1:nrow(High))
  {
    over_High[j,] <- High[j,-1]/High[j,1]
  }
  for(j in 1:nrow(High))
  {
    for(i in 1:(ncol(High)-1))
      {
      over_High_cosine[j,i] <- High[j,i+1]/(sqrt(High[j,1]*High[i,1]))
    }
  }
  
  KeepID <- which(Low[,1]>= 500)
  Low <- Low[KeepID,c(1, KeepID+1)]
  over_Low <- Low[,-1]
  over_Low_cosine <- Low[,-1]
  for(j in 1:nrow(Low))
  {
    over_Low[j,] <- Low[j,-1]/Low[j,1]
  }
  for(j in 1:nrow(Low))
  {
    for(i in 1:(ncol(Low)-1))
    {
      over_Low_cosine[j,i] <- Low[j,i+1]/(sqrt(Low[j,1])*sqrt(Low[i,1]))
    }
  }
  pdf(file=paste0("Result/Summary/",cell,"_high_overlap_count_to.pdf"),height = 12,width = 12)
  out <-  pheatmap(over_High,
           show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
           cex=1, clustering_distance_rows="euclidean", cex=1,
           clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE)
  dev.off()
  
  #get the significant overlap ratio for high,1 mains significant.
  over_High_O_1 <- over_High
  over_High_O_1[,]  <- 0
  new_cutoff <- rep(0,ncol(over_High))
  old_cutoff <- rep(0,ncol(over_High))
  error_TF_index <- NULL
  for(i in 1:ncol(over_High))
  {
    print(paste0(cell,"_High:",i,":",rownames(over_High)[i]))
    ratio_all <- as.numeric(over_High[i,])
    descdist(as.numeric(High[i,-1]))#to test the distribution #beta
    mean <- mean(as.numeric(ratio_all[as.numeric(ratio_all)<1]))
    var <- sd(as.numeric(ratio_all[as.numeric(ratio_all)<1]))
    f <- (mean*(1-mean)/var)-1
    a <- mean*f
    b <- f-a
    new_cutoff[i] <- tryCatch({ 
      t <- fitdistr(as.numeric(ratio_all[as.numeric(ratio_all)<1& as.numeric(ratio_all)>0]),"beta",start = list(shape1=a,shape2=b))
      qbeta(0.95,shape1=t$estimate[1],shape2=t$estimate[2])
      },error=function(e){
        tryCatch({
            t <- fitdistr(as.numeric(ratio_all[as.numeric(ratio_all)<1& as.numeric(ratio_all)>0]),"beta",start = list(shape1=1,shape2=9))
            qbeta(0.95,shape1=t$estimate[1],shape2=t$estimate[2])
            },error=function(e){
              ratio <- summary(ratio_all[which(ratio_all < 1)])
              ratio[4]+3*ratio[2]
              print(paste0(cell,"Use old_cutoff:",i))
            })
    })
    over_High_O_1[i,which(ratio_all > new_cutoff[i] & ratio_all < 1)]<-1
    over_High_O_1[i,i]<-0
    
    #   x_beta <- seq(0, 0.5, by = 0.01)           
    #   plot(dbeta(as.numeric(x_beta ), shape1=t$estimate[1],shape2=t$estimate[2] )) 
    #plot(dbeta(as.numeric(x_beta ), shape1=t$estimate[1],shape2=t$estimate[2] )) 
    #  new_cutoff[i] <- qbeta(0.99,shape1=t$estimate[1],shape2=t$estimate[2] )
    # },error=function(e){
    #   error_TF_index <- c(error_TF_index,i)
    #   print(paste0("Last error",i))
    # })
    # hist(as.numeric(High[i,which(as.numeric(High[i,])<1500)]))
    
    # hist(as.numeric(over_High[i,which(as.numeric(over_High[i,])<0.8)]),border = "#00BFFF",col ="#00BFFF",xlab = "Overlap_ratio", ylab = "density",breaks = 20,density = T)
    # x <- seq(0.001, 0.8, length.out = 100) 
    # y<- dbeta(x, shape1=t$estimate[1],shape2=t$estimate[2] ) 
    # lines(x, y*6, col = "red", lwd = 2)
    # abline(b=0,v=qbeta(0.95,shape1=t$estimate[1],shape2=t$estimate[2]),col="#8B008B")
    # text(0.62,15,round(qbeta(0.95,shape1=t$estimate[1],shape2=t$estimate[2]),4))
    
    # hist(as.numeric(over_High[i,which(as.numeric(over_High[i,])<0.4)]))
    # lambda1<-fitdistr(as.numeric(over_High[i,]), "Poisson")
    # new_cutoff[i] <- qpois(0.99,lambda = lambda1$estimate)
    # ratio <- summary(ratio_all[which(ratio_all < 1)])
    # old_cutoff[i]<- ratio[4]+3*ratio[2]#set the cutoff

    #outstanding_over[[i]] <- colnames(over_High)[which(ratio_all > cutoff)]
    
  }
  write.csv(over_High_O_1,file = paste0("Result/Summary/",cell,"_high_overlap_count_to_0_1.csv"),quote = F)
  pdf(file=paste0("Result/Summary/",cell,"_high_overlap_count_to_0_1.pdf"),height = 12,width = 12)
  out <-  pheatmap(over_High_O_1, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="none",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                   color=c("white","firebrick3"))
  dev.off()
  pdf(file=paste0("Result/Summary/",cell,"_high_overlap_count_to_0_1_Del0.pdf"),height = 12,width = 12)
  row_sum <- apply(over_High_O_1, 1, sum)
  col_sum <- apply(over_High_O_1, 2, sum)
  new_over_High_O_1 <- over_High_O_1[which(row_sum>0),which(col_sum>0)]
  out <-  pheatmap(new_over_High_O_1, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="none",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                   color=c("white","firebrick3"))
  out1<- pheatmap(new_over_High_O_1,  cluster_cols=T, cluster_rows=T, scale="none",
                 cex=1, clustering_distance_rows="euclidean", cex=1,
                 clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                 color=c("#E3E3E3","firebrick3"),legend = FALSE,treeheight_col = 0,treeheight_row = 0,
                 show_colnames = F,show_rownames=F)
  dev.off()
  #end the significant overlap ratio for high
  
  #get the significant cosine  ratio for high,1 mains significant.
  pdf(file=paste0("Result/Summary/",cell,"_high_overlap_count_to_cosine.pdf"),height = 12,width = 12)
  out <-  pheatmap(over_High_cosine,
                   show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE)
  dev.off()
  over_High_0_1_cosine <- over_High_cosine
  over_High_0_1_cosine[,]  <- 0
  descdist(over_High_cosine[over_High_cosine<1 & over_High_cosine >0])
  mean <- mean(over_High_cosine[over_High_cosine<1 & over_High_cosine >0])
  var <- sd(over_High_cosine[over_High_cosine<1 & over_High_cosine >0])
  f <- (mean*(1-mean)/var)-1
  a <- mean*f
  b <- f-a
  new_cutoff <- tryCatch({ 
    t <- fitdistr(over_High_cosine[over_High_cosine<1 & over_High_cosine >0],"beta",start = list(shape1=a,shape2=b))
    qbeta(0.95,shape1=t$estimate[1],shape2=t$estimate[2])
    },error=function(e){
      tryCatch({
        t <- fitdistr(over_High_cosine[over_High_cosine<1 & over_High_cosine >0],"beta",start = list(shape1=1,shape2=9))
        qbeta(0.95,shape1=t$estimate[1],shape2=t$estimate[2])
        },error=function(e){
          ratio <- summary(over_High_cosine[over_High_cosine<1])
          ratio[4]+3*ratio[2]
          print(paste0(cell,"Use old_cutoff:",ratio[4]+3*ratio[2]))
        })
    })
  for(i in 1:nrow(over_High_0_1_cosine))
  {
    for(j in 1:ncol(over_High_0_1_cosine))
    {
      if(over_High_cosine[i,j] >= new_cutoff & over_High_cosine[i,j]!=  1)
      {
        over_High_0_1_cosine[i,j] = 1
      }
    }
    over_High_0_1_cosine[i,i]<-0
  }
  write.csv(over_High_0_1_cosine,file = paste0("Result/Summary/",cell,"_high_overlap_count_to_0_1_consine.csv"),quote = F)
  pdf(file=paste0("Result/Summary/",cell,"_high_overlap_count_to_0_1_cosine.pdf"),height = 12,width = 12)
  out <-  pheatmap(over_High_0_1_cosine, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="none",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                   color=c("white","firebrick3"))
  dev.off()
  pdf(file=paste0("Result/Summary/",cell,"_high_overlap_count_to_0_1_consine_Del0.pdf"),height = 12,width = 12)
  row_sum <- apply(over_High_0_1_cosine, 1, sum)
  col_sum <- apply(over_High_0_1_cosine, 2, sum)
  new_over_High_0_1_cosine <- over_High_0_1_cosine[which(row_sum>0),which(col_sum>0)]
  out <-  pheatmap(new_over_High_0_1_cosine, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="none",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                   color=c("white","firebrick3"))
  out1_2<- pheatmap(new_over_High_0_1_cosine,  cluster_cols=T, cluster_rows=T, scale="none",
                  cex=1, clustering_distance_rows="euclidean", cex=1,
                  clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                  color=c("#E3E3E3","firebrick3"),legend = FALSE,treeheight_col = 0,treeheight_row = 0,
                  show_colnames = F,show_rownames=F)
  dev.off()
  #end high cosine
  
 
  #get the significant overlap ratio for low,1 mains significant.
  pdf(file=paste0("Result/Summary/",cell,"_low_overlap_count_to.pdf"),height = 18,width = 18)
  out <-  pheatmap(over_Low, 
                   show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE)
  dev.off()
  pdf(file=paste0("Result/Summary/",cell,"_low_overlap_count_to_0_1.pdf"),height = 12,width = 12)
  over_Low_O_1 <- over_Low
  over_Low_O_1[,]  <- 0
  new_cutoff <- rep(0,ncol(over_Low))
  old_cutoff <- rep(0,ncol(over_Low))
  error_TF_index <- NULL
  for(i in 1:ncol(over_Low))
  {
    print(paste0(cell,"_Low:",i,":",rownames(over_Low)[i]))
    ratio_all <- as.numeric(over_Low[i,])
    descdist(as.numeric(Low[i,-1]))#to test the distrubution#beta
    mean <- mean(as.numeric(ratio_all[as.numeric(ratio_all)<1]))
    var <- sd(as.numeric(ratio_all[as.numeric(ratio_all)<1]))
    f <- (mean*(1-mean)/var)-1
    a <- mean*f
    b <- f-a
    new_cutoff[i] <- tryCatch(
      { 
        t <- fitdistr(as.numeric(ratio_all[as.numeric(ratio_all)<1& as.numeric(ratio_all)>0]),"beta",start = list(shape1=a,shape2=b))
        qbeta(0.95,shape1=t$estimate[1],shape2=t$estimate[2])
      },error=function(e)
      {
        tryCatch({
          t <- fitdistr(as.numeric(ratio_all[as.numeric(ratio_all)<1& as.numeric(ratio_all)>0]),"beta",start = list(shape1=1,shape2=9))
          qbeta(0.95,shape1=t$estimate[1],shape2=t$estimate[2])
          },error=function(e)
          {
            ratio <- summary(ratio_all[which(ratio_all < 1)])
            ratio[4]+3*ratio[2]
            print(paste0(cell,"Use old_cutoff:",i))
          })
      })
    ratio <- summary(ratio_all[which(ratio_all < 1)])
    old_cutoff[i]<- ratio[4]+3*ratio[2]#set the cutoff

    over_Low_O_1[i,which(ratio_all > new_cutoff[i] & ratio_all !=  1)]<-1
    over_Low_O_1[i,i]<-0
  }
  write.csv(over_Low_O_1,file = paste0("Result/Summary/",cell,"_low_overlap_count_to_0_1.csv"),quote = F)
  out <-  pheatmap(over_Low_O_1, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="none",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                   color=c("white","firebrick3"))
  dev.off()
  pdf(file=paste0("Result/Summary/",cell,"_low_overlap_count_to_0_1_Del0.pdf"),height = 12,width = 12)
  row_sum <- apply(over_Low_O_1, 1, sum)
  col_sum <- apply(over_Low_O_1, 2, sum)
  new_over_Low_O_1 <- over_Low_O_1[which(row_sum>0),which(col_sum>0)]
  out <-  pheatmap(new_over_Low_O_1, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="none",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                   color=c("white","firebrick3"))
  out2<- pheatmap(new_over_Low_O_1, cluster_cols=T, cluster_rows=T, scale="none",
                  cex=1, clustering_distance_rows="euclidean", cex=1,
                  clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                  color=c("#E3E3E3","firebrick3"),legend = FALSE,treeheight_col = 0,treeheight_row = 0,
                  show_colnames = F,show_rownames=F)

  dev.off()
  #end the overlap ratio for low
  
  #get the significant cosine ratio for low,1 mains significant.
  remove(over_High_cosine)
  remove(over_High_0_1_cosine)
  pdf(file=paste0("Result/Summary/",cell,"_Low_overlap_count_to_cosine.pdf"),height = 12,width = 12)
  out <-  pheatmap(over_Low_cosine,
                   show_rownames=T, cluster_cols=T, cluster_rows=T, scale="row",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color=FALSE)
  dev.off()
  over_Low_0_1_cosine <- over_Low_cosine
  over_Low_0_1_cosine[,]  <- 0
  descdist(over_Low_cosine[over_Low_cosine<1 & over_Low_cosine >0])
  mean <- mean(over_Low_cosine[over_Low_cosine<1 & over_Low_cosine >0])
  var <- sd(over_Low_cosine[over_Low_cosine<1 & over_Low_cosine >0])
  f <- (mean*(1-mean)/var)-1
  a <- mean*f
  b <- f-a
  new_cutoff <- tryCatch(
    { 
      t <- fitdistr(over_Low_cosine[over_Low_cosine<1 & over_Low_cosine >0],"beta",start = list(shape1=a,shape2=b))
      qbeta(0.95,shape1=t$estimate[1],shape2=t$estimate[2])
    },error=function(e)
    {
      tryCatch({
        t <- fitdistr(over_Low_cosine[over_Low_cosine<1 & over_Low_cosine >0],"beta",start = list(shape1=1,shape2=9))
        qbeta(0.95,shape1=t$estimate[1],shape2=t$estimate[2])
        },error=function(e)
        {
          ratio <- summary(over_Low_cosine[over_Low_cosine<1])
          ratio[4]+3*ratio[2]
          print(paste0(cell,"Use old_cutoff:",ratio[4]+3*ratio[2]))
        })
    })
  for(i in 1:nrow(over_Low_0_1_cosine))
  {
    for(j in 1:ncol(over_Low_0_1_cosine))
    {
      if(over_Low_cosine[i,j] >= new_cutoff & over_Low_cosine[i,j]!= 1)
      {
        over_Low_0_1_cosine[i,j] = 1
      }
    }
    over_Low_0_1_cosine[i,i]<-0
  }
  write.csv(over_Low_0_1_cosine,file = paste0("Result/Summary/",cell,"_Low_overlap_count_to_0_1_consine.csv"),quote = F)
  pdf(file=paste0("Result/Summary/",cell,"_Low_overlap_count_to_0_1_cosine.pdf"),height = 12,width = 12)
  out <-  pheatmap(over_Low_0_1_cosine, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="none",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                   color=c("white","firebrick3"))
  dev.off()
  pdf(file=paste0("Result/Summary/",cell,"_Low_overlap_count_to_0_1_consine_Del0.pdf"),height = 12,width = 12)
  row_sum <- apply(over_Low_0_1_cosine, 1, sum)
  col_sum <- apply(over_Low_0_1_cosine, 2, sum)
  new_over_Low_0_1_cosine <- over_Low_0_1_cosine[which(row_sum>0),which(col_sum>0)]
  out <-  pheatmap(new_over_Low_0_1_cosine, show_rownames=T, cluster_cols=T, cluster_rows=T, scale="none",
                   cex=1, clustering_distance_rows="euclidean", cex=1,
                   clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                   color=c("white","firebrick3"))
  out2_2<- pheatmap(new_over_Low_0_1_cosine,  cluster_cols=T, cluster_rows=T, scale="none",
                  cex=1, clustering_distance_rows="euclidean", cex=1,
                  clustering_distance_cols="euclidean", clustering_method="complete", border_color="black",
                  color=c("#E3E3E3","firebrick3"),legend = FALSE,treeheight_col = 0,treeheight_row = 0,
                  show_colnames = F,show_rownames=F)
  dev.off()
  #end low cosine
  
  # require(ggplotify)
  # g1 = as.ggplot(out1)
  # #colnames(new_over_High_O_1)[out1$tree_col$order][1:10]
  # g2 = as.ggplot(out2)
  # #colnames(new_over_Low_O_1)[out2$tree_col$order][1:10]
  # pdf(file=paste0("Result/Write_paper/",cell,"_figure.pdf"),height = 4,width = 12)
  # cowplot::plot_grid(g1, g2, ncol=2)
  # dev.off()
  # 
  # g1 = as.ggplot(out1_2)
  # g2 = as.ggplot(out2_2)
  # pdf(file=paste0("Result/Write_paper/",cell,"_figure_consion.pdf"),height = 4,width = 12)
  # cowplot::plot_grid(g1, g2, ncol=2)
  # dev.off()
}

#it is used to plot fo publish paper
plot_barball<-function(cell,TF)
{
  cell="K562"
  #TF="ZBTB33"
  TF = "ATF3"
  High<- read.csv(file = paste0("Result/Summary/",cell,"_high_overlap_count_to.csv"),stringsAsFactors = F,row.names = 1)
  Low<- read.csv(file = paste0("Result/Summary/",cell,"_low_overlap_count_to.csv"),stringsAsFactors = F,row.names = 1)

  
  TF_name_list <-row.names(High)
  TF_name_list <- unlist(lapply(row.names(High),function(x)
  {
    temp <- unlist(strsplit(x,"-human"))[1]
    temp <- unlist(strsplit(temp,"_"))[2]
    return(temp)
  }))
  select_index <- which(TF_name_list ==TF)

  duplicated<- tapply(TF_name_list,as.factor(TF_name_list),function(x){length(x)})
  base_TF_name_list <- TF_name_list
  if(max(duplicated)>1)
  {
    for(i in 2:max(duplicated))
    {
      TF_name_list[which(duplicated(TF_name_list))] <- paste0(base_TF_name_list[which(duplicated(TF_name_list))],"_",i)
    }
  }
  row.names(High) <- TF_name_list
  colnames(High)<-c("total",TF_name_list )
  row.names(Low) <- TF_name_list
  colnames(Low)<-c("total",TF_name_list )
  
  over_High <- High[,-1]
  for(j in 1:nrow(High))
  {
    over_High[j,] <- High[j,-1]/High[j,1]
  }
  
  over_Low <- Low[,-1]
  for(j in 1:nrow(Low))
  {
    over_Low[j,] <- Low[j,-1]/Low[j,1]
  }

  over_High_O_1 <- over_High
  over_High_O_1[,]  <- 0
  for(i in 1:ncol(over_High))
  {
    ratio_all <- as.numeric(over_High[i,])
    ratio <- summary(ratio_all[which(ratio_all < 1)])
    cutoff <- ratio[4]+3*ratio[2]#set the cutoff
    #outstanding_over[[i]] <- colnames(over_High)[which(ratio_all > cutoff)]
    over_High_O_1[i,which(ratio_all > cutoff)]<-1
    over_High_O_1[i,i]<-0
  }
  over_Low_O_1 <- over_Low
  over_Low_O_1[,]  <- 0
  for(i in 1:ncol(over_Low))
  {
    ratio_all <- as.numeric(over_Low[i,])
    ratio <- summary(ratio_all[which(ratio_all < 1)])
    cutoff <- ratio[4]+3*ratio[2]
    #outstanding_over[[i]] <- colnames(over_High)[which(ratio_all > cutoff)]
    over_Low_O_1[i,which(ratio_all > cutoff)]<-1
    over_Low_O_1[i,i]<-0
  }
  

  h_index <- grep(TF,row.names(over_High_O_1))[1]
  l_index <- grep(TF,row.names(over_High_O_1))[1]
  table_t <- matrix(NA,nrow = 2,ncol = length(TF_name_list))
  table_t_0_1<- matrix(NA,nrow = 2,ncol = length(TF_name_list))
  table_t_count<- matrix(NA,nrow = 2,ncol = length(TF_name_list))
  for(i in 1:length(TF_name_list))
  {
    if(length(which(colnames(over_High)==TF_name_list[i]))!=0)
    {
      table_t[1,i] <- over_High[h_index,which(colnames(over_High)==TF_name_list[i])]
      table_t_0_1[1,i] <- over_High_O_1[h_index,which(colnames(over_High_O_1)==TF_name_list[i])]
      table_t_count[1,i] <- High[h_index,which(colnames(High)==TF_name_list[i])]
    }
    if(length(which(colnames(over_Low)==TF_name_list[i]))!=0)
    {
      table_t[2,i] <- over_Low[h_index,which(colnames(over_Low)==TF_name_list[i])]
      table_t_0_1[2,i] <- over_Low_O_1[h_index,which(colnames(over_Low_O_1)==TF_name_list[i])]
      table_t_count[2,i] <- Low[h_index,which(colnames(Low)==TF_name_list[i])]
    }
  }
  colnames(table_t) <- TF_name_list
  colnames(table_t_0_1) <- TF_name_list
  colnames(table_t_count) <- TF_name_list
  
  x <- c(rep("Meth",4),rep("Un_Me",4))#c(1,1,1,1,2,2,2,2)#c(rep("Meth",4),rep("Un_Me",4))
  y <- c(c("AGO1","ATF3","CTCF","RBFOX2"),c("AGO1","ATF3","CTCF","RBFOX2"))#c(1,2,3,4,1,2,3,4)#
  Ovarlap_ratio <- c(table_t[1,c("AGO1","ATF3","CTCF","RBFOX2")],
         table_t[2,c("AGO1","ATF3","CTCF","RBFOX2")])
  Ovarlap_count <- c( table_t_count[1,c("AGO1","ATF3","CTCF","RBFOX2")],
                      table_t_count[2,c("AGO1","ATF3","CTCF","RBFOX2")])
  Ovarlap_count_0_1 <- 
  TF <- rep("ZBTB33",8)
  data_ball<-data.frame(x,y,Ovarlap_ratio,Ovarlap_count,TF)#cpunt
  ggplot(data_ball,aes(x=x,y=y,size=Ovarlap_ratio))+geom_point(colour="#a30eed")+theme_bw()+
    theme(panel.grid.major = element_blank(),axis.text = element_text(colour = "black") )
  # ggplot(data_ball,aes(x=x,y=y,size=Ovarlap_count))+geom_point(colour="#a30eed",alpha=Ovarlap_ratio)+theme_bw()+
  #   theme(panel.grid.major = element_blank(),axis.text = element_text(colour = "black") )
  # # which(table_t[1,]>0.2): ATF3   CTCF_2   POLR2G   RBFOX2
  
  
  
  x <- c(rep("Meth",4),rep("Un_Me",4))#c(1,1,1,1,2,2,2,2)#c(rep("Meth",4),rep("Un_Me",4))
  y <- c(c("ATF4","CEBPB_2","CEBPG","JUND"),c("ATF4","CEBPB_2","CEBPG","JUND"))#c(1,2,3,4,1,2,3,4)#
  Ovarlap_ratio <- c(table_t[1,c("ATF4","CEBPB_2","CEBPG","JUND")],
                     table_t[2,c("ATF4","CEBPB_2","CEBPG","JUND")])
  TF <- rep("ATF3",8)
  data_ball_ATF3<-data.frame(x,y,Ovarlap_ratio,TF)
  ggplot(data_ball,aes(x=x,y=y,size=Ovarlap_ratio))+geom_point(colour="orange")+theme_bw()+
    theme(panel.grid.major = element_blank(),axis.text = element_text(colour = "black") )
  
  All_data <- data.frame(x=c(data_ball_ATF3$x,data_ball$x),
                         y = c(data_ball_ATF3$y,data_ball$y),
                         Ovarlap_ratio=c(data_ball_ATF3$Ovarlap_ratio,data_ball$Ovarlap_ratio),
                         TF=c(data_ball_ATF3$TF,data_ball$TF))
  All_data$y<- factor(All_data$y, levels=c("ATF4","CEBPB_2","CEBPG","JUND","AGO1","ATF3","CTCF","RBFOX2"))
  ggplot(All_data,aes(x=x,y=y,size=Ovarlap_ratio,color=factor(TF)))+geom_point()+theme_bw()+
    theme(panel.grid.major = element_blank(),axis.text = element_text(colour = "black") )+ 
    scale_color_manual(
      breaks = c("ATF3", "ZBTB33"), 
      values = c("darkorange", "purple")
    )
}

for(cell in c("GM12878","H1","HeLa-S3","K562","HepG2","A549","HEK293","SK-N-SH","MCF-7"))
{
  Summary.high.low.count_to(cell)
  Draw_heatmap(cell)
}




# library("cluster")
# library("dendextend")
# dist.eucl <- dist(over_High_O_1,method = "euclidean")#pearson
# res_hc <- hclust(d = dist.eucl, method = "complete")
# 
# ncoldata<- t(over_High_O_1)
# dist.eucl2 <- dist(t(over_High_O_1),method = "euclidean")#pearson
# res_hc2 <- hclust(d = dist.eucl2, method = "complete")
# cluster <- cutree(res_hc2, k = 2)
# factoextra::fviz_dend(res_hc2, k = 2,rect = TRUE)
# factoextra::fviz_cluster(list(data = t(over_High_O_1), cluster = cluster))
# ncoldata_1 <- ncoldata[names(cluster[cluster==1]),]
# rowdata_1 <- over_High_O_1[,names(cluster[cluster==1])]
# dist.eucl2_1 <- dist(ncoldata_1,method = "euclidean")#pearson
# res_hc2_1 <- hclust(d = dist.eucl2_1, method = "complete")
# plot(res_hc2_1 )
# 
# dend1 <- stats::as.dendrogram(res_hc) 
# dend2 <- stats::as.dendrogram(res_hc2)
# dend2_1 <- stats::as.dendrogram(res_hc2_1)
# dend_list <- dendextend::dendlist(dend1, dend2)
# dend_list <- dendextend::dendlist(dend1, dend2_1)
# pdf(file = "temp.pdf",height = 10,width = 6)
# tanglegram(dend1, dend2_1,lwd=1,k_labels =9,k_branches = 9,
#                        edge.lwd=1,
#                        lab.cex=1,
#                        margin_inner=6,
#                        highlight_distinct_edges=F,
#                        highlight_branches_lwd=F,
#                        )
# dev.off()
# dendextend::entanglement(dend1, dend2)
# dendextend::cor.dendlist(dend_list, method = "cophenetic")
# 
# 
# res_coph <- cophenetic(res_hc)#聚类树的共同距离和原始的距离矩阵的相似性来衡量聚类的好坏：
# cor(dist.eucl, res_coph)#0.8086
# cluster <- cutree(res_hc, k = 4)
# factoextra::fviz_dend(res_hc, k = 4,rect = TRUE)
# factoextra::fviz_cluster(list(data =over_High_O_1, cluster = cluster))


#To check the gap between the regions in GM12878
# regions <-  GRanges(seqnames =  t[,1],IRanges(start =  t [,2],end = t[,3]))
# regions[order(regions)]
# gap <- c()
# for(i in 1:(length(regions)-1))
# {
#   gap[i] <-  regions[i+1]@ranges@start-(regions[i]@ranges@start+regions[i]@ranges@width )
# }

