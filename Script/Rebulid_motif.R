PPM2PWM <- function(motif_file,
                    PWM_file)
{
  motif_PPM <- read.delim(motif_file)
  motif_PPM <- t( motif_PPM[,c(1:4)]) 
  NTperc = matrix(c(.295,.205,.205,.295));		# percentage of nucleotide in A, C, G, and T
  pseudo = NTperc %*% sqrt(colSums(motif_PPM))
  PWM = log2(motif_PPM / pseudo)
  save(PWM,file=PWM_file);
}

matrix.motif <- function (motif = NULL)
{
  matrix = matrix(0, nrow = 4, ncol = nchar(motif));
  if (gregexpr("[aA]", motif)[[1]][1]>0)
  {	matrix[1,gregexpr("[aA]", motif)[[1]]] = 1;	}
  if (gregexpr("[cC]", motif)[[1]][1]>0)
  {	matrix[2,gregexpr("[cC]", motif)[[1]]] = 1;	}
  if (gregexpr("[gG]", motif)[[1]][1]>0)
  {	matrix[3,gregexpr("[gG]", motif)[[1]]] = 1;	}
  if (gregexpr("[tT]", motif)[[1]][1]>0)
  {	matrix[4,gregexpr("[tT]", motif)[[1]]] = 1;	}
  matrix;
}

single.TFBS.sequence.score <- function(tfbs.seq.filename="~/Documents/Experiments/DNA_metht_TF/Result/Homer/Low_0_6/E14_ARID3A-human_ENCFF829MGY_GRCh38/E14_ARID3A-human_ENCFF829MGY_GRCh38.bed",
                                       PWM_file, 
                                       output.file=NULL,extent_bp=2)
{
  if(!("package:Biostrings" %in%  search()))
  {
    library(Biostrings);
  }
  load(PWM_file);#PWM
  len.pssm= ncol(PWM)
  tfbs <- read.delim(tfbs.seq.filename, header=FALSE)
  
  # high peak
  first.score <- array(-100, dim = nrow(tfbs));
  first.position <- array(-100, dim = nrow(tfbs));
  first.sequence <- array("", dim = nrow(tfbs));
  first.sequence_extent <- array("", dim = nrow(tfbs));
  first.strand <- array("", dim = nrow(tfbs));
  sequence <- as.character(tfbs[,22]);
  
  for(i in 1:length(sequence))
  {
    seq <- as.character(sequence[i]);
    reverse.seq <- as.character(reverseComplement(DNAString(sequence[i])));
    m.length <- nchar(seq)-len.pssm+1;
    seq.length <- nchar(seq);
    matrix.score <- matrix(-100, nrow = 2, ncol = m.length);
    
    if (m.length>0)
    {
      seq.matrix <- matrix.motif(seq);
      reverse.seq.matrix <- matrix.motif(as.character(reverse.seq));
      for (m in 1:m.length)
      {
        motif.matrix <- seq.matrix[,m:(m+len.pssm-1)];
        reverse.motif.matrix <- reverse.seq.matrix[,(seq.length-len.pssm-m+2):(seq.length-m+1)];
        matrix.score[1,m] <- round(sum(motif.matrix * PWM),2);
        matrix.score[2,m] <- round(sum(reverse.motif.matrix * PWM),2);
      } 
      
      position <- order(matrix.score, decreasing=TRUE)[1];
      if(position %% 2 != 0)
      {
        first.score[i] <- matrix.score[position];
        first.position[i] <- floor(position/2)+1;
        first.sequence_extent[i] <- substr(seq, first.position[i]-extent_bp, first.position[i]+len.pssm-1+extent_bp);
        first.sequence[i] <- substr(seq, first.position[i], first.position[i]+len.pssm-1);
        first.strand[i] <- "+";
      } else {
        first.score[i] <- matrix.score[position];
        first.position[i] <- floor(position/2);
        first.sequence_extent[i] <- substr(reverse.seq, seq.length-len.pssm-first.position[i]+2-extent_bp, seq.length-first.position[i]+1+extent_bp);
        first.sequence[i] <- substr(reverse.seq, seq.length-len.pssm-first.position[i]+2, seq.length-first.position[i]+1);
        first.strand[i] <- "-";
      }
    }
  }
  tfbs <- cbind(tfbs, first.score, first.position, first.sequence, first.strand);
  write.csv(tfbs, file=output.file, row.names=FALSE,col.names = FALSE);
}

call_single.TFBS.sequence.score<- function(name)
{
  #name<-'62_ADNP-human_ENCFF083UZC_GRCh38"
  common_name<-name
  path<-'Result/Homer/'#Result/Homer/High_0_6
  Low_motif_name<-paste0(path,'Low_0_6/',common_name,"/homerResults/motif1.motif")
  High_motif_name<-paste0(path,'High_0_6/',common_name,"/homerResults/motif1.motif")
  
  Low_PWM_file<-paste0(path,'Low_PWM/',common_name,'.Rdata')
  High_PWM_file<-paste0(path,'High_PWM/',common_name,'.Rdata')
  if(file.exists(High_motif_name))
  {
    PPM2PWM(High_motif_name,High_PWM_file) 
  }
  PPM2PWM(Low_motif_name,Low_PWM_file)
  
  Low_tfbs.seq.filename<-paste0(path,'Low_0_6/',common_name,"/",common_name,'.bed')
  High_tfbs.seq.filename<-paste0(path,'High_0_6/',common_name,"/",common_name,'.bed')
  
  Low_output.file<-paste0(path,'Low_bindingsite_bed/',common_name,'_match_bindingsites.bed')
  High_output.file<-paste0(path,'High_bindingsite_bed/',common_name,'_match_bindingsites.bed')
  single.TFBS.sequence.score(Low_tfbs.seq.filename,Low_PWM_file,Low_output.file)
  if(!file.exists(High_PWM_file))
  {
    print('!exists high')
    High_PWM_file<-Low_PWM_file
  }
  single.TFBS.sequence.score(High_tfbs.seq.filename,High_PWM_file,High_output.file)
} 
call_single.TFBS.sequence.score_cobind<- function(name)
{
  #name<-'62_ADNP-human_ENCFF083UZC_GRCh38"
  common_name<-name
  path<-'Result/Homer/'#Result/Homer/High_0_6
  motif_name<-paste0(path,'Co_bind/',common_name,"/homerResults/motif1.motif")
  if(!file.exists(motif_name))
  {
    print(paste0("No files ",motif_name)) 
  }else{
  PWM_file<-paste0(path,'Co_bind_PWM/',common_name,'.Rdata')
  PPM2PWM(motif_name,PWM_file)
  tfbs.seq.filename<-paste0(path,'Co_bind/',common_name,"/",common_name,'.bed')
  output.file<-paste0(path,'Co_bind_bindingsite_bed/',common_name,'_match_bindingsites.bed')
  single.TFBS.sequence.score(tfbs.seq.filename,PWM_file,output.file)
  }
} 
#step 1

#library(foreach)
#library(doParallel)
#cl <- makeCluster(4)
#registerDoParallel(cl)
# summary_file="Result/Summary/Strip_TF_all_involved_fileID.csv"
# TF_involved <- read.csv(file = summary_file)
# all_files <- paste0(TF_involved[,"Biosample.term.name"],"_",TF_involved[,"TF"],"-human_",TF_involved[,"File.accession"],"_GRCh38")
# #foreach(i=all_files) %dopar% call_single.TFBS.sequence.score(i)
# #stopCluster(cl)
# for(i in all_files)
# {
#   print(i)
#   call_single.TFBS.sequence.score(i)
# }

all_files <- list.files(path = "Result/Homer/Co_bind/")
for(i in all_files)
{
  print(i)
  call_single.TFBS.sequence.score_cobind(i)
}



#step2
CpG.num.and.ratio.in.TFBS.logo_3200 <- function(PWM_file="Result/Homer/High_0_6/SK-N-SH_BACH1-human_ENCFF245FMD_GRCh38/homerResults/motif1.motif_PWM.Rdata",
                                                bindingsite_file="Result/Homer/High_0_6/SK-N-SH_BACH1-human_ENCFF245FMD_GRCh38/SK-N-SH_BACH1-human_ENCFF245FMD_GRCh38_match_bindingsites.bed",
                                                Output_index= "Result/Homer/High_0_6/SK-N-SH_BACH1-human_ENCFF245FMD_GRCh38/SK-N-SH_BACH1-human_ENCFF245FMD_GRCh38_match_bindingsites",
                                                referdic = "~/Documents/ref/",methyl.filename="data/WGBS/SK-N-SH.GRCh38.WGBS.csv")
{
  load(PWM_file)
  motif.length <- ncol(PWM)
  if(!("package:Biostrings" %in%  search()))
  {
    library(Biostrings);
  }
  methyl.Peak <- read.csv(bindingsite_file);
  colnames(methyl.Peak)[1:22] <- c("Chrom" ,"Start","End" ,"Name" ,"Score","Strand","Signalvalue","Pvalue" , "Qvalue","Peak","CpG.num" ,"methyl.CpG.num","na.CpG.num","max.read.num","max.methyl.num",
                                   "min.read.num","min.methyl.num","all.read.num","all.methyl.num","avg.nona.methyl.ratio", "avg.na.methyl.ratio","Seq")
  methyl.Peak <- methyl.Peak[order(methyl.Peak[,"Start"]), ];
  methyl.Peak <- methyl.Peak[order(methyl.Peak[,"Chrom"]), ];
  motif.CpG.num.filename <- paste0(Output_index,".motif.CpG.num.",as.character(3200+motif.length),"bp.csv");
  motif.CpG.ratio.filename <-paste0(Output_index,".motif.CpG.ratio.",as.character(3200+motif.length),"bp.csv");
  
  tmp.chr <- "";
  motif.CpG.num <- matrix(0, nrow=nrow(methyl.Peak), ncol=3200+motif.length);###??
  motif.CpG.ratio <- matrix(0, nrow=nrow(methyl.Peak), ncol=3200+motif.length);###??
  
  #dic <- paste0(referdic,"/chroms/");
  
  methylall <- read.csv(methyl.filename);
  tmp.chr <- ""
  for(i in 1:nrow(methyl.Peak))
  {
    
    if(tmp.chr != as.character(methyl.Peak$Chrom[i]))
    {
      tmp.chr <- as.character(methyl.Peak$Chrom[i]);
      s <- readDNAStringSet(paste0(referdic, tmp.chr, '.fa'));
      methyl.chr <- methylall[which(as.character(methylall[,1])==tmp.chr),]
      methyl.ratio <- methyl.chr$methy.read/methyl.chr$all.read;
      methyl.chr <- data.frame(methyl.chr$position, methyl.ratio);
      colnames(methyl.chr) <- c("position", "methyl.ratio");
      methyl.chr <- methyl.chr[order(methyl.chr[,1]), ];
    }
    # retrieve motif region		
    if(methyl.Peak$first.strand[i]=="+")
    {
      region.start <- methyl.Peak$Start[i]+methyl.Peak$first.position[i]-1-1600;
      region.end <- methyl.Peak$Start[i]+methyl.Peak$first.position[i]-2+1600+motif.length;
      if (region.end>s@ranges@width){
        next
      }
      # CpG number 
      motifSeq <- toString(subseq(s, start=region.start, end=region.end+1));
      CpG.pos <- unlist(gregexpr("CG", motifSeq));
      if (CpG.pos[1]!=-1){motif.CpG.num[i, CpG.pos] <- 1;}
      # CpG methylation ratio
      motif.methyl <- methyl.chr[methyl.chr$position>=region.start-1 & methyl.chr$position<=region.end-1, ];
      if(nrow(motif.methyl)>0){motif.CpG.ratio[i, motif.methyl$position-region.start+2] <- motif.methyl$methyl.ratio;}
      
    }
    if(methyl.Peak$first.strand[i]=="-")
    {
      region.start <- methyl.Peak$Start[i]+methyl.Peak$first.position[i]-1-1600;
      region.end <- methyl.Peak$Start[i]+methyl.Peak$first.position[i]-2+1600+motif.length;
      if (region.end>s@ranges@width){
        next
      }
      # CpG number
      motifSeq <- toString(subseq(s, start=region.start-1, end=region.end));
      motifSeq <- as.character(reverseComplement(DNAStringSet(motifSeq)));
      CpG.pos <- unlist(gregexpr("CG", motifSeq));
      if (CpG.pos[1]!=-1){motif.CpG.num[i, CpG.pos] <- 1;}
      # CpG methylation ratio
      motif.methyl <- methyl.chr[methyl.chr$position>=region.start-2 & methyl.chr$position<=region.end-2, ];
      if(nrow(motif.methyl)>0){motif.CpG.ratio[i, region.end-motif.methyl$position-1] <- motif.methyl$methyl.ratio;}
    }
    cat(i, "\t", as.character(methyl.Peak$Chrom[i]), "\t", as.character(methyl.Peak$Start[i]), "\n");
  }
  write.csv(motif.CpG.num, file=motif.CpG.num.filename, row.names=FALSE);
  write.csv(motif.CpG.ratio, file=motif.CpG.ratio.filename, row.names=FALSE);
}

call_CpG.num.and.ratio.in.TFBS.logo <- function(name)
{
  #name<-'62_ADNP-human_ENCFF083UZC_GRCh38.motif'
  common_name<-name
  path<-'Result/Homer/'
  referdic<-'/home/zqluoximei/ref_hg38/hg38/'
  cell<-unlist(strsplit(name,"_"))[1]
  methyl.filename<-paste0("data/WGBS/",cell,".GRCh38.WGBS.csv")
  Low_PWM_file<-paste0(path,'Low_PWM/',common_name,'.Rdata')
  High_PWM_file<-paste0(path,'High_PWM/',common_name,'.Rdata')

  Low_bindingsite_file<-paste0(path,'Low_bindingsite_bed/',common_name,'_match_bindingsites.bed')
  High_bindingsite_file<-paste0(path,'High_bindingsite_bed/',common_name,'_match_bindingsites.bed')

  Low_Output_index<-paste0(path,'Low_bindingsite_bed/',common_name,'_match_bindingsites')
  High_Output_index<-paste0(path,'High_bindingsite_bed/',common_name,'_match_bindingsites')
  
  if(!file.exists(High_PWM_file))
  {
    print('exists')
    High_PWM_file<-Low_PWM_file
    
  }
  CpG.num.and.ratio.in.TFBS.logo_3200(High_PWM_file,High_bindingsite_file,High_Output_index, referdic,methyl.filename)
  CpG.num.and.ratio.in.TFBS.logo_3200(Low_PWM_file,Low_bindingsite_file,Low_Output_index, referdic,methyl.filename)
  
}
call_CpG.num.and.ratio.in.TFBS.logo_cobind <- function(name)
{
  #name<-'62_ADNP-human_ENCFF083UZC_GRCh38.motif'
  common_name<-name
  path<-'Result/Homer/'
  referdic<-'/home/zqluoximei/ref_hg38/hg38/'
  cell<-unlist(strsplit(name,"_"))[1]
  methyl.filename<-paste0("data/WGBS/",cell,".GRCh38.WGBS.csv")
  PWM_file<-paste0(path,'Co_bind_PWM/',common_name,'.Rdata')
  
  
  bindingsite_file<-paste0(path,'Co_bind_bindingsite_bed/',common_name,'_match_bindingsites.bed')
  Output_index<-paste0(path,'Co_bind_bindingsite_bed/',common_name,'_match_bindingsites')
  CpG.num.and.ratio.in.TFBS.logo_3200(PWM_file,bindingsite_file,Output_index, referdic,methyl.filename)
}

# summary_file="Result/Summary/Strip_TF_all_involved_fileID.csv"
# TF_involved <- read.csv(file = summary_file)
# all_files <- paste0(TF_involved[,"Biosample.term.name"],"_",TF_involved[,"TF"],"-human_",TF_involved[,"File.accession"],"_GRCh38")
# for(i in 1:length(all_files))
# {
#   temp <-all_files[i]
#   print(temp)
#   call_CpG.num.and.ratio.in.TFBS.logo (temp)
# }

all_files <- list.files(path = "Result/Homer/Co_bind/")
for(i in 1:length(all_files))
{
  temp <-all_files[i]
  print(temp)
  call_CpG.num.and.ratio.in.TFBS.logo_cobind(temp)
}