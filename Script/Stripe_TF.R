Identify_Specific_Strip_TF<- function(cell="HeLa-S3")
{
  library(ggplot2)
  library(reshape2)
  library(pheatmap)
  High<- read.csv(file = paste0("Result/Summary/",cell,"_high_overlap_count_to_0_1.csv"),stringsAsFactors = F,row.names = 1)
  Low<- read.csv(file = paste0("Result/Summary/",cell,"_low_overlap_count_to_0_1.csv"),stringsAsFactors = F,row.names = 1)
  
  over_High_O_1 <- High
  row_sum <- apply(over_High_O_1, 1, sum)
  length(which(row_sum>0))
  High_col_ratio  <- apply(over_High_O_1, 2, sum)/length(which(row_sum>0))

  over_Low_O_1 <- Low
  row_sum <- apply(over_Low_O_1, 1, sum)
  Low_col_ratio <- apply(over_Low_O_1, 2, sum)/length(which(row_sum>0))
  
  all_index <-unique(c(names( High_col_ratio), names( Low_col_ratio)))
  t <- cbind(all_index,Low = Low_col_ratio[match(all_index,names( Low_col_ratio))],High = High_col_ratio[match(all_index,names( High_col_ratio))])
  t <- t[which(t[,2]>0|t[,3]>0),]
  Low_Stripe_TF <- t[which(t[,2] > 0.5),1]
  High_Stripe_TF <- t[which(t[,3] > 0.5),1]
  Special_low <- t[which(t[,2] > 0.5 & as.numeric(t[,2])/as.numeric(t[,3]) > 2),1]
  Special_hi <- t[which(t[,3] > 0.5 & as.numeric(t[,3])/as.numeric(t[,2]) > 2),1]

  library(ggplot2)
  type = rep("Com",nrow(t))
  type[which(t[,2] > 0.5 & as.numeric(t[,2])/as.numeric(t[,3]) > 2)] ="Spe_Low"
  type[which(t[,3] > 0.5 & as.numeric(t[,3])/as.numeric(t[,2]) > 2)] = "Spe_Hi"
  data <- data.frame(Low = as.numeric(t[,2]),High = as.numeric(t[,3]),Type = type)
  
  # p <- ggplot(data) + geom_point(aes(x = Low, y = High,color= Type ))+theme_bw() + theme(legend.position =  "bottom",text = element_text(color = "black")) +   
  # scale_color_manual(values = c("#FEB048", "#ed0e0e", "#9570D3"))+xlim(0,1)+ylim(0,1)+expand_limits(x=0,y=0)
  p <- ggplot(data) + geom_point(aes(x = Low, y = High,color= Type ))+theme_bw() + theme(legend.position ="bottom",text = element_text(color = "black")) +   
    scale_color_manual(values = c("#FEB048",  "#9570D3","#ed0e0e"))+xlim(0,1)+ylim(0,1)+expand_limits(x=0,y=0)
  
  ggsave(paste0("Result/Summary/Strip_TF_",cell,"_plot.pdf"),p,dpi=300,width = 2.5,height = 2.5)
  
  annotation.fit_add_names <- read.csv("data/annotation/annotation.fit_add_names.csv")
  TF_file_name <-  annotation.fit_add_names[which( annotation.fit_add_names[,"Biosample.term.name"]==cell),c("File.accession","Biosample.term.name","Unique_TF_name")]
  TF_candicate <- c(High_Stripe_TF,Low_Stripe_TF)
  TF_involved <- TF_file_name[TF_file_name[,"Unique_TF_name"] %in%TF_candicate, ]

  return(list(High_Stripe_TF,Low_Stripe_TF,Special_hi,Special_low,TF_involved))
  
}

get_summmary_All_STF <- function()
{
  All_TF_involved  <- NULL
  All_STF <- matrix(data = NA,nrow = 9,ncol = 7)
  rownames(All_STF) <- c("GM12878","H1","HeLa-S3","K562","HepG2","A549","HEK293","SK-N-SH","MCF-7")
  for(cell in c("GM12878","H1","HeLa-S3","K562","HepG2","A549","HEK293","SK-N-SH","MCF-7"))
  {
  temp_Strip_TF<- Identify_Specific_Strip_TF(cell)
  
  temp_low_raw_STF<- paste0(temp_Strip_TF[[2]],collapse = "/")
  temp_low_unique_STF<- unique(unlist(lapply(temp_Strip_TF[[2]], function(x){
    unlist(strsplit(x,"_"))[1]
  })))
  
  temp_hi_raw_STF<- paste0(temp_Strip_TF[[1]],collapse = "/")
  temp_hi_unique_STF<- unique(unlist(lapply(temp_Strip_TF[[1]], function(x){
    unlist(strsplit(x,"_"))[1]
  })))
  
  All_STF[cell,1] <- temp_low_raw_STF
  All_STF[cell,2] <- paste0(temp_low_unique_STF,collapse = "/")
  All_STF[cell,3] <- temp_hi_raw_STF
  All_STF[cell,4] <- paste0(temp_hi_unique_STF,collapse = "/")
  All_STF[cell,5] <- paste0(intersect(temp_low_unique_STF,temp_hi_unique_STF),collapse = "/")
  All_STF[cell,6] <- paste0(temp_Strip_TF[[4]],collapse = "/")
  All_STF[cell,7] <- paste0(temp_Strip_TF[[3]],collapse = "/")
  All_TF_involved <- rbind(All_TF_involved,temp_Strip_TF[[5]])
  }
  
  colnames(All_STF) <- c("low_raw_STF","low_unique_STF","hi_raw_STF","hi_unique_STF","overlap","special_low","special_hi")
   write.csv(All_STF,file = "Result/Summary/Strip_TF_new.csv",quote = F)
   write.csv(All_TF_involved,file = "Result/Summary/Strip_TF_all_involved_fileID.csv",quote = F)
}
#get_summmary_All_STF()

Go_term_Gene <- function(tfbs.seq.filename="GM12878_IKZF1-human_ENCFF819VMH_GRCh38.csv")
{
  Peaks.seq <- read.csv(paste0("Result/TFBS.Peak.methyl/",tfbs.seq.filename),stringsAsFactors = F);
  high.methyl.Peak <- Peaks.seq[Peaks.seq[,"avg.na.methyl.ratio"] >= 0.6 & !is.na(Peaks.seq[,"avg.na.methyl.ratio"]), ];
  high_regions <- GRanges(seqnames = high.methyl.Peak[,1],IRanges(start = high.methyl.Peak[,2],end = high.methyl.Peak[,3]))
  
  
  library("TxDb.Hsapiens.UCSC.hg38.knownGene")
  genes <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(biomaRt)
  library("GenomicFeatures")
  library(clusterProfiler)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  bm <- getBM(attributes = c("external_gene_name",'entrezgene_id'), values=names(genes),filters ='entrezgene_id', mart = mart)
  names(genes) <- bm$external_gene_name[match(genes$gene_id,bm$entrezgene_id)]
　hg38_gene_promoter<-promoters(genes,downstream=0)
  overlap_index  <- findOverlaps(hg38_gene_promoter,high_regions)
  Overlap_gene <- hg38_gene_promoter[overlap_index@from]
  

  library("org.Hs.eg.db")
  library("topGO")
  library(Rgraphviz)
  temp <- enrichGO(gene     = unlist(Overlap_gene@elementMetadata@listData),
           OrgDb    = org.Hs.eg.db,
           ont="BP",
           pvalueCutoff=1, 
           qvalueCutoff=1)
  BP_result <- temp@result
  pdf(paste0("Result/Go_Term/High_",tfbs.seq.filename,".BP.pdf"))
  plotGOgraph(temp)
  dev.off()
  
  temp <- enrichGO(gene     = unlist(Overlap_gene@elementMetadata@listData),
                   OrgDb    = org.Hs.eg.db,
                   ont="MF",
                   pvalueCutoff=1, 
                   qvalueCutoff=1)
  MF_result <- temp@result
  pdf(paste0("Result/Go_Term/High_",tfbs.seq.filename,".MF.pdf"))
  plotGOgraph(temp)
  dev.off()
  
  temp <- enrichGO(gene     = unlist(Overlap_gene@elementMetadata@listData),
                   OrgDb    = org.Hs.eg.db,
                   ont="CC",
                   pvalueCutoff=1, 
                   qvalueCutoff=1)
  CC_result <- temp@result
  pdf(paste0("Result/Go_Term/High_",tfbs.seq.filename,".CC.pdf"))
  plotGOgraph(temp)
  dev.off()
  
  Go_result <- rbind(BP_result[1:10,],MF_result[1:10,],CC_result[1:10,] )
  write.csv(Go_result,file = paste0("Result/Go_Term/High_",tfbs.seq.filename),quote = F,row.names = F)
}

# Go_term_Gene("GM12878_IKZF1-human_ENCFF819VMH_GRCh38.csv")
# Go_term_Gene("H1_CTCF-human_ENCFF692RPA_GRCh38.csv")
# Go_term_Gene("H1_ZNF143-human_ENCFF235ROG_GRCh38.csv")
# Go_term_Gene("K562_ZBTB33-human_ENCFF917RIN_GRCh38.csv")
# Go_term_Gene("HepG2_ZNF687-human_ENCFF581GZR_GRCh38.csv")
# Go_term_Gene("SK-N-SH_GATA3-human_ENCFF833ZAC_GRCh38.csv")
# all_files<- list.files("Result/Go_Term/")
# all_Go <- NULL
# for(i in 1:6)
# {
#   all_Go<- c(all_Go,read.csv(paste0("Result/Go_Term/",all_files[i]))[,1])
# }
# t <- unique(all_Go)

Strip_TF_feature <-function()
{
  Homo_sapiens_TF <- read.delim("~/Documents/RcodeLXM/Co_factor_methyl/data/annotation/Homo_sapiens_TF.txt")
  Homo_sapiens_TF_cofactors <- read.delim("~/Documents/RcodeLXM/Co_factor_methyl/data/annotation/Homo_sapiens_TF_cofactors.txt")

  Strip_TF_new <- read.csv("~/Documents/RcodeLXM/Co_factor_methyl/Result/Summary/Strip_TF_new.csv")
  Sp_Low <- NULL
  Sp_High <- NULL
  for(i in 1:nrow(Strip_TF_new))
  {
    Sp_Low <- c(Sp_Low,unlist(strsplit(Strip_TF_new[i,"special_low"],"/")))
    Sp_High <- c(Sp_High,unlist(strsplit(Strip_TF_new[i,"special_hi"],"/")))
  }
Sp_Low  <- unlist(lapply(Sp_Low, function(x){unlist(strsplit(x,"_"))[1]}))#43
Sp_High  <- unlist(lapply(Sp_High, function(x){unlist(strsplit(x,"_"))[1]}))#5
Sp_Low <- unique(Sp_Low)#34
Sp_High  <- unique(Sp_High)#4
Sp_Low_Family <-  Homo_sapiens_TF [ Homo_sapiens_TF [,"Symbol"]%in%Sp_Low, c("Symbol","Family")]#only 22 ,lost 10
Sp_Low_Family <- rbind(Sp_Low_Family ,Homo_sapiens_TF_cofactors [ Homo_sapiens_TF_cofactors [,"Symbol"]%in%Sp_Low, c("Symbol","Family")])
Sp_High_Family <-  Homo_sapiens_TF [ Homo_sapiens_TF [,"Symbol"]%in%Sp_High, c("Symbol","Family")]#only 22 ,lost 10
 Homo_sapiens_TF_cofactors [ Homo_sapiens_TF_cofactors [,"Symbol"]=="HNF4A", c("Symbol","Family")]#only 22 ,lost 10

summary(as.factor(Sp_Low_Family [,2]))
summary(as.factor(Sp_High_Family[,2]))
ancerdrivers <- read.delim("data/annotation/NCG_cancerdrivers_systemslevelproperties.tsv")
ancerdrivers[ancerdrivers[,"symbol"] %in% Sp_Low,]#14
ancerdrivers[ancerdrivers[,"symbol"] %in% Sp_High,]#1

healthydriver <- read.delim("data/annotation/NCG_healthydrivers_annotation_supporting_evidence.tsv")
length(which(Sp_Low %in% healthydriver[,"symbol"]))#2
length(which(Sp_High%in% healthydriver[,"symbol"]))#0
healthydriver[healthydriver[,"symbol"] %in% Sp_Low,]#2
}






