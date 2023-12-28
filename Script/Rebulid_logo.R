get_bindsite_methyl_3200 <- function(bindingsite_file,ratio_file,CpG_file,output,pading=3)
{
  library(Biostrings);
  methyl.Peak <- read.csv(bindingsite_file);
  colnames(methyl.Peak)[1:22] <- c("Chrom" ,"Start","End" ,"Name" ,"Score","Strand","Signalvalue","Pvalue" , "Qvalue","Peak","CpG.num" ,"methyl.CpG.num","na.CpG.num","max.read.num","max.methyl.num",
                                   "min.read.num","min.methyl.num","all.read.num","all.methyl.num","avg.nona.methyl.ratio", "avg.na.methyl.ratio","Seq")
  methyl.Peak <- methyl.Peak[order(methyl.Peak[,"Start"]), ];
  methyl.Peak <- methyl.Peak[order(methyl.Peak[,"Chrom"]), ];
  
  ratio <- read.csv(ratio_file,stringsAsFactors =F )
  CpG_num <- read.csv(CpG_file,stringsAsFactors =F )
  
  CpG.position.in.bindingsite <- rep("",nrow(ratio))
  CpG.methyl.in.bindingsite <- rep("",nrow(ratio))
  pading_first.sequence <- rep("",nrow(ratio))
  
  CpG.position.in.bindingsite_rev <- rep("",nrow(ratio))
  CpG.methyl.in.bindingsite_rev <- rep("",nrow(ratio))
  pading_first.sequence_rev <- rep("",nrow(ratio))
  
  for(i in 1:nrow(methyl.Peak))
  {
    if(methyl.Peak[i,"first.strand"]=="+")
    {
      pading_first.sequence[i] <- substr(as.character(methyl.Peak[i,"Seq"]),as.numeric(as.character(methyl.Peak[i,"first.position"]))-pading,as.numeric(as.character(methyl.Peak[i,"first.position"]))+nchar(as.character(methyl.Peak[i,"first.sequence"]))+pading-1)
      if(nchar(pading_first.sequence[i]) < nchar(as.character(methyl.Peak[i,"first.sequence"]))+2*pading )
      {
        pading_first.sequence[i] <- paste0(c(pading_first.sequence[i],rep("N",nchar(as.character(methyl.Peak[i,"first.sequence"]))+2*pading-nchar(pading_first.sequence[i]))),collapse = "")
      }
      index <- which(CpG_num[i,(1601-pading):(ncol(CpG_num)-1600+pading)]>0)
      if(length(index) >0 )
      {
        CpG.position.in.bindingsite[i] <- paste(index,collapse ="|")
        CpG.methyl.in.bindingsite[i] <- paste(ratio[i,index+1600-pading],collapse ="|")
      }
      #nutest
      pading_first.sequence_rev[i] <- as.character(reverseComplement(DNAStringSet(as.character(pading_first.sequence[i])))) 
      index <- which(CpG_num[i,((1601-pading)-1):(ncol(CpG_num)-1600+pading-1)]>0)
      if(length(index) >0 )
      {
        temp_ratio <-ratio[i,index+1600-pading-1]
        temp_index <- nchar(pading_first.sequence[i])-(index)+1 
        CpG.position.in.bindingsite_rev[i] <- paste(rev(temp_index),collapse ="|")
        CpG.methyl.in.bindingsite_rev[i] <- paste(rev(temp_ratio),collapse ="|")
      }
      
    }else
    {
      temp_seq <- substr(as.character(methyl.Peak[i,"Seq"]),as.numeric(as.character(methyl.Peak[i,"first.position"]))-pading,as.numeric(as.character(methyl.Peak[i,"first.position"]))+nchar(as.character(methyl.Peak[i,"first.sequence"]))+pading-1)
      pading_first.sequence[i] <- as.character(reverseComplement(DNAStringSet(as.character(temp_seq)))) 
      if(nchar(pading_first.sequence[i]) < nchar(as.character(methyl.Peak[i,"first.sequence"]))+2*pading )
      {
        pading_first.sequence[i] <- paste0(c(rep("N",nchar(as.character(methyl.Peak[i,"first.sequence"]))+2*pading-nchar(pading_first.sequence[i])),pading_first.sequence[i]),collapse = "")
      }
      index <- which(CpG_num[i,((1601-pading)-1):(ncol(CpG_num)-1600+pading-1)]>0)
      #index <- which(CpG_num[i,(((1601-pading)-1)-14):((1601-pading)-1)]>0)
      if(length(index) >0 )
      {
        temp_ratio <-ratio[i,index+1600-pading-1]
        temp_index <- (index)-1 
        CpG.position.in.bindingsite[i] <- paste(temp_index,collapse ="|")
        CpG.methyl.in.bindingsite[i] <- paste(temp_ratio,collapse ="|")
      }
      
      pading_first.sequence_rev[i] <- as.character(reverseComplement(DNAStringSet(as.character(pading_first.sequence[i])))) 
      index <- which(CpG_num[i,((1601-pading)-1):(ncol(CpG_num)-1600+pading-1)]>0)
      if(length(index) >0 )
      {
        temp_index <- nchar(pading_first.sequence[i])-index+1
        CpG.position.in.bindingsite_rev[i] <- paste(rev(temp_index) ,collapse ="|")
        CpG.methyl.in.bindingsite_rev[i] <- paste(rev(ratio[i,index+1600-pading-1]),collapse ="|")
      }
    }
  }
  temp <- cbind(methyl.Peak,CpG.methyl.in.bindingsite,pading_first.sequence,CpG.position.in.bindingsite,CpG.methyl.in.bindingsite,
                pading_first.sequence_rev,CpG.position.in.bindingsite_rev,CpG.methyl.in.bindingsite_rev)
  
  write.csv(temp,file = output,quote = F,row.names = F)
}
get5letters<-function(bindingsite_file_methyl="Result/Homer/High_0_6/H1-hESC_BACH1-human_ENCFF429LKK_GRCh38/H1-hESC_BACH1-human_ENCFF429LKK_GRCh38_match_bindingsites.methyl.csv")
{
  tempbsolute_bindingStart<- read.csv(bindingsite_file_methyl,stringsAsFactors = F)
  methylsite<-tempbsolute_bindingStart[,"CpG.position.in.bindingsite"]
  methyllevelsite<-tempbsolute_bindingStart[,"CpG.methyl.in.bindingsite.1"]
  sequence <- tempbsolute_bindingStart[,"pading_first.sequence"]
  for(i in 1:nrow(tempbsolute_bindingStart))
  {
    sites<-as.integer(unlist(strsplit(as.character(methylsite[i]),split="\\|")))
    methyllevel<-as.numeric(unlist(strsplit(as.character(methyllevelsite[i]),split="\\|")))
    tempsequence<-NULL
    if(length(sites)>0)
    {
      sites <- sites[which(methyllevel>=0.6)]
    }
    if(length(sites)>0)
    {
      if(length(sites)==1)
      {
        tempsequence<-paste0(tempsequence,substr(sequence[i],0,sites[1]-1))
        tempsequence<-paste0(tempsequence,"E")
        tempsequence<-paste0(tempsequence,substr(sequence[i],sites[1]+1,nchar(sequence[i])))
        sequence[i]<-tempsequence
        
      }
      if(length(sites)==2)
      {
        tempsequence<-paste0(tempsequence,substr(sequence[i],0,sites[1]-1))
        tempsequence<-paste0(tempsequence,"E")
        tempsequence<-paste0(tempsequence,substr(sequence[i],sites[1]+1,sites[2]-1))
        tempsequence<-paste0(tempsequence,"E")
        tempsequence<-paste0(tempsequence,substr(sequence[i],sites[2]+1,nchar(sequence[i])))
        sequence[i]<-tempsequence
      }
      if(length(sites)>=3)
      {
        tempsequence<-paste0(tempsequence,substr(sequence[i],0,sites[1]-1))
        tempsequence<-paste0(tempsequence,"E")
        for(j in 2:(length(sites)))
        {
          tempsequence<-paste0(tempsequence,substr(sequence[i],sites[j-1]+1,sites[j]-1))
          tempsequence<-paste0(tempsequence,"E")
        }
        tempsequence<-paste0(tempsequence,substr(sequence[i],sites[length(sites)]+1,nchar(sequence[i])))
        sequence[i]<-tempsequence
      }
    }
    
  }
  rm(methylsite)
  rm(methyllevelsite)
  
  methylsite_rev<-tempbsolute_bindingStart[,"CpG.position.in.bindingsite_rev"]
  methyllevelsite_rev<-tempbsolute_bindingStart[,"CpG.methyl.in.bindingsite_rev"]
  sequence_rev <- tempbsolute_bindingStart[,"pading_first.sequence_rev"]
  for(i in 1:nrow(tempbsolute_bindingStart))
  {
    sites<-as.integer(unlist(strsplit(as.character(methylsite_rev[i]),split="\\|")))
    methyllevel<-as.numeric(unlist(strsplit(as.character(methyllevelsite_rev[i]),split="\\|")))
    tempsequence<-NULL
    if(length(sites)>0)
    {
      sites <- sites[which(methyllevel>=0.6)]
    }
    
    
    if(length(sites)>0)
    {
      if(length(sites)==1)
      {
        tempsequence<-paste0(tempsequence,substr(sequence_rev[i],0,sites[1]-1))
        tempsequence<-paste0(tempsequence,"E")
        tempsequence<-paste0(tempsequence,substr(sequence_rev[i],sites[1]+1,nchar(sequence_rev[i])))
        sequence_rev[i]<-tempsequence
        
      }
      if(length(sites)==2)
      {
        tempsequence<-paste0(tempsequence,substr(sequence_rev[i],0,sites[1]-1))
        tempsequence<-paste0(tempsequence,"E")
        tempsequence<-paste0(tempsequence,substr(sequence_rev[i],sites[1]+1,sites[2]-1))
        tempsequence<-paste0(tempsequence,"E")
        tempsequence<-paste0(tempsequence,substr(sequence_rev[i],sites[2]+1,nchar(sequence_rev[i])))
        sequence_rev[i]<-tempsequence
      }
      if(length(sites)>=3)
      {
        tempsequence<-paste0(tempsequence,substr(sequence_rev[i],0,sites[1]-1))
        tempsequence<-paste0(tempsequence,"E")
        for(j in 2:(length(sites)))
        {
          tempsequence<-paste0(tempsequence,substr(sequence_rev[i],sites[j-1]+1,sites[j]-1))
          tempsequence<-paste0(tempsequence,"E")
        }
        tempsequence<-paste0(tempsequence,substr(sequence_rev[i],sites[length(sites)]+1,nchar(sequence_rev[i])))
        sequence_rev[i]<-tempsequence
      }
    }
    
  }
  return(list(sequence,sequence_rev))
}
getPWM5letters <-function(first.sequence,Rdata.name,png_out_index)
{
  library(ggseqlogo)
  library(ggplot2)
  t<-nchar(as.character(first.sequence[1]))
  Acount<-rep(0,t)
  Ccount<-rep(0,t)
  Gcount<-rep(0,t)
  Tcount<-rep(0,t)
  Ecount<-rep(0,t)
  for(i in 1:length(first.sequence))
  {
    tempchr=toupper(as.character(first.sequence[i]))
    for(j in 1:t)
    {
      if(substr(tempchr,j,j)=="A" )
      {
        Acount[j]<-Acount[j]+1
      }
      else
      {
        if(substr(tempchr,j,j)=="C")
        {
          Ccount[j]<-Ccount[j]+1
        }
        else
        {
          if(substr(tempchr,j,j)=="G")
          {
            Gcount[j]<-Gcount[j]+1
          }
          else
          {
            if(substr(tempchr,j,j)=="T")
            {
              Tcount[j]<-Tcount[j]+1
            }
            else
            {
              if(substr(tempchr,j,j)=="E")
              {
                Ecount[j]<-Ecount[j]+1
              }
            }
          }
        }
        
      }
    }
  }
  composite_motif_count <-cbind(Acount,Ccount,Gcount,Tcount,Ecount)
  
  motif_qurey  <- t(composite_motif_count/length(first.sequence))
  rownames(motif_qurey)<-c("A","C","G","T","E")
  csl2 <- make_col_scheme(chars = c("A","C","G","T","E"),group=c("1","2","3","4","5"),cols =c("#08d61d","#0d09ed","#edb409","#ed0909","#cc00ff"))
  gp <- ggplot()+geom_logo(motif_qurey,method = "prob",col_scheme = csl2,show_guide=FALSE)+labs(y="")+
    theme(panel.background = element_rect(fill = "transparent",color = NA))+theme(axis.text = element_blank())+theme(axis.ticks = element_blank())
  ggsave(paste0(png_out_index,".png"),gp,height = 4,width=20)
  
  # motif_qurey_rev_temp <- motif_qurey[c(4,3,2,1,5),]
  # motif_qurey_rev_temp[3,] <-  motif_qurey_rev_temp[3,]+motif_qurey_rev_temp[5,]
  # motif_qurey_rev_temp <- motif_qurey_rev_temp[,c(12:1)]
  # motif_qurey_rev_temp[5,] <- c(motif_qurey_rev_temp[5,c(2:12)],0)
  # motif_qurey_rev <- motif_qurey_rev_temp
  # rownames( motif_qurey_rev)<-c("A","C","G","T","E")
  # gp_rev <- ggplot()+geom_logo(motif_qurey_rev,method = "prob",col_scheme = csl2,show_guide=FALSE)+labs(y="")+
  #   theme(panel.background = element_rect(fill = "transparent",color = NA))+theme(axis.text = element_blank())+theme(axis.ticks = element_blank())
  # ggsave("2.png",gp_rev,height = 4,width=20)
  
  write.table(composite_motif_count,Rdata.name,col.names = T,row.names = F,sep = "\t",quote = F)
}
#for Strip TF
call_get_bindsite_methyl_3200 <- function(name)
{

  common_name<-name
  path<-'Result/Homer/'
  High_PWM_file<-paste0(path,'High_PWM/',common_name,'.Rdata')
  Low_PWM_file<-paste0(path,'Low_PWM/',common_name,'.Rdata')
  if(!file.exists(High_PWM_file))
  {
    print('exists')
    High_PWM_file<-Low_PWM_file
  }
  load(High_PWM_file)
  Cpg_index_high = paste0(path,'High_bindingsite_bed/',common_name,'_match_bindingsites',".motif.CpG.num.",3200+ncol(PWM),"bp.csv")
  ratio_index_high = paste0(path,'High_bindingsite_bed/',common_name,'_match_bindingsites',".motif.CpG.ratio.",3200+ncol(PWM),"bp.csv")
  
  load(Low_PWM_file)
  Cpg_index_low = paste0(path,'Low_bindingsite_bed/',common_name,'_match_bindingsites',".motif.CpG.num.",3200+ncol(PWM),"bp.csv")
  ratio_index_low = paste0(path,'Low_bindingsite_bed/',common_name,'_match_bindingsites',".motif.CpG.ratio.",3200+ncol(PWM),"bp.csv")
  
  bindingsite_file_high =paste0(path,'High_bindingsite_bed/',common_name,'_match_bindingsites.bed')
  bindingsite_file_low =paste0(path,'Low_bindingsite_bed/',common_name,'_match_bindingsites.bed')
  
  output_high <- paste0(path,'High_bindingsite_bed/',common_name,"_match_bindingsites.methyl.csv")
  output_low <- paste0(path,'Low_bindingsite_bed/',common_name,"_match_bindingsites.methyl.csv")
  get_bindsite_methyl_3200(bindingsite_file_high,ratio_index_high,Cpg_index_high,output_high)
  get_bindsite_methyl_3200(bindingsite_file_low,ratio_index_low,Cpg_index_low,output_low)
  
  first.sequence <- get5letters(output_high)
  Rdata.name <- paste0(path,'E_Logo_Rdata/',common_name,"_high_bindingsites.motif.E.csv")
  png_out_index <- paste0(path,"E_Logo_figures/",common_name,"_high_bindingsites.motif.E")
  Rdata.name_rev <- paste0(path,'E_Logo_Rdata/',common_name,"_high_bindingsites.motif.E_rev.csv")
  getPWM5letters(first.sequence[[1]],Rdata.name,png_out_index)
  getPWM5letters(first.sequence[[2]],Rdata.name_rev ,paste0(png_out_index,"_rev"))
  
  first.sequence <- get5letters(output_low)
  Rdata.name <- paste0(path,'E_Logo_Rdata/',common_name,"_low_bindingsites.motif.E.csv")
  png_out_index <- paste0(path,"E_Logo_figures/",common_name,"_low_bindingsites.motif.E")
  Rdata.name_rev <- paste0(path,'E_Logo_Rdata/',common_name,"_low_bindingsites.motif.E_rev.csv")
  getPWM5letters(first.sequence[[1]],Rdata.name,png_out_index)
  getPWM5letters(first.sequence[[2]],Rdata.name_rev ,paste0(png_out_index,"_rev"))
  
}
#for co-binding TF pairs
call_get_bindsite_methyl_3200_cobind <- function(name)
{
  common_name<-name
  path<-'Result/Homer/'
  PWM_file<-paste0(path,'Co_bind_PWM/',common_name,'.Rdata')
  load(PWM_file)
  Cpg_index = paste0(path,'Co_bind_bindingsite_bed/',common_name,'_match_bindingsites',".motif.CpG.num.",3200+ncol(PWM),"bp.csv")
  ratio_index = paste0(path,'Co_bind_bindingsite_bed/',common_name,'_match_bindingsites',".motif.CpG.ratio.",3200+ncol(PWM),"bp.csv")
  bindingsite_file =paste0(path,'Co_bind_bindingsite_bed/',common_name,'_match_bindingsites.bed')
  output <- paste0(path,'Co_bind_bindingsite_bed/',common_name,"_match_bindingsites.methyl.csv")
  get_bindsite_methyl_3200(bindingsite_file,ratio_index,Cpg_index,output,pading=12)

  first.sequence <- get5letters(output)
  Rdata.name <- paste0(path,'Co_bind_E_Logo_Rdata/',common_name,".motif.E.csv")
  png_out_index <- paste0(path,"Co_bind_E_Logo_figures/",common_name,".motif.E")
  Rdata.name_rev <- paste0(path,'Co_bind_E_Logo_Rdata/',common_name,".motif.E_rev.csv")
  getPWM5letters(first.sequence[[1]],Rdata.name,png_out_index)
  getPWM5letters(first.sequence[[2]],Rdata.name_rev ,paste0(png_out_index,"_rev"))
}

# summary_file="Result/Summary/Strip_TF_all_involved_fileID.csv"
# TF_involved <- read.csv(file = summary_file)
# all_files <- paste0(TF_involved[,"Biosample.term.name"],"_",TF_involved[,"TF"],"-human_",TF_involved[,"File.accession"],"_GRCh38")
# for(i in 1:length(all_files))
# {
#   temp <-all_files[i]
#   print(temp)
#   call_get_bindsite_methyl_3200(temp)
# }

# all_files <- list.files(path = "Result/Homer/Co_bind/")
# for(i in 1:length(all_files))
# {
#   temp <-all_files[i]
#   print(temp)
#   call_get_bindsite_methyl_3200_cobind(temp)
# }


Make_Logo_merge_Strip_TF_figure<- function()#to make the figure to publish the paper
{
  summary_file="Result/Summary/Strip_TF_all_involved_fileID.csv"
  TF_involved <- read.csv(file = summary_file)
  Strip_TF_new <- read.csv("~/Documents/RcodeLXM/Co_factor_methyl/Result/Summary/Strip_TF_new.csv")
  Sp_Low <- NULL;Sp_High <- NULL;Sp_Low_Cell <- NULL; Sp_High_Cell <- NULL
  for(i in 1:nrow(Strip_TF_new))
  {
    Sp_Low <- c(Sp_Low,unlist(strsplit(Strip_TF_new[i,"special_low"],"/")))
    Sp_Low_Cell <- c(Sp_Low_Cell,rep(Strip_TF_new[i,"X"],length(unlist(strsplit(Strip_TF_new[i,"special_low"],"/")))))
    Sp_High <- c(Sp_High,unlist(strsplit(Strip_TF_new[i,"special_hi"],"/")))
    Sp_High_Cell <- c(Sp_High_Cell,rep(Strip_TF_new[i,"X"],length(unlist(strsplit(Strip_TF_new[i,"special_hi"],"/")))))
  }
  TF_involved_Low <- TF_involved[match(paste0(Sp_Low,"_",Sp_Low_Cell),
                                     paste0(TF_involved[,"Unique_TF_name"] ,"_",TF_involved[,"Biosample.term.name"])),]
  TF_involved_High <- TF_involved[match(paste0(Sp_High,"_",Sp_High_Cell),
                                      paste0(TF_involved[,"Unique_TF_name"] ,"_",TF_involved[,"Biosample.term.name"])),]
  File_Low_File<- paste0(TF_involved_Low[,"Biosample.term.name"],"_",
                                  TF_involved_Low[,"TF"],"-human_",
                                  TF_involved_Low[,"File.accession"],"_GRCh38_low_bindingsites.motif.E.png")
  File_High_File<- paste0(TF_involved_High[,"Biosample.term.name"],"_",
                       TF_involved_High[,"TF"],"-human_",
                       TF_involved_High[,"File.accession"],"_GRCh38_high_bindingsites.motif.E.png")
  #make Logo merged figure for low
  library(png);library(grid);library(gridExtra)
  folder_path <- "Result/Homer/E_Logo_figures/"
  image_list <- paste0(folder_path, File_Low_File)#by hand
  for(i in c(9,10,12,13,14,16,21,22,25,29,31,32,33,35,36,40,42,43))
  {
    image_list[i] <- paste0(unlist(strsplit(image_list[i],split = "E.png"))[1],"E_rev.png")
  }
  new_index <- c(4,9:14,16,21:23,31,35:37,39,40,
                 1,3,5,17,18,24,41,6,7)
  All_index <- c(1:43)
  image_list<-image_list[c(new_index,All_index[!All_index %in% new_index])]
  image_list_high <- gsub("_low_bindingsites.motif", "_high_bindingsites.motif",image_list)
  for(i in c(1,11,28,31))#add E_rev
  {
    image_list_high[i] <- paste0(unlist(strsplit(image_list_high[i],split = "E.png"))[1],"E_rev.png")
  }
  for(i in c(5,9,10,32,36,42,43))#mv E_rev
  {
    image_list_high[i] <- paste0(unlist(strsplit(image_list_high[i],split = "E_rev.png"))[1],"E.png")
  }
  #First merge with in the TF, then all TFs to one PNG figure
  for (i in 1:43)
    {
    temp_merge <-list()
    temp_img_low <- readPNG(image_list[i])
    temp_img_high <- readPNG(image_list_high[i])
    raster_img <- rasterGrob(temp_img_low, interpolate = TRUE)
    temp_merge<- c(temp_merge,list(raster_img))
    raster_img <- rasterGrob(temp_img_high, interpolate = TRUE)
    temp_merge<- c(temp_merge,list(raster_img))
    temp_name <- gsub("Result/Homer/E_Logo_figures/", "Result/Homer/Logo_merge/",image_list[i])
    temp_name <- paste0(unlist(strsplit(temp_name,split = "_GRCh38"))[1],".png")
    ggsave(temp_name,width=2, height=1,
         marrangeGrob(grobs = temp_merge,nrow = 2,ncol=1,top=NULL))
  }
  combined_plot <- list()
  for (i in 1:43) {
    temp_merge <-list()
    temp_name <- gsub("Result/Homer/E_Logo_figures/", "Result/Homer/Logo_merge/",image_list[i])
    temp_name <- paste0(unlist(strsplit(temp_name,split = "_GRCh38"))[1],".png")
    img<- readPNG(temp_name)
    raster_img <- rasterGrob(img, interpolate = TRUE)
    combined_plot <- c(combined_plot, list(raster_img))
  }
  library(ggplot2)
  ggsave("Result/Homer/Logo_merge/Stripe_Low_TF.pdf",width=8, height=14, marrangeGrob(grobs = combined_plot,nrow = 11,ncol=4,top=NULL))
  ggsave("Result/Homer/Logo_merge/Stripe_Low_TF.png",width=8, height=14, marrangeGrob(grobs = combined_plot,nrow = 11,ncol=4,top=NULL))
  image_list <- lapply(image_list,function(x){
    temp <- strsplit(x,split = "/")[[1]][4]
    temp <-strsplit(temp ,split = "_GRCh38")[[1]][1]
    temp <- gsub("-human", "",temp)
  })
  write.csv(unlist(image_list),file = "Result/Homer/Logo_merge/Stripe_Low_TF.csv")

  #make Logo merged figure for Hi
  library(png);library(grid);library(gridExtra);library(ggplot2)
  folder_path <- "Result/Homer/E_Logo_figures/"
  image_list <- paste0(folder_path, File_High_File)#by hand
  image_list_low <- gsub("_high_bindingsites.motif", "_low_bindingsites.motif",image_list)
  for(i in c(1,5))#add E_rev
  {
    image_list[i] <- paste0(unlist(strsplit(image_list[i],split = "E.png"))[1],"E_rev.png")
  }
  for(i in c(2,4))#add E_rev
  {
    image_list_low[i] <- paste0(unlist(strsplit(image_list_low[i],split = "E.png"))[1],"E_rev.png")
  }
  combined_plot <- list()
  for (i in 1:5) {
  temp_merge <-list()
  temp_img_low <-readPNG(image_list_low[i]) 
  temp_img_high <- readPNG(image_list[i])
  raster_img <- rasterGrob(temp_img_low, interpolate = TRUE)
  temp_merge<- c(temp_merge,list(raster_img))
  raster_img <- rasterGrob(temp_img_high, interpolate = TRUE)
  temp_merge<- c(temp_merge,list(raster_img))
  temp_name <- gsub("Result/Homer/E_Logo_figures/", "Result/Homer/Logo_merge/",image_list[i])
  temp_name <- paste0(unlist(strsplit(temp_name,split = "_GRCh38"))[1],".png")
  ggsave(temp_name,width=2, height=1,
         marrangeGrob(grobs = temp_merge,nrow = 2,ncol=1,top=NULL))
  temp_img <- readPNG(temp_name)
  raster_img <- rasterGrob(temp_img, interpolate = TRUE)
  combined_plot <- c(combined_plot, list(raster_img))
}
  ggsave("Result/Homer/Logo_merge/Stripe_High_TF.pdf",width=8, height=2.5, marrangeGrob(grobs = combined_plot,nrow = 2,ncol=4,top=NULL))
  ggsave("Result/Homer/Logo_merge/Stripe_High_TF.png",width=8, height=2.5, marrangeGrob(grobs = combined_plot,nrow = 2,ncol=4,top=NULL))
  image_list <- lapply(image_list,function(x){
    temp <- strsplit(x,split = "/")[[1]][4]
    temp <-strsplit(temp ,split = "_GRCh38")[[1]][1]
    temp <- gsub("-human", "",temp)
  })
  write.csv(unlist(image_list),file = "Result/Homer/Logo_merge/Stripe_High_TF.csv")
}
Make_Logo_merge_Co_bind_figure<- function()#to make the figure to publish the paper
{
  High_Low <- read.csv(file = "Result/Pair_group/ASummary_High_Low.csv",stringsAsFactors = F)
  High_Low <- High_Low[High_Low[,"cell_count"]>=2,]
  image_list <- NULL
  for(i in 1:nrow(High_Low ))
  {
    temp <- High_Low[i,]
    TF_1 <-  unlist(strsplit( temp[,"index"]," _ "))[1]
    TF_2 <-  unlist(strsplit( temp[,"index"]," _ "))[2]
    for(cell in unlist(strsplit( temp[,"cell_list"],"_")))
    {
      image_list <- c(image_list,paste0(cell,"_",TF_1,"_",TF_2,".motif.E.png"))
    }
  }
  #make Logo merged figure for High
  library(png);library(grid);library(gridExtra)
  folder_path <- "Result/Homer/Co_bind_E_Logo_figures/"
  image_list <- paste0(folder_path, image_list)#by hand
  temp_merge <-list()
  for (i in 1:30)
  {
    temp_img <- readPNG(image_list[i])
    raster_img <- rasterGrob(temp_img, interpolate = TRUE)
    temp_merge<- c(temp_merge,list(raster_img))
  }
  ggsave("Result/Homer/Co_High_Low.png",width=2, height=5,
           marrangeGrob(grobs = temp_merge,nrow = 15,ncol=2,top=NULL))
  #######################
  High_Low <- read.csv(file = "Result/Pair_group/ASummary_Low_High.csv",stringsAsFactors = F)
  High_Low <- High_Low[High_Low[,"cell_count"]>=2,]
  image_list <- NULL
  for(i in 1:nrow(High_Low ))
  {
    temp <- High_Low[i,]
    TF_1 <-  unlist(strsplit( temp[,"index"]," _ "))[1]
    TF_2 <-  unlist(strsplit( temp[,"index"]," _ "))[2]
    for(cell in unlist(strsplit( temp[,"cell_list"],"_")))
    {
      image_list <- c(image_list,paste0(cell,"_",TF_1,"_",TF_2,".motif.E.png"))
    }
  }
  #make Logo merged figure for High
  library(png);library(grid);library(gridExtra)
  folder_path <- "Result/Homer/Co_bind_E_Logo_figures/"
  image_list <- paste0(folder_path, image_list)#by hand
  temp_merge <-list()
  for (i in 1:45)
  {
    temp_img <- readPNG(image_list[i])
    raster_img <- rasterGrob(temp_img, interpolate = TRUE)
    temp_merge<- c(temp_merge,list(raster_img))
  }
  ggsave("Result/Homer/Co_Low_High.png",width=2, height=5,
         marrangeGrob(grobs = temp_merge,nrow = 15,ncol=2,top=NULL))
  
}
