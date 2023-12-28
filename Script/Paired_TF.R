duplicated_TFpair <- function(TF1,TF2)
{
  TF_index <- c(TF1,TF2)
  TF_number <- unique(TF_index)
  TF_1 <-  match(TF1,TF_number)
  TF_2 <-  match(TF2,TF_number)
  all_specific_index <- NULL
  for(i in 1:length(TF1))
  {
    if(TF_1[i] > TF_2[i])
    {
      all_specific_index[i] <- paste(TF1[i] ,"_",TF2[i])
    }
    else
    {
      all_specific_index[i] <- paste(TF2[i] ,"_",TF1[i])
    }
  }
  return(duplicated(all_specific_index))
}
# minp is set based on the figures published in the paper. 
# If the TFs are too many and the interactions are too dense. We set the minp higher.
get_pair_and_module <- function(cell,minp=3,vertex_label_cex=1.3,vertex_size=25)
{
  over_High_O_1<- read.csv(file = paste0("Result/Summary/",cell,"_high_overlap_count_to_0_1.csv"),stringsAsFactors = F,row.names = 1)
  over_High_O_1_cosine <- read.csv(file = paste0("Result/Summary/",cell,"_high_overlap_count_to_0_1_consine.csv"),stringsAsFactors = F,row.names = 1)
  over_Low_O_1<- read.csv(file = paste0("Result/Summary/",cell,"_low_overlap_count_to_0_1.csv"),stringsAsFactors = F,row.names = 1)
  over_Low_O_1_cosine <- read.csv(file = paste0("Result/Summary/",cell,"_low_overlap_count_to_0_1_consine.csv"),stringsAsFactors = F,row.names = 1)

  for(i in 1:nrow(over_High_O_1))
  {
    for(j in 1:ncol(over_High_O_1))
    {
      if(over_High_O_1[i,j] ==1 &over_High_O_1_cosine[i,j]==1)
      {
        over_High_O_1[i,j] =1
      }else
      {
        over_High_O_1[i,j]=0
      }
    }
  }
  for(i in 1:nrow(over_Low_O_1))
  {
    for(j in 1:ncol(over_Low_O_1))
    {
      if(over_Low_O_1[i,j] ==1 &over_Low_O_1_cosine[i,j]==1)
      {
        over_Low_O_1[i,j] =1
      }else
      {
        over_Low_O_1[i,j]=0
      }
    }
  }
  print(paste0(cell," has ",sum(over_High_O_1)," high-methylated interactions and ",sum(over_Low_O_1)," unmethylated interactions!"))
  TF1 <- NULL
  TF2 <- NULL
  for(i in 1:nrow(over_High_O_1))
  {
    for(j in i:ncol(over_High_O_1)){
    if(over_High_O_1[i,j]==1 & over_High_O_1[j,i]==1 ){
      r1 <-unique(unlist(lapply(row.names(over_High_O_1)[i], function(x){
        unlist(strsplit(x,"_"))[1]
      })))
      r2 <-unique(unlist(lapply(colnames(over_High_O_1)[j], function(x){
        unlist(strsplit(x,"_"))[1]
      })))
      if(r1!=r2)
      {
        TF1 <- c(TF1,r1)
        TF2 <- c(TF2,r2)
      }
    }}}
  index <- duplicated_TFpair(TF1,TF2)
  TF1 <-  TF1[!index]
  TF2 <-  TF2[!index]
  write.csv(data.frame(TF1,TF2),file = paste0("Result/Pair_group/High_",cell,".paired.csv"),quote = F,row.names = F)

  library("igraph")
  net_PC <- graph_from_data_frame(data.frame(TF1,TF2))
  hist(igraph::degree(net_PC),col="lightblue",xlim=c(0,50),
       xlab="Vertex Degree",ylab="frequency",
       main="Degree Distribution")
  Most_connected_subgraph <- max_cliques(net_PC, min=minp)#find Most connected subgraph
  
  subtable_point<- unique(names(unlist(Most_connected_subgraph)))
  subtable <- data.frame(TF1,TF2)[which(TF1 %in% subtable_point & TF2 %in% subtable_point),]
  Most_connected <- graph_from_data_frame(subtable)
  pdf(file = paste0("Result/Pair_group/High_",cell,".Most_connected_subgraph.pdf"),width =6,height = 6)
  par(mfrow=c(1,1), mar=c(0.1,0.1,0.1,0.1)) 
  l <- layout_components(Most_connected) 
  plot.igraph(Most_connected,edge.arrow.size=.4, edge.curved=.0,edge.color="gray50",vertex.color="orange",
              vertex.label.cex=vertex_label_cex,vertex.label.color="black",vertex.label.dist=0,vertex.size=vertex_size,layout=l,edge.arrow.mode=0)
  dev.off()
  save(Most_connected_subgraph,file = paste0("Result/Pair_group/High_",cell,".Most_connected_subgraph.Rdata"))
  

  TF1 <- NULL
  TF2 <- NULL
  for(i in 1:nrow(over_Low_O_1))
  {for(j in i:ncol(over_Low_O_1)){
    if(over_Low_O_1[i,j]==1 & over_Low_O_1[j,i]==1 ){
      r1 <-unique(unlist(lapply(row.names(over_Low_O_1)[i], function(x){
        unlist(strsplit(x,"_"))[1]
      })))
      r2 <-unique(unlist(lapply(colnames(over_Low_O_1)[j], function(x){
        unlist(strsplit(x,"_"))[1]
      })))
      if(r1!=r2)
      {
        TF1 <- c(TF1,r1)
        TF2 <- c(TF2,r2)
      }
    }}}
  
  index <- duplicated_TFpair(TF1,TF2)
  TF1 <-  TF1[!index]
  TF2 <-  TF2[!index]
  write.csv(data.frame(TF1,TF2),file = paste0("Result/Pair_group/Low_",cell,".paired.csv"),quote = F,row.names = F)
  
  library("igraph")
  net_PC <- graph_from_data_frame(data.frame(TF1,TF2))
  hist(igraph::degree(net_PC),col="lightblue",xlim=c(0,50),
       xlab="Vertex Degree",ylab="frequency",
       main="Degree Distribution")
  Most_connected_subgraph <- max_cliques(net_PC, min=minp)#find Most connected subgraph
  
  subtable_point<- unique(names(unlist(Most_connected_subgraph)))
  subtable <- data.frame(TF1,TF2)[which(TF1 %in% subtable_point & TF2 %in% subtable_point),]
  Most_connected <- graph_from_data_frame(subtable)
  pdf(file = paste0("Result/Pair_group/Low_",cell,".Most_connected_subgraph.pdf"),width = 6,height = 6)
  par(mfrow=c(1,1), mar=c(0.1,0.1,0.1,0.1))
  l <- layout_components(Most_connected) 
  plot.igraph(Most_connected,edge.arrow.size=.4, edge.curved=.0,edge.color="gray50",vertex.color="lightblue",
              vertex.label.cex=vertex_label_cex,vertex.label.color="black",vertex.label.dist=0,vertex.size=vertex_size,layout=l,edge.arrow.mode=0)
  dev.off()
  save(Most_connected_subgraph,file = paste0("Result/Pair_group/Low_",cell,".Most_connected_subgraph.Rdata"))
}

summary_count_High_Low_and_specific <- function()
{
  Summ_table = data.frame(0,nrow=9,ncol = 7)
  i<-1
  for(cell in c("GM12878","H1","HeLa-S3","K562","HepG2","A549","HEK293","SK-N-SH","MCF-7"))
  {
    over_High_O_1<- read.csv(file = paste0("Result/Summary/",cell,"_high_overlap_count_to_0_1.csv"),stringsAsFactors = F,row.names = 1)
    over_Low_O_1<- read.csv(file = paste0("Result/Summary/",cell,"_low_overlap_count_to_0_1.csv"),stringsAsFactors = F,row.names = 1)
    High_TF <-unique(unlist(lapply(colnames(over_High_O_1), function(x){
      unlist(strsplit(x,"_"))[1]
    })))
    Low_TF <-unique(unlist(lapply(colnames(over_Low_O_1), function(x){
      unlist(strsplit(x,"_"))[1]
    })))
    
    High_pair <- read.csv(paste0("Result/Pair_group/High_",cell,".paired.csv"))
    Low_pair <- read.csv(paste0("Result/Pair_group/Low_",cell,".paired.csv"))
    High_pair_S <- read.csv(paste0("Result/Pair_group/High_",cell,".paired_specific.csv"))
    Low_pair_S <- read.csv(paste0("Result/Pair_group/Low_",cell,".paired_specific.csv"))
    Summ_table[i,1] <- cell
    Summ_table[i,2] <- length(High_TF)
    Summ_table[i,3] <- length(Low_TF)
    Summ_table[i,4] <- nrow(High_pair)
    Summ_table[i,5] <- nrow(Low_pair)
    Summ_table[i,6] <- nrow(High_pair_S)
    Summ_table[i,7] <- nrow(Low_pair_S)
    i<-i+1
  }
  colnames(Summ_table) <- c("cell","Count_Me_TF","Count_Un-Me_TF","Methylation","Un-Me","Specific Me","Specific Un-me")
}
compare_High_Low <- function(cell)
{
  High_pair <- read.csv(paste0("Result/Pair_group/High_",cell,".paired.csv"))
  Low_pair <- read.csv(paste0("Result/Pair_group/Low_",cell,".paired.csv"))
  High_pair_index <- paste0(High_pair[,1],"_",High_pair [,2])
  Low_pair_index <- c(paste0( Low_pair[,1],"_", Low_pair [,2]),paste0( Low_pair[,2],"_", Low_pair [,1]))
  write.csv(High_pair[which(!High_pair_index %in% Low_pair_index),],file = paste0("Result/Pair_group/High_",cell,".paired_specific.csv"),quote = F)
  
  Low_pair_index <- paste0( Low_pair[,1],"_", Low_pair [,2])
  High_pair_index <- c(paste0(High_pair[,1],"_",High_pair [,2]),paste0(High_pair[,2],"_",High_pair [,1]))
  write.csv(Low_pair[which(!Low_pair_index %in% High_pair_index),],file = paste0("Result/Pair_group/Low_",cell,".paired_specific.csv"),quote = F)
  return(c(nrow(High_pair),nrow(Low_pair)))
}
summary_compare_High_Low <- function()
{
  all_specific <- NULL
  for(cell in c("GM12878","H1","HeLa-S3","K562","HepG2","A549","HEK293","SK-N-SH","MCF-7"))
  {
    temp <-read.csv( file = paste0("Result/Pair_group/High_",cell,".paired_specific.csv"))
    all_specific<- rbind(all_specific, cbind(temp,rep(cell,nrow(temp))))
  }
  TF_index <- c(all_specific[,2],all_specific[,3])
  TF_number <- unique(TF_index)
  TF_1 <-  match(all_specific[,2],TF_number)
  TF_2 <-  match(all_specific[,3],TF_number)
  all_specific_index <- NULL
  for(i in 1:length(all_specific[,2]))
  {
    if(TF_1[i] > TF_2[i])
    {
      all_specific_index[i] <- paste(all_specific[i,2] ,"_",all_specific[i,3])
    }
    else
    {
      all_specific_index[i] <- paste(all_specific[i,3] ,"_",all_specific[i,2])
    }
  }
  cell_list <- NULL 
  cell_count <- NULL
  index <-unique(all_specific_index)
  for(i in 1:length(index))
  {
    cell_count[i] <- length(which(index[i]==all_specific_index))
    cell_list[i] <- paste(all_specific[which(all_specific_index==index[i]),4],collapse ="_" )
  }
  cell_list[cell_count > 1]
  temp <- cbind(index,cell_list,cell_count)
  write.csv(temp,file = "Result/Pair_group/ASummary_High_Low.csv")
}
summary_compare_Low_High <- function()
{
  all_specific <- NULL
  for(cell in c("GM12878","H1","HeLa-S3","K562","HepG2","A549","HEK293","SK-N-SH","MCF-7"))
  {
    temp <-read.csv( file = paste0("Result/Pair_group/Low_",cell,".paired_specific.csv"))
    all_specific<- rbind(all_specific, cbind(temp,rep(cell,nrow(temp))))
  }
  TF_index <- c(all_specific[,2],all_specific[,3])
  TF_number <- unique(TF_index)
  TF_1 <-  match(all_specific[,2],TF_number)
  TF_2 <-  match(all_specific[,3],TF_number)
  all_specific_index <- NULL
  for(i in 1:length(all_specific[,2]))
  {
    if(TF_1[i] > TF_2[i])
    {
      all_specific_index[i] <- paste(all_specific[i,2] ,"_",all_specific[i,3])
    }
    else
    {
      all_specific_index[i] <- paste(all_specific[i,3] ,"_",all_specific[i,2])
    }
  }
  cell_list <- NULL 
  cell_count <- NULL
  index <-unique(all_specific_index)
  for(i in 1:length(index))
  {
    cell_count[i] <- length(which(index[i]==all_specific_index))
    cell_list[i] <- paste(all_specific[which(all_specific_index==index[i]),4],collapse ="_" )
  }
  cell_list[cell_count > 1]
  temp <- cbind(index,cell_list,cell_count)
  write.csv(temp,file = "Result/Pair_group/ASummary_Low_High.csv")
  
}
call_function <- function()
{
  all_cell <- c("GM12878","A549","HEK293","SK-N-SH","MCF-7")
  for(i in 1:5)
  {
    cell <- all_cell[i]
    get_pair_group(cell,3)
  }
  get_pair_group("HeLa-S3",2)
  get_pair_group("H1",2)
  get_pair_group("K562",6,0.8,16)
  get_pair_group("HepG2",5,0.8,16)
  count_TF_pair <- matrix(0,nrow = 9,ncol = 2)#c(nrow(High_pair),nrow(Low_pair))
  all_cell <- c("GM12878","H1","HeLa-S3","K562","HepG2","A549","HEK293","SK-N-SH","MCF-7")
  for(i in 1:9 )
  {
    cell <- all_cell[i]
    temp <- compare_High_Low(cell)
    count_TF_pair[i,1]<-temp[1]
    count_TF_pair[i,2]<-temp[2]
  }
  row.names(count_TF_pair) <- all_cell
  colnames(count_TF_pair) <- c("High","Low")
  summary_compare_High_Low()
  summary_compare_Low_High()
}
#call_function()