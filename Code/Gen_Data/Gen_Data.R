

Gen_Data<-function(data_name,samples)
{

  library(bnlearn)

  
  random_seed<-123    # random seed
  # If you don¡¯t change it, the generated data will be the same every time
  
  
  rds_dir<-paste("RDS\\",data_name,".rds",sep = "")
  
  bn <- readRDS(rds_dir)  
  
  
  
  # ------------ generate data
  
  
  set.seed(random_seed)     # use random seed
  
  dat<- rbn(bn, samples) 
  
  
  for (i in 1:length(dat)) 
  {
    dat[,i]<-as.numeric(dat[,i])
  }
  
  data_file<-paste("..\\Data\\",data_name,"_",samples,".txt",sep = "")
  
  
  write.table(dat-1,data_file,quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  
  
  
  # ------------ generate graph
  
  
  adjgraph<-matrix(0,length(bn),length(bn))
  
  for(i in 1:length(bn)) 
  {
    c <- bn[[i]]$children
    
    if( length(c) > 0 ) 
    {
      for (j in 1:length(c))
      {
        k<-which(names(bn)==c[j])
        adjgraph[i,k] = 1
      }
    }
  }
  
  
  graph_file<-paste("..\\Data\\",data_name,"_graph.txt",sep = "")
  
  
  write.table(adjgraph,graph_file,quote = FALSE,row.names = FALSE,col.names = FALSE)
  
}


