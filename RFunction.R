library('move')
library('geodist')

rFunction = function(rad=NULL, dur=NULL, dur_unit="days", data, ...) {
  
  Sys.setenv(tz="UTC") 
  
  if (is.null(rad))
  {
    logger.info("Your cluster radius is not supplied. We use default 200 m.")
    rad <- 200
  }
  if (is.null(dur))
  {
    logger.info(paste0("Your minimum cluster duration is not supplied. We use 14 (",dur_unit,")  as default."))
    dur <- 14
  }
  
  #tried to include recurse pachage here to pre-filter only revisited locations, but the runtime of getRecursions() was too long
  
  # cluster for all locations (not by ID)
  coos <- coordinates(data)
  dista <- geodist_vec(x1=coos[,1],y1=coos[,2],measure="vincenty") #unit=m, "geodesic" is probably better, but takes even longer
  
  #clu <- hclust(as.dist(dista),method="ward.D2") #measure in dendrogram is not distance
  clu <- hclust(as.dist(dista),method="average")
  #plot(as.dendrogram(clu), ylim = c(0,1000))
  #abline(h=400,col=2)
  memb <- cutree(clu,h=2*rad) #group membership for each location
  
  data@data <- cbind(data@data,"clusterID"=memb)
  cluID_all <- unique(memb)
  
  cluID <- apply(matrix(cluID_all), 1, function(x) ifelse(as.numeric(difftime(max(timestamps(data)[data@data$clusterID==x]),min(timestamps(data)[data@data$clusterID==x]),unit=dur_unit))>=dur, x, NA))
  
  cluID <- cluID[!is.na(cluID)]
  
  if (length(cluID)>0)
  {
    result <- data[data@data$clusterID %in% cluID] #these are all locations that are in a cluster with difftime>dur
    midlon <- apply(matrix(cluID), 1, function(x) mean(coordinates(data[data@data$clusterID==x])[,1])) 
    midlat <- apply(matrix(cluID), 1, function(x) mean(coordinates(data[data@data$clusterID==x])[,2])) 
    
    result@data$cluster.mid.long <- apply(matrix(result@data$clusterID), 1, function(x) midlon[which(cluID==x)])
    result@data$cluster.mid.lat <- apply(matrix(result@data$clusterID), 1, function(x) midlat[which(cluID==x)])
    result.df <- data.frame(as.data.frame(result),coordinates(result))
    clu.ix <- which(names(result.df) %in% c("clusterID","cluster.mid.long","cluster.mid.lat"))
    result.df <- data.frame(result.df[,clu.ix],result.df[,-clu.ix])
    write.csv(result.df,file=paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Points_With_Clusters.csv"),row.names=FALSE)
    #write.csv(result.df,file="Points_With_Clusters.csv",row.names=FALSE) 
    
    #cluster table
    n.locs <- apply(matrix(cluID), 1, function(x) length(which(result@data$clusterID==x)))
    n.ids <- apply(matrix(cluID), 1, function(x) length(unique(result.df$trackId[result.df$clusterID==x])))
    id.names <- apply(matrix(cluID), 1, function(x) paste(unique(result.df$trackId[result.df$clusterID==x]),collapse=", "))
    timestamp.start <- apply(matrix(cluID), 1, function(x) paste(as.character(min(timestamps(result[result@data$clusterID==x]))),"UTC"))
    timestamp.end <- apply(matrix(cluID), 1, function(x) paste(as.character(max(timestamps(result[result@data$clusterID==x]))),"UTC"))
    duration <- as.numeric(difftime(as.POSIXct(timestamp.end), as.POSIXct(timestamp.start),units=dur_unit))
    # add local timestamps!!
    
    clu_tab <- data.frame("cluster.ID"=cluID,"mid.long"=midlon,"mid.lat"=midlat,timestamp.start,timestamp.end,duration,n.locs,n.ids,id.names)
    names(clu_tab)[names(clu_tab)=="duration"] <- paste0("duration (",dur_unit,")")
    write.csv(clu_tab,file=paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Cluster_Table.csv"),row.names=FALSE)
    #write.csv(clu_tab,file="Cluster_Table.csv",row.names=FALSE)
  } else result <- NULL
  
  # could add n.locs and duration for different individuals
  # add another option to use package recurse
  
  return(result)
}
