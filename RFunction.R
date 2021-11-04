library('move')
library('geodist')
library('lubridate')
library('lutz')
library('sf')

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
  
  # take out "remove" locations from data if there are any
  remo <- FALSE
  if (any(namesIndiv(data)=="remove"))
  {
    data.split <- move::split(data)
    ix <- which(names(data.split)=="remove")
    remove <- data.split[[ix]] #move object
    data <- moveStack(data.split[-ix])
    remo <- TRUE
    logger.info(paste0("Your data set contains", length(remove), "locations with the ID 'remove'. Clusters close (< rad) to those locations will be removed from your results."))
  }

  #tried to include recurse package here to pre-filter only revisited locations, but the runtime of getRecursions() was too long
  
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
    midlon <- apply(matrix(cluID), 1, function(x) mean(coordinates(data[data@data$clusterID==x])[,1])) 
    midlat <- apply(matrix(cluID), 1, function(x) mean(coordinates(data[data@data$clusterID==x])[,2])) 
    
    #take out clusters in rad radius around "remove"
    if (remo==TRUE) 
    {
      remo_dist <- geodist_vec(x1=coordinates(remove)[,1],y1=coordinates(remove)[,2],x2=midlon,y2=midlat,measure="vincenty")
      if (any(remo_dist<rad))
      {
        out <- which(remo_dist<rad,arr.ind=TRUE)[,2]
        cluID <- cluID[-out]
        midlon <- midlon[-out]
        midlat <- midlat[-out]
        logger.info(paste0(out," clusters were removed from your results, because they were close (< rad) to the provided locations with ID 'remove'."))
      }
    }
    
    result <- data[data@data$clusterID %in% cluID] #these are all locations that are in a (non-remove) cluster with difftime>dur
    result@data <- result@data[,!sapply(result@data, function(x) all(is.na(x)))]
    result@data$cluster.mid.long <- apply(matrix(result@data$clusterID), 1, function(x) midlon[which(cluID==x)])
    result@data$cluster.mid.lat <- apply(matrix(result@data$clusterID), 1, function(x) midlat[which(cluID==x)])
    
    tz_info_result <- tz_lookup_coords(coordinates(result)[,2], coordinates(result)[,1], method = "accurate")
    result@data$timestamp.local <- apply(data.frame(timestamps(result),tz_info_result), 1, function(x) as.character(lubridate::with_tz(x[1], x[2])))
    
    result.df <- data.frame(as.data.frame(result),coordinates(result))
    clu.ix <- which(names(result.df) %in% c("clusterID","cluster.mid.long","cluster.mid.lat"))
    result.df <- data.frame(result.df[,clu.ix],result.df[,-clu.ix])
    write.csv(result.df,file=paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Points_With_Clusters.csv"),row.names=FALSE)
    #write.csv(result.df,file="Points_With_Clusters.csv",row.names=FALSE) 
    
    #cluster table
    n.locs <- apply(matrix(cluID), 1, function(x) length(which(result@data$clusterID==x)))
    n.ids <- apply(matrix(cluID), 1, function(x) length(unique(result.df$trackId[result.df$clusterID==x])))
    id.names <- apply(matrix(cluID), 1, function(x) paste(unique(result.df$trackId[result.df$clusterID==x]),collapse=", "))
    
    id.locs <- id.durs <- character(length(cluID))
    for (i in seq(along=cluID))
    {
      idsi <- as.character(unique(result.df$trackId[result.df$clusterID==cluID[i]]))
      id.locs[i] <- paste(apply(matrix(idsi), 1, function(x) length(result.df$trackId[result.df$trackId==x & result.df$clusterID==cluID[i]])),collapse=", ")
      id.durs[i] <- paste(apply(matrix(idsi), 1, function(x) round(as.numeric(difftime(max(result.df$timestamp[result.df$trackId==x & result.df$clusterID==cluID[i]]),min(result.df$timestamp[result.df$trackId==x & result.df$clusterID==cluID[i]]),units="days")),2)),collapse=", ")
    }
    timestamp.start <- apply(matrix(cluID), 1, function(x) paste(as.character(min(timestamps(result[result@data$clusterID==x]))),"UTC"))
    timestamp.end <- apply(matrix(cluID), 1, function(x) paste(as.character(max(timestamps(result[result@data$clusterID==x]))),"UTC"))
    duration <- as.numeric(difftime(as.POSIXct(timestamp.end), as.POSIXct(timestamp.start),units=dur_unit))
    
    tz_info_clu<- tz_lookup_coords(midlat, midlon, method = "accurate")
    timestamp.start.local <- apply(data.frame(timestamp.start,tz_info_clu), 1, function(x) as.character(lubridate::with_tz(x[1], x[2])))
    timestamp.end.local <- apply(data.frame(timestamp.end,tz_info_clu), 1, function(x) as.character(lubridate::with_tz(x[1], x[2])))
    
    clu_tab <- data.frame("cluster.ID"=cluID,"mid.long"=midlon,"mid.lat"=midlat,timestamp.start,timestamp.end,timestamp.start.local,timestamp.end.local,duration,n.locs,n.ids,id.names,id.locs,id.durs)
    names(clu_tab)[names(clu_tab)=="duration"] <- paste0("duration (",dur_unit,")")
    names(clu_tab)[names(clu_tab)=="id.durs"] <- paste0("id.durs (",dur_unit,")")
    write.csv(clu_tab,file=paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Cluster_Table.csv"),row.names=FALSE)
    #write.csv(clu_tab,file="Cluster_Table.csv",row.names=FALSE)
  } else result <- NULL
  
  # the use of package recurse or adehabitatLT does only work on tracks, but not for clusters by all animals...
  
  return(result)
}
