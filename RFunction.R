library('move')
library('geodist')
library('lubridate') #version x.7.y!
library('lutz')
library('sf')
library('sp')
library('rgeos')

rFunction = function(meth="buff", rad=NULL, dur=NULL, dur_unit="days", data, ...) {
  Sys.setenv(tz="UTC")
  names(data) <- make.names(names(data),allow_=FALSE)
  
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
    ix <- which(namesIndiv(data)=="remove")
    remove <- data[[ix]] #move object
    data <- data[[-ix]]
    remo <- TRUE
    logger.info(paste("Your data set contains", length(remove), "locations with the ID 'remove'. Clusters close (< rad) to those locations will be removed from your results."))
    logger.info(paste("Your remaining data set has", length(data), "locations of", length(namesIndiv(data)),"individuals."))
  }

  #tried to include recurse package here to pre-filter only revisited locations, but the runtime of getRecursions() was too long
  
  # cluster for all locations (not by ID)
  coos <- coordinates(data)
  
  if (meth=="buff")
  {
    data_eq <- spTransform(data,CRSobj=paste0("+proj=aeqd +lat_0=",mean(coos[,2])," +lon_0=",mean(coos[,1])," +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
    data_eq_buffer <- buffer(data_eq,rad)
    data_eq_buffer_disag <- disaggregate(data_eq_buffer) #this function needs the rgeos package to work with holes in the polygons
    
    data_eq_extr <- extract(data_eq_buffer_disag,SpatialPoints(data_eq))
    #plot(data_eq_buffer_disag,col=rainbow(26))
    #points(SpatialPoints(data_eq),col=data_eq_extr$poly.ID,pch=20,cex=5)
    memb <- data_eq_extr$poly.ID
  } else if (meth=="hclust")
  {
    dista <- geodist_vec(x1=coos[,1],y1=coos[,2],measure="vincenty") #unit=m, "geodesic" is probably better, but takes even longer
    
    #clu <- hclust(as.dist(dista),method="ward.D2") #measure in dendrogram is not distance
    clu <- hclust(as.dist(dista),method="average")
    #plot(as.dendrogram(clu), ylim = c(0,1000))
    #abline(h=400,col=2)
    memb <- cutree(clu,h=2*rad) #group membership for each location
  }

  memb <- as.character(memb) #changed to character for splitting possibility
  cluID_all <- unique(memb)
  data@data <- cbind(data@data,"clusterID"=memb)
  
  # split clusters with gaps larger than "maxgap" (multi-individual tracks in cluster) - NEED TO TEST!!!
  for (i in seq(along=cluID_all)) #for loop ok?
  {
    x <- cluID_all[i]
    datax <- data[data@data$clusterID==x]
    gapix <- which(difftime(timestamps(datax)[-1],timestamps(datax)[-length(datax)],unit=gap_unit)>maxgap)
    if (length(gapix)>0)
    {
      ends <- c(gapix,length(datax))
      for (i in rev(seq(along=gapix))) data@data$clusterID[data@data$clusterID==x][(ends[i]+1):ends[i+1]] <- paste0(x,".",i) #keeps orig. name for first component
    }
  }
  
  cluID_all <- unique(data@data$clusterID) #feed new clusters into list of all

  #remove clusters with duration below "dur"
  cluID <- apply(matrix(cluID_all), 1, function(x) ifelse(as.numeric(difftime(max(timestamps(data)[data@data$clusterID==x]),min(timestamps(data)[data@data$clusterID==x]),unit=dur_unit))>=dur, x, NA))
  
  cluID <- cluID[!is.na(cluID)]
  
  if (length(cluID)>0) #include here to calc. number of locations/bursts not in cluster TODO!
  {
    #midlon <- apply(matrix(cluID), 1, function(x) mean(coordinates(data[data@data$clusterID==x])[,1])) 
    #midlat <- apply(matrix(cluID), 1, function(x) mean(coordinates(data[data@data$clusterID==x])[,2])) 
    
    centrloc <- t(apply(matrix(cluID), 1, function(x) coordinates(data[data@data$clusterID==x])[min(which(rowMeans(geodist_vec(x1=coordinates(data[data@data$clusterID==x])[,1],y1=coordinates(data[data@data$clusterID==x])[,2],measure="vincenty"))==min(rowMeans(geodist_vec(x1=coordinates(data[data@data$clusterID==x])[,1],y1=coordinates(data[data@data$clusterID==x])[,2],measure="vincenty"))))),]))
    
    centrlon <- centrloc[,1]
    centrlat <- centrloc[,2]
    
    #take out clusters in rad radius around "remove"
    if (remo==TRUE) 
    {
      remo_dist <- geodist_vec(x1=coordinates(remove)[,1],y1=coordinates(remove)[,2],x2=centrlon,y2=centrlat,measure="vincenty")
      if (any(remo_dist<rad))
      {
        out <- which(remo_dist<rad,arr.ind=TRUE)[,2]
        cluID <- cluID[-out]
        centrlon <- centrlon[-out]
        centrlat <- centrlat[-out]
        logger.info(paste0(length(out)," clusters were removed from your results, because they were close (< ",rad," m) to the provided locations with ID 'remove'."))
      }
    }
    
    result <- data[data@data$clusterID %in% cluID] #these are all locations that are in a (non-remove) cluster with difftime>dur
    result@data <- result@data[,!sapply(result@data, function(x) all(is.na(x)))]
    result@data$clu.centr.long <- apply(matrix(result@data$clusterID), 1, function(x) centrlon[which(cluID==x)])
    result@data$clu.centr.lat <- apply(matrix(result@data$clusterID), 1, function(x) centrlat[which(cluID==x)])
    
    tz_info_result <- tz_lookup_coords(coordinates(result)[,2], coordinates(result)[,1], method = "accurate")
    result@data$timestamp.local <- apply(data.frame(timestamps(result),tz_info_result), 1, function(x) as.character(lubridate::with_tz(x[1], x[2])))
    result@data$local.timezone <- tz_info_result 
    result@data$date.local <- format(as.POSIXct(result@data$timestamp.local),format="%Y-%m-%d")
    result@data$time.local <- format(as.POSIXct(result@data$timestamp.local),format="%H:%M:%S")
    
    coo <- data.frame(coordinates(result))
    names(coo) <- c("location.long","location.lat")
    result.df <- data.frame(as.data.frame(result),coo)
    names(result.df) <- make.names(names(result.df),allow_=FALSE)
    clu.ix <- which(names(result.df) %in% c("clusterID","clu.centr.long","clu.centr.lat"))
    result.df <- data.frame(result.df[,clu.ix],result.df[,-clu.ix])
    heightix <- min(grep("height",names(result.df)))
    heightname <- names(result.df)[heightix]
    
    #cluster table
    n.locs <- apply(matrix(cluID), 1, function(x) length(which(result@data$clusterID==x)))
    n.ids <- apply(matrix(cluID), 1, function(x) length(unique(result.df$trackId[result.df$clusterID==x])))
    id.names <- apply(matrix(cluID), 1, function(x) paste(unique(result.df$trackId[result.df$clusterID==x]),collapse=", "))
    id.tags <- apply(matrix(cluID), 1, function(x) paste(unique(result.df$tag.local.identifier[result.df$clusterID==x]),collapse=", "))
    
    alldata.df <- as.data.frame(data)
    
    id.locs <- id.durs <- id.locsout <- id.locsBETout <- character(length(cluID))
    for (i in seq(along=cluID))
    {
      idsi <- as.character(unique(result.df$trackId[result.df$clusterID==cluID[i]]))
      id.locs[i] <- paste(apply(matrix(idsi), 1, function(x) length(result.df$trackId[result.df$trackId==x & result.df$clusterID==cluID[i]])),collapse=", ")
      id.durs[i] <- paste(apply(matrix(idsi), 1, function(x) round(as.numeric(difftime(max(result.df$timestamp[result.df$trackId==x & result.df$clusterID==cluID[i]],na.rm=TRUE),min(result.df$timestamp[result.df$trackId==x & result.df$clusterID==cluID[i]],na.rm=TRUE),units=dur_unit)),1)),collapse=", ")
      id.locsout[i] <- paste(apply(matrix(idsi), 1, function(x) length(which(alldata.df$timestamp[alldata.df$trackId==x] >= min(timestamps(data[data@data$clusterID==cluID[i]])) & alldata.df$timestamp[alldata.df$trackId==x] <= max(timestamps(data[data@data$clusterID==cluID[i]])))) - length(result.df$trackId[result.df$trackId==x & result.df$clusterID==cluID[i]])),collapse=", ") #locations outside of cluster in complete cluster interval
      id.locsBETout[i] <- paste(apply(matrix(idsi), 1, function(x) length(which(alldata.df$timestamp[alldata.df$trackId==x] >= min(alldata.df$timestamp[alldata.df$clusterID==cluID[i] & alldata.df$trackId==x]) & alldata.df$timestamp[alldata.df$trackId==x] <= max(alldata.df$timestamp[alldata.df$clusterID==cluID[i] & alldata.df$trackId==x]))) - length(result.df$trackId[result.df$trackId==x & result.df$clusterID==cluID[i]])),collapse=", ") #locations outside of cluster in indiv specific cluster interval
    }
    timestamp.start <- apply(matrix(cluID), 1, function(x) paste(as.character(min(timestamps(result[result@data$clusterID==x]))),"UTC"))
    timestamp.end <- apply(matrix(cluID), 1, function(x) paste(as.character(max(timestamps(result[result@data$clusterID==x]))),"UTC"))
    duration <- as.numeric(difftime(as.POSIXct(timestamp.end), as.POSIXct(timestamp.start),units=dur_unit))
    
    cluster.diameter.m <- apply(matrix(cluID), 1, function(x) max(geodist_vec(x1=coordinates(data[data@data$clusterID==x])[,1],y1=coordinates(data[data@data$clusterID==x])[,2],measure="vincenty"),na.rm=TRUE))
    realised.centr.radius.m <- apply(matrix(cluID), 1, function(x) max(geodist_vec(x1=coordinates(data[data@data$clusterID==x])[,1],y1=coordinates(data[data@data$clusterID==x])[,2],x2=centrlon[which(cluID==x)],y2=centrlat[which(cluID==x)],measure="vincenty"),na.rm=TRUE))
    
    tz_info_clu<- tz_lookup_coords(centrlat, centrlon, method = "accurate")
    timestamp.start.local <- apply(data.frame(timestamp.start,tz_info_clu), 1, function(x) as.character(lubridate::with_tz(x[1], x[2])))
    timestamp.end.local <- apply(data.frame(timestamp.end,tz_info_clu), 1, function(x) as.character(lubridate::with_tz(x[1], x[2])))
    
    clu_tab <- data.frame("cluster.ID"=cluID,n.locs,n.ids,id.tags,id.locs,id.durs,"centr.long"=centrlon,"centr.lat"=centrlat,timestamp.start.local,timestamp.end.local,"local.timezone"=tz_info_clu,duration,timestamp.start,timestamp.end,id.names,cluster.diameter.m,realised.centr.radius.m,id.locsout,id.locsBETout)
    names(clu_tab)[names(clu_tab)=="duration"] <- paste0("duration (",dur_unit,")")
    names(clu_tab)[names(clu_tab)=="id.durs"] <- paste0("id.durs (",dur_unit,")")
    
    o <- order(clu_tab$n.ids,clu_tab$n.locs,decreasing=TRUE)
    clu_tab <- clu_tab[o,]
    
    write.csv(clu_tab,file=paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Cluster_Table.csv"),row.names=FALSE)
    #write.csv(clu_tab,file="Cluster_Table.csv",row.names=FALSE)
    
    # finish points with clusters table, add n.ids and n.locs for Email Alert App
    result@data$n.ids <- apply(matrix(result@data$clusterID), 1, function(x) n.ids[which(cluID==x)])
    result@data$n.locs <- apply(matrix(result@data$clusterID), 1, function(x) n.locs[which(cluID==x)])
    result.df <- cbind(result.df,"n.ids"=result@data$n.ids,"n.locs"=result@data$n.locs)
    
    ixcsv <- which(c("clusterID","tag.local.identifier","n.ids","n.locs","timestamp.local","location.long","location.lat","date.local","time.local","local.timezone","trackId","ground.speed","heading",heightname,"clu.centr.long","clu.centr.lat") %in% names(result.df)) #fix if e.g. some data sets done have ground.speed or heading
    result.df.csv <- result.df[,c("clusterID","tag.local.identifier","n.ids","n.locs","timestamp.local","location.long","location.lat","date.local","time.local","local.timezone","trackId","ground.speed","heading",heightname,"clu.centr.long","clu.centr.lat")[ixcsv]]
    
    names(result.df.csv)[names(result.df.csv)=="trackId"] <- c("animalID")
    names(result.df.csv)[names(result.df.csv)=="tag.local.identifier"] <- c("tagID")
    result@data$animalID <- result.df.csv$animalID
    result@data$tagID <- result.df.csv$tagID
     
    #for utm locations we would need a separate App
    
    write.csv(result.df.csv,file=paste0(Sys.getenv(x = "APP_ARTIFACTS_DIR", "/tmp/"),"Points_With_Clusters.csv"),row.names=FALSE)
    #write.csv(result.df.csv,file="Points_With_Clusters.csv",row.names=FALSE) 
    
    selnames <- c("clusterID","tagID","n.ids","n.locs","timestamp.local","location.long","location.lat","date.local","time.local","local.timezone","animalID","ground.speed","heading",heightname,"clu.centr.long","clu.centr.lat")[ixcsv]
    result@data <- data.frame(result@data,coo)
    sel <- which(names(result@data) %in% selnames)
    result@data <- data.frame(result@data[,selnames],result@data[-sel])
   
    #force moveStack if only one ID (can lead to strange error)
    if (is(result,'Move')) {
      result <- moveStack(result,forceTz="UTC")
    }
    
  } else result <- NULL
  
  # the use of package recurse or adehabitatLT does only work on tracks, but not for clusters by all animals...
  return(result)
}
