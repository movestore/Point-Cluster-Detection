library('move2')
library('geodist')
library('lubridate')
library('lutz')
library('sf')
library('terra')

# to do:
# 1. redo all with sf
# 2. add to output (and csvs): "most common location", i.e. for each cluster cut off at 4 or 5 decimals (round), make table and select location with most occurances (radius 5-10 or 20-30m should be ok) --> rather with extra radius and hclust (clusters within clusters)

rFunction = function(meth="buff", rad=NULL, dur=NULL, minloc=NULL, dur_unit="days", maxgap=1, gap_unit="days", clu_transm="all", new_dur=24, data, ...) {
  Sys.setenv(tz="UTC")
  time_now <- Sys.time()
  if (("timestamp" %in% (names(data)))==FALSE) data$timestamp <- mt_time(data)
  
  if (is.null(rad))
  {
    logger.info("Your cluster radius is not supplied. We use default 200 m.")
    rad <- 200
  }
  if (is.null(dur))
  {
    logger.info(paste0("Your minimum cluster duration is not supplied. We use 1 (",dur_unit,")  as default."))
    dur <- 1
  }
  if (is.null(minloc))
  {
    logger.info(paste0("Your minimum number of locations in a cluster is not supplied. We use 1 as default."))
    minloc <- 1
  }

  #tried to include recurse package here to pre-filter only revisited locations, but the runtime of getRecursions() was too long
  
  # cluster for all locations (not by ID)
  
  meth_buff <- function(data, rad)
  {
    coos <- st_coordinates(data)
    data_eq <- st_transform(data,crs=paste0("+proj=aeqd +lat_0=",mean(coos[,2])," +lon_0=",mean(coos[,1])," +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
    data_eq_buffer <- st_buffer(data_eq,dist=rad,allow_holes=TRUE)
    data_eq_buffer_uni <- st_union(data_eq_buffer)
    data_eq_buffer_unic <- st_cast(data_eq_buffer_uni,"POLYGON")
    
    #old sp code - keep for comparing with old results
    #data_eq <- spTransform(data,CRSobj=paste0("+proj=aeqd +lat_0=",mean(coos[,2])," +lon_0=",mean(coos[,1])," +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"))
    #data_eq_buffer <- buffer(data_eq,rad,doEdge=TRUE) #takens 20 min for 250000 locs
    #data_eq_buffer_disag <- disaggregate(data_eq_buffer) #this function needs the rgeos package to work with holes in the polygons (merges overlaying buffers) # this function comes from sp and had to be rewritten to sf::st_cast (check if does the same!)
    
    pts_eq <- st_coordinates(data_eq)
    len_pts <- dim(pts_eq)[1]
    data_eq_extr <- data.frame("point.ID"=numeric(len_pts),"poly.ID"=numeric(len_pts))
    stp <- 2000
    if (len_pts>stp)
    {
      vonbis <- matrix(c(seq(1,len_pts,by=stp),(seq(1,len_pts,by=stp))[-1]-1,len_pts),nc=2)
      #apply(vonbis, 1, function(x) data_eq_extr[x[1]:x[2],] <- as.data.frame(extract(data_eq_buffer_disag,pts_eq[x[1]:x[2]]))) #does not work
      for (i in seq(along=vonbis[,1])) data_eq_extr[vonbis[i,1]:vonbis[i,2],] <- extract(vect(data_eq_buffer_unic),pts_eq[vonbis[i,1]:vonbis[i,2],])
    } else data_eq_extr[1:len_pts,] <- extract(vect(data_eq_buffer_unic),pts_eq) #halleluja :)
    
    #plot(data_eq_buffer_disag,col=rainbow(26))
    #points(SpatialPoints(data_eq),col=data_eq_extr$poly.ID,pch=20,cex=5)
    memb <- data_eq_extr$poly.ID
    return(memb)
  }
  
  meth_hclust <- function(data,rad)
  {
    coos <- st_coordinates(data)
    if(dim(coos)[1]>20000) logger.info("Your data set is likely too large for this calculation. Try the buffer method, subset your data or use less locations/individuals for the analysis.")
    
    dista <- geodist_vec(x1=coos[,1],y1=coos[,2],measure="vincenty") #unit=m, "geodesic" is probably better, but takes even longer
    
    #clu <- hclust(as.dist(dista),method="ward.D2") #measure in dendrogram is not distance
    clu <- hclust(as.dist(dista),method="average")
    #plot(as.dendrogram(clu), ylim = c(0,1000))
    #abline(h=400,col=2)
    memb <- cutree(clu,h=2*rad) #group membership for each location
    return(memb)
  }
  
  if (meth=="buff")
  {
    memb <- meth_buff(data,rad)
  } else if (meth=="hclust") #cannot handle very much data, 
  {
    memb <- meth_hclust(data,rad)
  }

  memb <- as.character(memb) #changed to character for splitting possibility
  cluID_all <- unique(memb)
  data$clusterID <- memb
  
  # split clusters with gaps larger than "maxgap" (multi-individual tracks in cluster), but careful: the timestamps were not ordered if there are more than 1 idv - now adapted!
  for (i in seq(along=cluID_all))
  {
    #print(i)
    xx <- cluID_all[i]
    ixx <- which(data$clusterID==xx)
    datax <- data[ixx,]
    
    o <- order(mt_time(datax))
    ixo <- ixx[o]
    gapixo <- which(difftime((mt_time(datax)[o])[-1],(mt_time(datax)[o])[-dim(datax)[1]],unit=gap_unit)>maxgap)
    #gapix <- which(difftime(timestamps(datax)[-1],timestamps(datax)[-length(datax)],unit=gap_unit)>maxgap) #did not combine timestamps of multiple IDs
    if (length(gapixo)>0)
    {
      endso <- c(gapixo,dim(datax)[1])
      for (j in seq(along=gapixo)) 
      {
        ixo_endsoj <- ixo[(endso[j]+1):endso[j+1]]
        data$clusterID[ixo_endsoj] <- paste0(xx,".",j) #keeps orig. name for first component
      }
    }
  }
  
  cluID_all <- unique(data$clusterID) #feed new clusters into list of all
  
  #remove clusters with duration below "dur" or number of locs within below "minloc" (both conditions must be met for a cluster to be kept!)
  cluID <- apply(matrix(cluID_all), 1, function(x) ifelse(as.numeric(difftime(max(mt_time(data)[data$clusterID==x]),min(mt_time(data)[data$clusterID==x]),unit=dur_unit))>=dur & length(mt_time(data)[data$clusterID==x])>=minloc, x, NA))
  
  cluID <- cluID[!is.na(cluID)]
  
  ### if asked for reduced transmission, here cluster with data minus new locations and compare (note that cannot user clusterIDs, but need locations)
  ## strategy: look through the clusters of "all", for each check if all old locations are in a same cluster of the "old" cluster set, if so, then mark to remove from cluID list

  if (clu_transm=="new")
  {
    logger.info(paste("You have selected to only calculate and transmit clusters that cannot be found with the 'old' data, excluding the past",new_dur,"hours."))
    
    dataold <- data[mt_time(data) < time_now-hours(new_dur),]
    if (meth=="buff")
    {
      membold <- meth_buff(dataold,rad)
    } else if (meth=="hclust") #cannot handle very much data, 
    {
      membold <- meth_hclust(dataold,rad)
    }
    
    membold <- as.character(membold) #changed to character for splitting possibility
    cluID_allold <- unique(membold)
    dataold$clusterID <- membold
    
    #must match both data sets with trackid and timestamp
    timetrack_all <- paste(mt_time(data),mt_track_id(data),sep=",")
    timetrack_old <- paste(mt_time(dataold),mt_track_id(dataold),sep=",")
    
    old_cluix <- numeric() # indices of "old" clusters that can be removed from "all" set
    for (k in seq(along=cluID)) #now look at each of the "all" clusters
    {
      locix_k <- which(memb==cluID[k])
      sel_old <- difftime(time_now,mt_time(data)[locix_k],units="hours") > new_dur #which are old locations
      if (sum(sel_old)>0) 
      {
        locix_ko <- locix_k[sel_old]
        locix_ko_old <- which(timetrack_old %in% timetrack_all[locix_ko])
        if (length(unique(membold[locix_ko_old]))==1) old_cluix <- c(old_cluix,k)
      }
    }
    cluID <- cluID[-old_cluix]
    logger.info(paste(length(old_cluix),"clusters were removed from the complete list, as they were fully included in clusters calculated from the data excluding the previous",new_dur,"hours.",length(cluID),"'new' clusters are in the remaining list."))
  }
# from here all is based on this list of cluIDs
 
  if (length(cluID)>0) #include here to calc. number of locations/bursts not in cluster - see below
  {
    #midlon <- apply(matrix(cluID), 1, function(x) mean(coordinates(data[data@data$clusterID==x])[,1])) 
    #midlat <- apply(matrix(cluID), 1, function(x) mean(coordinates(data[data@data$clusterID==x])[,2])) 
    
    centrloc <- t(apply(matrix(cluID), 1, function(x) st_coordinates(data[data$clusterID==x,])[min(which(rowMeans(geodist_vec(x1=st_coordinates(data[data$clusterID==x,])[,1],y1=st_coordinates(data[data$clusterID==x,])[,2],measure="vincenty"))==min(rowMeans(geodist_vec(x1=st_coordinates(data[data$clusterID==x,])[,1],y1=st_coordinates(data[data$clusterID==x,])[,2],measure="vincenty"))))),]))

    centrlon <- centrloc[,1]
    centrlat <- centrloc[,2]
    
    #include most used location (cut locs to 4 decimals and select most frequent loc. If more than one with exactly max usage, select the one closest to centrloc)
    maxuselon <- maxuselat <- numeric(length(cluID)) #same order as cluID and centrloc
    for (j in seq(along=cluID))
    {
      locs4 <- data.frame(round(st_coordinates(data[data$clusterID==cluID[j],]),digits=4))
      locs4a <- aggregate(cbind(locs4[0],numdup=1), locs4, length)
      maxuse0 <- which(locs4a$numdup==max(locs4a$numdup))
      if (length(maxuse0>1)) 
      {
        diss <- apply(matrix(maxuse0), 1, function(x) geodist_vec(x1=as.numeric(locs4a[x,1:2]),y1=centrloc[j,],measure="vincenty")[2])
        maxuse <- maxuse0[which(diss==min(diss))[1]]
      } else maxuse <- maxuse0
      maxuselon[j] <- locs4a[maxuse,1]
      maxuselat[j] <- locs4a[maxuse,2]
    }
    
    # load remove locations file and take out clusters in rad radius around them
    #remo_file_path <- "./data/local_app_files/uploaded-app-files/" 
    #remo_file_path <- getAppFilePath("remo_sites") #returns NULL if no fallback
    remo_file_name <- getAuxiliaryFilePath("remo_sites") #returns NULL if no fallback

    #if ("remove.csv" %in% list.files(remo_file_path)) #no sure what this would have done without fallback
    if (!is.null(remo_file_name))
    {
     remo <- read.csv(remo_file_name,header=TRUE)
     if ("latitude" %in% names(remo) & "longitude" %in% names(remo))
     {
       logger.info(paste("You have uploaded a file with",length(remo$latitude),"locations to remove from the set of clusters."))
       remo_dist <-  geodist_vec(x1=remo$longitude,y1=remo$latitude,x2=centrlon,y2=centrlat,measure="vincenty")
       if (any(remo_dist<rad))
       {
         out <- which(remo_dist<rad,arr.ind=TRUE)[,2]
         cluID <- cluID[-out]
         centrlon <- centrlon[-out]
         centrlat <- centrlat[-out]
         maxuselon <- maxuselon[-out]
         maxuselat <- maxuselat[-out]
         logger.info(paste0(length(out)," clusters were removed from your results, because they were close (< ",rad," m) to the provided remove locations'."))
       } else logger.info("None of your provided remove locations were close to a cluster.")
     } else logger.info ("Your csv file does not include the required columns longitude and latitude. Your required remove locations cannot be taken into account.")
    }
    
    result <- data[data$clusterID %in% cluID,] #these are all locations that are in any (non-remove) cluster with difftime>dur
    result <- result[,!sapply(result, function(x) all(is.na(x))),] #remove columns with no info
    result$clu.centr.long <- apply(matrix(result$clusterID), 1, function(x) centrlon[which(cluID==x)])
    result$clu.centr.lat <- apply(matrix(result$clusterID), 1, function(x) centrlat[which(cluID==x)])
    result$clu.maxuse.long <- apply(matrix(result$clusterID), 1, function(x) maxuselon[which(cluID==x)])
    result$clu.maxuse.lat <- apply(matrix(result$clusterID), 1, function(x) maxuselat[which(cluID==x)])
    
    tz_info_result <- tz_lookup_coords(st_coordinates(result)[,2], st_coordinates(result)[,1], method = "accurate")
    result$timestamp.local <- apply(data.frame(mt_time(result),tz_info_result), 1, function(x) as.character(lubridate::with_tz(x[1], x[2])))
    result$local.timezone <- tz_info_result 
    result$date.local <- format(as.POSIXct(result$timestamp.local),format="%Y-%m-%d")
    result$time.local <- format(as.POSIXct(result$timestamp.local),format="%H:%M:%S")
    
    result$location.long <- st_coordinates(result)[,1]
    result$location.lat <- st_coordinates(result)[,2]
    result <- result[,!sapply(result, function(x) all(is.na(x)))] #take NA out
    #names(result.df) <- make.names(names(result.df),allow_=FALSE)
    clu.ix <- which(names(result) %in% c("clusterID","clu.centr.long","clu.centr.lat","clu.maxuse.long","clu.maxuse.lat"))
    result <- result[,c(clu.ix,seq(along=result)[-clu.ix])]
    heightix <- grep("height",names(result))
    if(length(heightix)>0) heightname <- names(result)[min(heightix)] else heightname <- NA
    
    #add tag_local_identifier to results
    if("tag_local_identifier" %in% names(mt_track_data(result))) result <- mt_as_event_attribute(result,"tag_local_identifier",.keep=TRUE)
    #if("tag.local.identifier" %in% names(mt_track_data(result))) result <- mt_as_event_attribute(result,"tag.local.identifier",.keep=TRUE)

    # add track resolutions from file
    fixrates_file_name <- getAuxiliaryFilePath("track_fixrates") #returns NULL if no fallback
    
    if (!is.null(fixrates_file_name))
    {
      fixr <- read.csv(fixrates_file_name,header=TRUE)
      if ("trackid" %in% names(fixr) & "fixrate" %in% names(fixr))
      {
        logger.info(paste("You have uploaded a file with",length(fixr$fixrate),"trackids to which fix rates (unit=per hour; must be regular) will be added in the output rds."))
        for (j in seq(along=fixr$trackid)) 
        {
          if (fixr$trackid[j] %in% mt_track_id(result)) 
          {
            mt_track_data(result)[mt_track_id(result)==fixr$trackid[j],] <- fixr$fixrate[j]
            logger.info(paste("fixrate of track",fixr$trackid[j],"added to rds output."))
          } else logger.info(paste("trackid",fixr$trackid[j],"is not available in the data set."))
        result <- mt_as_event_attribute(result,"fixrate",.keep=TRUE)
        } else logger.info("None of your provided trackids matched to the input data set, please check if their spelling is correct.")
      } else logger.info ("Your csv file does not include the required columns trackid and fixrate. These cannot be added to your results object.")
    }

    #cluster table
    n.locs <- apply(matrix(cluID), 1, function(x) length(which(result$clusterID==x)))
    n.ids <- apply(matrix(cluID), 1, function(x) length(unique(mt_track_id(result)[result$clusterID==x])))
    id.names <- apply(matrix(cluID), 1, function(x) paste(unique(mt_track_id(result)[result$clusterID==x]),collapse=", "))
    id.tags <- rep(NA,length(n.ids))
    if("tag_local_identifier" %in% names(result)) id.tags <- apply(matrix(cluID), 1, function(x) paste(unique(result[result$clusterID==x,]$tag_local_identifier),collapse=", "))
    #if("tag.local.identifier" %in% names(result)) id.tags <- apply(matrix(cluID), 1, function(x) paste(unique(result[result$clusterID==x,]$tag.local.identifier),collapse=", "))

    # add cumlag.locs (sum of timelags of each location in cluster) - addition from Sep 2022
    #result$timelag <- unlist(lapply(timeLag(result, units=dur_unit), function(x) c(as.vector(x), NA)))
    result$timelag <- mt_time_lags(result,units=dur_unit)
    cumlag.locs <-  apply(matrix(cluID), 1, function(x) sum((result$timelag[which(result$clusterID==x)])[-sum(result$clusterID==x)],na.rm=TRUE)) #take out timelag of last location as that is movement to next cluster (previously done by unlist...)
    
    alldata.df <- as.data.frame(data)
    
    id.locs <- id.durs <- id.locsout <- id.locsBETout <- id.nrevs <- id.durBETrevs <- id.maxdistBETrevs <- character(length(cluID))
    n.revs <- n.locsout <- numeric(length(cluID))
    for (i in seq(along=cluID))
    {
      idsi <- as.character(unique(mt_track_id(result)[result$clusterID==cluID[i]])) #individuals that use this cluster
      #if(length(idsi)>1) print(i)
      id.locs[i] <- paste(apply(matrix(idsi), 1, function(x) length(mt_track_id(result)[mt_track_id(result)==x & result$clusterID==cluID[i]])),collapse=", ")
      id.durs[i] <- paste(apply(matrix(idsi), 1, function(x) round(as.numeric(difftime(max(mt_time(result)[mt_track_id(result)==x & result$clusterID==cluID[i]],na.rm=TRUE),min(mt_time(result)[mt_track_id(result)==x & result$clusterID==cluID[i]],na.rm=TRUE),units=dur_unit)),1)),collapse=", ")
      id.locsout[i] <- paste(apply(matrix(idsi), 1, function(x) length(which(mt_time(data)[mt_track_id(data)==x] >= min(mt_time(data)[data$clusterID==cluID[i]]) & mt_time(data)[mt_track_id(data)==x] <= max(mt_time(data)[data$clusterID==cluID[i]]))) - length(mt_track_id(result)[mt_track_id(result)==x & result$clusterID==cluID[i]])),collapse=", ") #locations outside of cluster in complete cluster interval
      id.locsBETout[i] <- paste(apply(matrix(idsi), 1, function(x) length(which(mt_time(data)[mt_track_id(data)==x] >= min(mt_time(data)[data$clusterID==cluID[i] & mt_track_id(data)==x]) & mt_time(data)[mt_track_id(data)==x] <= max(mt_time(data)[data$clusterID==cluID[i] & mt_track_id(data)==x]))) - length(mt_track_id(result)[mt_track_id(result)==x & result$clusterID==cluID[i]])),collapse=", ") #locations outside of cluster in indiv specific cluster interval
      
      n.locsout[i] <- sum(as.numeric(trimws(strsplit(as.character(id.locsout[i]),",")[[1]])),na.rm=TRUE) # sum of all locations outside of cluster (all idvs. that use cluster, complete cluster duration)
      
      id.nrevisits <- id.durBETrevisits <- id.maxdistBETrevisits <- character(0)
      for (j in seq(along=idsi))
      {
        clutrackj <- data[mt_track_id(data)==idsi[j] & mt_time(data) >= min(mt_time(data)[data$clusterID==cluID[i] & mt_track_id(data)==idsi[j]]) & mt_time(data) <= max(mt_time(data)[data$clusterID==cluID[i] & mt_track_id(data)==idsi[j]]),]
        #print(clutrackj[,c("trackId","timestamp","clusterID")])
        clutr_p1 <- clutrackj[-1,]
        clutr_m1 <- clutrackj[-dim(clutrackj)[1],]

        lastix <- which(clutrackj$clusterID[-dim(clutrackj)[1]]==cluID[i] & clutr_p1$clusterID!=cluID[i]) #when leave cluster
        firstix <- which(clutrackj$clusterID[-1]==cluID[i] & clutr_m1$clusterID!=cluID[i]) +1 #when enter cluster anew
          
        #need to append with comma
        if (length(firstix)>0)
        {
          #if (length(firstix)>1) print(paste(i,j))
          id.nrevisits <- c(id.nrevisits,length(firstix))
          id.durBETrevisits <- c(id.durBETrevisits,paste(difftime(mt_time(clutrackj)[firstix],mt_time(clutrackj)[lastix],units="hours"),collapse=","))
          id.maxdistBETrevisits <- c(id.maxdistBETrevisits,paste(apply(matrix(seq(along=firstix)), 1, function(x) round(max(geodist_vec(x1=centrlon[i],y1=centrlat[i],x2=st_coordinates(clutrackj)[(lastix[x]+1):(firstix[x]-1),1],y2=st_coordinates(clutrackj)[(lastix[x]+1):(firstix[x]-1),2])),digits=2)),collapse=","))
        } else
        {
          id.nrevisits <- c(id.nrevisits,NA)
          id.durBETrevisits <- c(id.durBETrevisits,NA)
          id.maxdistBETrevisits <- c(id.maxdistBETrevisits,NA)                
        }
        n.revs[i] <- sum(as.numeric(id.nrevisits),na.rm=TRUE) # sum of revisits (all idvs.)
        id.nrevs[i] <- paste(id.nrevisits,collapse=", ")  # number of individual revisits (zero if only once at cluster)
        id.durBETrevs[i] <- paste(id.durBETrevisits,collapse=", ") # duration outside of cluster before revisit (timediff between last loc in cluster until first loc back in cluster) - NA if no revisits
        id.maxdistBETrevs[i] <- paste(id.maxdistBETrevisits,collapse=", ") #max dist to cluster centre during gap before revisit
      }
    }
    timestamp.start <- apply(matrix(cluID), 1, function(x) paste(as.character(min(mt_time(result)[result$clusterID==x])),"UTC"))
    timestamp.end <- apply(matrix(cluID), 1, function(x) paste(as.character(max(mt_time(result)[result$clusterID==x])),"UTC"))
    duration <- as.numeric(difftime(as.POSIXct(timestamp.end), as.POSIXct(timestamp.start),units=dur_unit))

    cluster.diameter.m <- apply(matrix(cluID), 1, function(x) max(geodist_vec(x1=st_coordinates(data[data$clusterID==x,])[,1],y1=st_coordinates(data[data$clusterID==x,])[,2],measure="vincenty"),na.rm=TRUE))
    realised.centr.radius.m <- apply(matrix(cluID), 1, function(x) max(geodist_vec(x1=st_coordinates(data[data$clusterID==x,])[,1],y1=st_coordinates(data[data$clusterID==x,])[,2],x2=centrlon[which(cluID==x)],y2=centrlat[which(cluID==x)],measure="vincenty"),na.rm=TRUE))
    
    tz_info_clu<- tz_lookup_coords(centrlat, centrlon, method = "accurate")
    timestamp.start.local <- apply(data.frame(timestamp.start,tz_info_clu), 1, function(x) as.character(lubridate::with_tz(x[1], x[2])))
    timestamp.end.local <- apply(data.frame(timestamp.end,tz_info_clu), 1, function(x) as.character(lubridate::with_tz(x[1], x[2])))
    
    clu_tab <- data.frame("cluster.ID"=cluID,n.locs,cumlag.locs,n.ids,id.tags,id.locs,id.durs,"maxuse.long"=maxuselon,"maxuse.lat"=maxuselat,"centr.long"=centrlon,"centr.lat"=centrlat,timestamp.start.local,timestamp.end.local,"local.timezone"=tz_info_clu,duration,timestamp.start,timestamp.end,id.names,cluster.diameter.m,realised.centr.radius.m,n.locsout,id.locsout,id.locsBETout,n.revs,id.nrevs,id.durBETrevs,id.maxdistBETrevs)
    
    names(clu_tab)[names(clu_tab)=="cumlag.locs"] <- paste0("cumlag.locs (",dur_unit,")")
    names(clu_tab)[names(clu_tab)=="duration"] <- paste0("duration (",dur_unit,")")
    names(clu_tab)[names(clu_tab)=="id.durs"] <- paste0("id.durs (",dur_unit,")")
    
    o <- order(clu_tab$n.ids,clu_tab$n.locs,decreasing=TRUE)
    clu_tab <- clu_tab[o,]
    
    write.csv(clu_tab,file=appArtifactPath("Cluster_Table.csv"),row.names=FALSE)

    # finish points with clusters table, add n.ids and n.locs for Email Alert App (and n.revs and n.locsout for Earth Ranger)
    result$n.ids <- apply(matrix(result$clusterID), 1, function(x) n.ids[which(cluID==x)])
    result$n.locs <- apply(matrix(result$clusterID), 1, function(x) n.locs[which(cluID==x)])
    result$cumlag.locs <- apply(matrix(result$clusterID), 1, function(x) cumlag.locs[which(cluID==x)])
    result$n.revs <- apply(matrix(result$clusterID), 1, function(x) n.revs[which(cluID==x)])
    result$n.locsout <- apply(matrix(result$clusterID), 1, function(x) n.locsout[which(cluID==x)])

    #inlcude only the listed variables into result csv, in given order
    expnames <- c("test","clusterID","tag_local_identifier","n.ids","n.locs","cumlag.locs","n.locsout","n.revs","timestamp.local","location.long","location.lat","date.local","time.local","local.timezone","individual_name_deployment_id","ground_speed","heading","clu.centr.long","clu.centr.lat","clu.maxuse.long","clu.maxuse.lat") 
    result.df <- data.frame(result)
    if (!is.na(heightname))
    {
      ixcsv <- as.numeric(apply(matrix(c(expnames,heightname)), 1, function(x) which(names(result) == x)))
      if (any(is.na(ixcsv))) ixcsv <- ixcsv[-which(is.na(ixcsv))]
      result.csv <- result.df[,ixcsv]
    } else
    {
      ixcsv <- as.numeric(apply(matrix(expnames), 1, function(x) which(names(result) == x)))
      if (any(is.na(ixcsv))) ixcsv <- ixcsv[-which(is.na(ixcsv))]
      result.csv <- result.df[,ixcsv]
    }
    
    names(result.csv)[names(result.csv)=="individual_name_deployment_id"] <- c("animalID")
    names(result.csv)[names(result.csv)=="tag_local_identifier"] <- c("tagID")
    result$animalID <- result.csv$animalID
    result$tagID <- result.csv$tagID
     
    #for utm locations we would need a separate App
    
    write.csv(result.csv,file=appArtifactPath("Points_With_Clusters.csv"),row.names=FALSE)

    if (!is.na(heightname)) selnames <- c("clusterID","tagID","n.ids","n.locs","cumlag.locs","n.locsout","n.revs","timestamp.local","location.long","location.lat","date.local","time.local","local.timezone","animalID","ground.speed","heading",heightname,"clu.centr.long","clu.centr.lat","clu.maxuse.long","clu.maxuse.lat") else selnames <- c("clusterID","tagID","n.ids","n.locs","cumlag.locs","n.locsout","n.revs","timestamp.local","location.long","location.lat","date.local","time.local","local.timezone","animalID","ground.speed","heading","clu.centr.long","clu.centr.lat","clu.maxuse.long","clu.maxuse.lat")
    
    sel <- as.numeric(apply(matrix(selnames), 1, function(x) which(names(result) == x)))
    if (any(is.na(sel))) sel <- sel[-which(is.na(sel))]
    result <- result[,c(sel,seq(along=result)[-sel])]
    
  } else result <- data[0,]
  
  # the use of package recurse or adehabitatLT does only work on tracks, but not for clusters by all animals...
  return(result)
}