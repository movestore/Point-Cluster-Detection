# Multiple Animal Cluster Detection

MoveApps

Github repository: *github.com/movestore/Point-Cluster_Detection*

## Description
Detection of point clusters, where possibly more than one animal returns to within a specified time interval. Provides a table of each cluster with the times, duration, number of locations and animals. Clusters close to locations provided in a file 'remove.csv' can be excluded. By setting, only 'new' clusters can be transmitted.

## Documentation
This App uses two types of clustering (buffer overlap or hierarchical clustering) for the detection of point clusters where one or more animals return to repeatedly within a specified time frame. 

For buffer clustering all locations are transformed into a circular spatial polygon of radius `cluster radius`. Then all overlapping polygons (i.e. of which the center locations are less than `2 * cluster radius` apart) are joined into clusters. Only clusters that were used for at least the specified number of hours/days/weeks (dur parameter) are returned. Clusters with time gaps between successive (indiv. independent) locations above the user defined maximum (maxgap parameter) are split. Resulting partial clusters are kept if they fulfill the minimum duration requirement (dur parameter).

For hierarchical clustering the `average` method is used, i.e. the clusters are defined at the minimum average distance between all locations of the clusters. Clusters are selected to have at least a radius of `cluster radius` or be `2 * cluster radius` apart. Only clusters that were used for at least the specified number of hours/days/weeks are returned.

If one has uploaded a local file called 'remove.csv' of locations (column names must be `longitute` and `latitude`!), this App will automatically exclude those clusters that have centre points less than `rad` metre from the locations of 'remove.csv'. That way, fixed stations like nests or other points of attraction can be excluded from the results of this cluster analysis.

According to the maxgap parameter, clusters with usage gaps larger than the defined value are split.

A cluster overview table is returned as a .csv artefact to download. It includes for each cluster the most central location (minimum distance to all other locations), the most used location (if rounded to 4 decimals (+/-11m) which location is most used; note: if there are several candidates, then the one closest to the central location), timestamps of first and last location (UTC and local time), duration, cluster diameter, realised cluster radius (related to most central location), number of locations, cumulative timelag of locations, number of animals, the names of those animals, their tag numbers and their respective duration and number of locations in the cluster. See definition of all parameters in the list(s) below.

The output of the App includes the locations that could be attributed to a cluster (that fulfilled the minimum duration requirement), including the clusterID, number of animals using the cluster, number of locations in and outside the cluster (in the time that the cluster was used and by the individuals that used the clsuter) and number of revisits. This dataset is also returned as a .csv artefact to download, including local timestamps.

If selected, only 'new' clusters can be returned. 'New' means, that those clusters were not found if excluding all locations from the last Y hours (in relation to time or App run). Those clusters from the complete cluster set, which were completely embedded in one cluster of the 'old' data set are not returned. This setting is most useful for regular, scheduled runs of this App and integration into EarthRanger, so as to minimise data traffic.

### Input data
move2 location object

### Output data
move2 location object


### Products (Artefacts)
`Cluster_Table.csv`: Overview of properties of detected point clusters.

`Points_With_Clusters.csv`: Result data set as .csv, with all locations in clusters.

### Table column definitions

#### `Cluster_Table.csv`

clusterID:				ID of cluster

centr.long:				longitude of most central positon in cluster (minimum average distance to all other locations)

centr.lat:				latitutde of most central position in cluster (minimum average distance to all other locations)
	
maxuse.long:			longitude of maximum used positon in cluster (most locations if rounded to 4 decimals)

maxuse.lat:				latitutde of maximum used position in cluster (most locations if rounded to 4 decimals)

timestamp.start:		time of first location in cluster (UTC, format=YYYY-MM-DD HH:MM:SS)

timestamp.end:			time of last location in cluster (UTC, format=YYYY-MM-DD HH:MM:SS)

timestamp.start.local:	local timestamp of first location in cluster (format=YYYY-MM-DD HH:MM:SS)

timestamp.end.local:	local timestamp of last location in cluster (format=YYYY-MM-DD HH:MM:SS)

local.timezone:			name of local timezone (incl. closest large city)

duration (hours):		time duration between first and last location in cluster (unit=hours)

n.locs:					number of locations in cluster (independent of individual)

cumlag.locs:			cumulative timelag (duration to next logged location) of all locations in cluster to which this location belongs (unit as `dur_unit`; useful if input data are tracks with different resolutions; warning: data gaps can cause problemes here, always compare with n.locs

n.ids:					number of different individuals using the cluster

id.names:				animalIDs of all indiviuals using the cluster (comma separated)

id.tags:				tagIDs of the tags of all individuals using the cluster (comma separated)

id.locs:				number of locations of each individual in the cluster (comma separated, same order as id.names)

id.durs (hours):		time duration between first and last location in cluster of each individual (unit=hours, same order as id.names, zero if only one location of the individual in the cluster)

cluster.diameter.m:		maximum pairwise distance between all locations of the cluster (unit=metres)

realised.centr.radius.m:	realised cluster radius towards most central position, i.e. maximum distance of cluster locations to most central position in cluster (see above, unit=metres)

n.locsout:				sum of locations outside of the cluster of all individuals that used the cluster in the time interval that the cluster was used. This can be used to calculate fertility of the cluster (fertility = n.locs - n.locsout)

id.locsout: 			number of locations of each individual (see names/tags in id.names/tags) that are outside of the cluster in the time interval that the cluster is used by any individual (i.e. can be edgle locations before/after the respective individual has started/ended using the cluster

id.locsBETout: 			number of locations of each individual (see names/tags in id.names/tags) that are outside of the cluster in the time interval that the cluster is used by that specific individual (i.e. between the first and last visit of the indiviual to that cluster)

n.revs: 				number of revisits of cluster (i.e. number of times individuals leave cluster and come back, summed over all individuals that use the cluster

id.nrevs: 				individual number of revisits (i.e. number of time an individual leaves the cluster and comes back). Comma separated by individual, see id.names/id.tags for order

id.durBETrevs: 			duration of time outside of cluster for each leave/revisit event. Comma separated by individual (with space) and revisit event (without space). See id.names/id.tags for order

id.maxdistBETrevs: 		maximum distance from cluster centre location during each leave/revisit event. Comma separated by individual (with space) and revisit event (without space). See id.names/id.tags for order

#### `Points_With_Clusters.csv`

clusterID:			ID of cluster in which the location lies

animalID:			animalID of the animal to which the location belongs

timestamp.local:	local timezone timestamp of the location (format=YYYY-MM-DD HH:MM:SS)

date.local:			date of local timezone timestamp (format: YYYY-MM-DD)

time.local:			time of local timezone timestamp (format: HH:MM:SS)

local.timezone:		name of local timezone (incl. closest large city)

location.long:		longitude of location as measured by tag (GPS)

location.lat:		latitude of location as measured by tag (GPS)

ground.speed:		ground speed as measured by tag (GPS, unit: m/s, see you data set)

heading:			heading as measured by tag (GPS, unit: degrees, see you data set)

height.raw:			raw altitude as measured by tag (GPS, unit:m, see you data set)

clu.centr.long:		longitude of most central positon in cluster  to which this location belongs

clu.centr.lat:		latitutde of most central position in cluster to which this location belongs

clu.maxuse.long:	longitude of maximum used positon in cluster to which this location belongs

clu.maxuse.lat:		latitutde of maximum used position in cluster to which this location belongs

n.ids:				number of individuals using the cluster to which this location belongs

n.locs:				number of locations forming the cluster to which this location belongs

cumlag.locs:			cumulative timelag (duration to next logged location) of all locations in cluster to which this location belongs (unit as `dur_unit`; useful if input data are tracks with different resolutions; warning: data gaps can cause problemes here, always compare with n.locs

n.locsout:			number of locations outside the cluster to which this location belongs (restricted to full duration of cluser use and animals that used the cluster)

n.revs: 			number of revisits of cluster (i.e. number of times individuals leave cluster and come back, summed over all individuals that use the cluster


### Settings 
**Cluster Method (`meth`):** Method to cluster points, either `buff` or `hclust`. See details above in Documentation. Default `buff`.

**Cluster radius (`rad`):** Radius within which locations have to lie to be defined as a cluster. Unit = metre. Default 200 m.

**Minimum duration at a cluster (`dur`):** Duration that a cluster has to be repeatadly visited. Unit below. Default 1.

**Minimum number of locations building a cluster (`minloc`):** Minimum number of locations that a cluster needs to be made up of to be included into the results. Default 1.

**Time duration unit (`dur_unit`):** Duration unit for variable `dur`. Can be `hours`, `days` or `weeks`. Default `days`.

**Maximim gap duration between locations in a cluster (`maxgap`):** Maximum tolerated usage gap duration between successive locations in a cluster. Individuals to who the loations belong are irrelevant. Clusters are split at gaps exceeding this maximum gap duration. Unit below. Default 1.

**Gap duration unit (`gap_unit`):** Duration unit for variable `maxgap`. Can be `hours`, `days` or `weeks`. Default `days`.

**Cluster Transmission Mode (`clu_transm`):** Setting to allow only returning 'new' clusters that were not detected with the data set excluding the previous Y hours. Defaults to returning all clusters.

**Y - number of hours representing 'new' (`new_dur`):** Number of hours that data set shall be compared with to extract only 'new' clusters. Defaults to 24 (h).

**Locations of nests, roost or other sites to remove. (`remo_sites`):** A csv-file with locations, close to which clusters shall be excluded. Defaults to no file.

### Most common errors
Please post an issues [here](https://github.com/movestore/Point-Cluster-Detection/issues) if you encounter recurring errors or problems.

### Null or error handling:
**Setting `meth`:** Default `buff` allows no NULL.

**Setting `rad`:** Radius NULL defaults to 200 m. Too small radii might lead to small clusters, please include location inaccuracies here.

**Setting `dur`:** Duration NULL defaults to 1 (day). Too large durations might lead to few clusters.

**Setting `minloc`:** Number NULL defaults to 1. Too large values might lead to few clusters.

**Setting `dur_unit`:** Duration defaults to `days`. Only regular time units can be used (see above).

**Setting `maxgap`:** Duration NULL defaults to 1 (day). Too small maximum gaps might lead to fewer and smaller clusters.

**Setting `gap_unit`:** Duration defaults to `days`. Only regular time units can be used (see above).

**Setting `new_dur`:** If this number is too small, most clusters will be removed.

**Data:** All locations that are in a (any) cluster are returned to the next App. If no clusters are found in your data set NULL is returned, likel with an error.
