# Point Cluster Detection

MoveApps

Github repository: *github.com/movestore/Point-Cluster_Detection*

## Description
Detection of point clusters, where possibly more than one animal returns to within a specified time interval. Provides a table of each cluster with the times, duration, number of locations and animals. Clusters close to locations of ID "remove" can be excluded.

## Documentation
This App uses two types of clustering (buffer overlap or hierarchical clustering) for the detection of point clusters where one or more animals return to repeatedly within a specified time frame. 

For buffer clustering all locations are transformed into a circular spatial polygon of radius `cluster radius`. Then all overlapping polygons (i.e. of which the center locations are less than `2 * cluster radius` apart) are joined into clusters. Only clusters that were used for at least the specified number of hours/days/weeks (dur parameter) are returned. Clusters with time gaps between successive (indiv. independent) locations above the user defined maximum (maxgap parameter) are split. Resulting partial clusters are kept if they fulfill the minimum duration requirement (dur parameter).

For hierarchical clustering the `average` method is used, i.e. the clusters are defined at the minimum average distance between all locations of the clusters. Clusters are selected to have at least a radius of `cluster radius` or be `2 * cluster radius` apart. Only clusters that were used for at least the specified number of hours/days/weeks are returned.

If one has uploaded a file of locations with individual.local.identifier="remove" in a preceding App, this App will automatically exclude those locations from the analysis and after the clustering exclude those clusters that have centre points less than `rad` metre from the locations of "remove". That way, fixed stations like nests or other points of attraction can be excluded from the results of this cluster analysis.

According to the maxgap parameter, clusters with usage gaps larger than the defined value are split.

A cluster overview table is returned as a .csv artefact to download. It includes for each cluster the most central location (minimum distance to all other locations), timestamps of first and last location (UTC and local time), duration, cluster diameter, realised cluster radius (related to most central location), number of locations, number of animals, the names of those animals, their tag numbers and their respective duration and number of locations in the cluster. See definition of all parameters in the list(s) below.

The output of the App includes the locations that could be attributed to a cluster (that fulfilled the minimum duration requirement), including the clusterID, number of animals using the cluster, number of locations in and outside the cluster (in the time that the cluster was used and by the individuals that used the clsuter) and number of revisits. This dataset is also returned as a .csv artefact to download, including local timestamps.

### Input data
moveStack in Movebank format

### Output data
moveStack in Movebank format


### Products (Artefacts)
`Cluster_Table.csv`: Overview of properties of detected point clusters.

`Points_With_Clusters.csv`: Result data set as .csv, with all locations in clusters.

### Table column definitions

#### `Cluster_Table.csv`

clusterID:				ID of cluster

centr.long:				longitude of most central positon in cluster (minimum average distance to all other locations)

centr.lat:				latitutde of most central position in cluster (minimum average distance to all other locations)

timestamp.start:			time of first location in cluster (UTC, format=YYYY-MM-DD HH:MM:SS)

timestamp.end:			time of last location in cluster (UTC, format=YYYY-MM-DD HH:MM:SS)

timestamp.start.local:	local timestamp of first location in cluster (format=YYYY-MM-DD HH:MM:SS)

timestamp.end.local:		local timestamp of last location in cluster (format=YYYY-MM-DD HH:MM:SS)

local.timezone:			name of local timezone (incl. closest large city)

duration (hours):		time duration between first and last location in cluster (unit=hours)

n.locs:					number of locations in cluster (independent of individual)

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

n.ids:				number of individuals using the cluster to which this location belongs

n.locs:				number of locations forming the cluster to which this location belongs

n.locsout:			number of locations outside the cluster to which this location belongs (restricted to full duration of cluser use and animals that used the cluster)

n.revs: 			number of revisits of cluster (i.e. number of times individuals leave cluster and come back, summed over all individuals that use the cluster


### Parameters 
`meth`: Method to cluster points, either `buff` or `hclust`. See details above in Documentation. Default `buff`.

`rad`: Radius within which locations have to lie to be defined as a cluster. Unit = metre. Default 200 m.

`dur`: Duration that a cluster has to be repeatadly visited. Unit below. Default 14.

`dur_unit`: Duration unit for variable `dur`. Can be `hours`, `days` or `weeks`. Default `days`.

`maxgap`: Maximum tolerated usage gap duration between successive locations in a cluster. Individuals to who the loations belong are irrelevant. Clusters are split at gaps exceeding this maximum gap duration. Unit below. Default 1.

`gap_unit`: Duration unit for variable `maxgap`. Can be `hours`, `days` or `weeks`. Default `days`.

### Null or error handling:
**Parameter `meth`:** Default `buff` allows no NULL.

**Parameter `rad`:** Radius NULL defaults to 200 m. Too small radii might lead to small clusters, please include location inaccuracies here.

**Parameter `dur`:** Duration NULL defaults to 14 (days). Too large durations might lead to few clusters.

**Parameter `dur_unit`:** Duration defaults to `days`. Only regular time units can be used (see above).

**Parameter `maxgap`:** Duration NULL defaults to 1 (day). Too small maximum gaps might lead to fewer and smaller clusters.

**Parameter `gap_unit`:** Duration defaults to `days`. Only regular time units can be used (see above).

**Data:** All locations that are in a (any) cluster are returned to the next App. If no clusters are found in your data set NULL is returned, likel with an error.
