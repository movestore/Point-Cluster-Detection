# Point Cluster Detection

MoveApps

Github repository: *github.com/movestore/Point-Cluster_Detection*

## Description
Detection of point clusters, where possibly more than one animal returns to within a specified time interval. Provides a table of each cluster with the times, duration, number of locations and animals. Clusters close to locations of ID "remove" can be excluded.

## Documentation
This App uses two types of clustering (buffer overlap or hierarchical clustering) for the detection of point clusters where one or more animals return to repeatedly within a specified time frame. 

For buffer clustering all locations are transformed into a circular spatial polygon of radius `cluster radius`. Then all overlapping polygons (i.e. of which the center locations are less than `2 * cluster radius` apart) are joined into clusters. Only clusters that were used for at least the specified number of hours/days/weeks are returned.

For hierarchical clustering the `average` method is used, i.e. the clusters are defined at the minimum average distance between all locations of the clusters. Clusters are selected to have at least a radius of `cluster radius` or be `2 * cluster radius` apart. Only clusters that were used for at least the specified number of hours/days/weeks are returned.

If one has uploaded a file of locations with individual.local.identifier="remove" in a preceding App, this App will automatically exclude those locations from the analysis and after the clustering exclude those clusters that have centre points less than `rad` metre from the locations of "remove". That way, fixed stations like nests or other points of attraction can be excluded from the results of this cluster analysis.

A cluster overview table is returned as a .csv artefact to download. It includes for each cluster the most central location (minimum distance to all other locations), timestamps of first and last location (UTC and local time), duration, cluster diameter, realised cluster radius (related to most central location), number of locations, number of animals, the names of those animals, their tag numbers and their respective duration and number of locations in the cluster.

The output of the App includes the locations that could be attributed to a cluster (that fulfilled the minimum duration requirement). This dataset is also returned as a .csv artefact to download, including local timestamps.

### Input data
moveStack in Movebank format

### Output data
moveStack in Movebank format


### Artefacts
`Cluster_Table.csv`: Overview of properties of detected point clusters.

`Points_With_Clusters.csv`: Result data set as .csv, with all locations in clusters.

### Parameters 
`meth`: Method to cluster points, either `buff` or `hclust`. See details above in Documentation. Default `buff`.

`rad`: Radius within which locations have to lie to be defined as a cluster. Unit = metre. Default 200 m.

`dur`: Duration that a cluster has to be repeatadly visited. Unit below. Default 14.

`dur_unit`: Duration unit for variable `dur`. Can be `hours`, `days` or `weeks`. Default `days`.

### Null or error handling:
**Parameter `meth`:** Default `buff` allows no NULL.

**Parameter `rad`:** Radius NULL defaults to 200 m. Too small radii might lead to small clusters, please include location inaccuracies here.

**Parameter `dur`:** Duration NULL defaults to 14 (days). Too large durations might lead to few clusters.

**Parameter `dur_unit`:** Duration defaults to `days`. Only regular time units can be used (see above).

**Data:** All locations that are in a (any) cluster are returned to the next App. If no clusters are found in your data set NULL is returned, likel with an error.
