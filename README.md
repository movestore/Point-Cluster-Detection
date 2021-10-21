# Point Cluster Detection

MoveApps

Github repository: *github.com/movestore/Point-Cluster_Detection*

## Description
Detection of point clusters, where possibly more than one animal returns to within a specified time interval. Provides a table of each cluster with the times, duration, number of locations and animals.

## Documentation
This App uses hierarchical clustering for the detection of point clusters where one or more animals return to repeatedly within a specified time frame. For clustering the Ward methods is used and the clusters are defined at the minimum distance between cluster centres of `2 * cluster radius`. Only clusters that were used for at least the specified number of hours/days/weeks are returned.

A cluster overview table is returned as a .csv artefact to download. It included for each cluster the mid location, timestamps of first and last location, duration, number of locations, number of animals and the names of those animals.

The output of the App includes only the locations that could be attributed to a cluster that fulfilled the minimum duration requirement. This dataset is also returned as a .csv artefact to download.

### Input data
moveStack in Movebank format

### Output data
moveStack in Movebank format


### Artefacts
`Cluster_Table.csv`: Overview of properties of detected point clusters.

`Points_With_Clusters.csv`: Result data set as .csv, with all locations in clusters.

### Parameters 
`rad`: Radius within which locations have to lie to be defined as a cluster. Unit = metre. Default 200 m.

`dur`: Duration that a cluster has to be repeatadly visited. Unit below. Default 14.

`dur_unit`: Duration unit for variable `dur`. Can be `hours`, `days` or `weeks`. Default `days`.

### Null or error handling:
**Parameter `rad`:** Radius NULL defaults to 200 m. Too small radii might lead to small clusters, please include location inaccuracies here.

**Parameter `dur`:** Duration NULL defaults to 14 (days). Too large durations might lead to few clusters.

**Parameter `dur_unit`:** Duration defaults to `days`. Only regular time units can be used (see above).

**Data:** All locations that are in a (any) cluster are returned to the next App. If no clusters are found in your data set NULL is returned, likel with an error.
