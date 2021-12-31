# Clean data objects

This folder contains, for each dataset, the following R .rds objects:

* **'*.SE.rds': ** A SummarizedExperiment object with annotated samples, raw and normalized quantifications
* **'*.DEA.rds': ** A data.frame containing the results of the various differential expression analyses relative to the object
* **'*.DEA.enrichments.rds': ** A list of data.frames containing the results of the various enrichment analyses performed on the differential expression analyses

The code to generate these objects can be found in the [indiv_processing](../indiv_processing/) folder.