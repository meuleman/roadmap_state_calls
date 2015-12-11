Code used to process ChIP-seq histone tail modification data and call chromatin states, using existing Roadmap Epigenomics ChromHMM chromatin state models.

Largely based on code provided by Anshul Kundaje and the ENCODE3 (draft?) processing pipelines:
https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit#heading=h.9ecc41kilcvq

Processing steps:
1. Remove un/mis-mapped, low quality reads, mark & remove duplicates, compute library complexity QC, filter for mappability.
2. Pool reads across replicates and (if needed) sub-sample to 30M reads, use phantomPeakQualTools to obtain QC info.
3. Binarize data using ChromHMM BinarizeBed, create segmentation based on existing Roadmap models and output states per 200bp bins.

