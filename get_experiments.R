library(jsonlite)
library(httr)
library(gplots)

################################################################################################################
consumer_key = "QBJZ2HWO";
consumer_secret = "gb2prbadp2vuwucd";
secret <- RCurl::base64(paste(consumer_key, consumer_secret, sep=":"));
################################################################################################################

### This is where we obtain most of the metadata
URLs <- list("released_hg19"    = "type=experiment&status=released&assembly=hg19",
             "released_GRCh38"    = "type=experiment&status=released&assembly=GRCh38",
             # "preliminary" = "type=experiment&status=preliminary&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens",
             # "proposed"    = "type=experiment&status=proposed&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens",
             # "submitted"   = "type=experiment&status=submitted&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens",
             "all"         = "type=experiment&replicates.library.biosample.donor.organism.scientific_name=Homo%20sapiens");

for (i in 1:length(URLs)) {
    url <- paste("https://www.encodeproject.org/search/?", URLs[[i]], "&format=json&limit=all", sep="");
    nam <- names(URLs)[i];
    print(nam);
    if (length(grep(paste("^experiments_", nam, ".RData$", sep=""), dir())) == 0) {
        call <- GET(url, config(c("Authorization" = paste("Basic", secret))))
        obj <- fromJSON(rawToChar(call$content))

        # Parse some of the data and save.
        assays <- obj$facets[3,"terms"][[1]]
        dates <- obj$facets[22,"terms"][[1]]
        filetypes <- obj$facets[15,"terms"][[1]];
        experiments <- obj[["@graph"]]
        save(assays, dates, filetypes, experiments, file=paste("experiments_", nam, ".RData", sep=""))
    }
}

