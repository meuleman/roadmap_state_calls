library(jsonlite)
library(httr)
library(gplots)

args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied.")
} else {        
    eval(parse(text=args[[1]])) # parse argument: type 
}

################################################################################################################

consumer_key = "QBJZ2HWO";
consumer_secret = "gb2prbadp2vuwucd";
secret <- RCurl::base64(paste(consumer_key, consumer_secret, sep=":"));

################################################################################################################

links_dir <- paste("file_links", type, sep="/");
dir.create(links_dir, showWarnings=FALSE, recursive=TRUE);

################################################################################################################
### Load experiment data

load(file=paste("experiments_", type, ".RData", sep=""))

##################################################################################################################

get_hrefs <- function(sel, filetype="bigWig") {
    hrefs <- data.frame();
    for (dataset in experiments[["@id"]][sel]) {
        url_files <- paste("https://www.encodeproject.org/search/?type=file&dataset=", 
                           dataset, "&file_format=", filetype, "&frame=object&format=json&limit=all", sep="");
        call <- GET(url_files, config(c("Authorization" = paste("Basic", secret))))
        obj_files <- fromJSON(rawToChar(call$content))

        files <- obj_files[["@graph"]]
        hrefs <- rbind(hrefs, data.frame(href=files$href, file_size=files$file_size));
    }
    invisible(hrefs);
}

assay <- "ChIP-seq";
sel <- which(experiments$assay_term_name == assay);
hrefs_bigWig <- tapply(sel, list(experiments$target$label[sel], experiments$biosample_term_name[sel]), 
                       function(x) {
                           print(paste(experiments$target$label[x[1]], experiments$biosample_term_name[x[1]], sep=" - "))
                           get_hrefs(x, filetype="bigWig");
                       });

num_cell_types <- colSums(apply(hrefs_bigWig, 1, sapply, function(x) sum(!is.null(x))))
ord <- order(num_cell_types, decreasing=TRUE);
hrefs_bigWig <- hrefs_bigWig[ord,];
save(hrefs_bigWig, file=paste(links_dir, "hrefs.RData", sep="/"));

###

for (i in 1:nrow(hrefs_bigWig)) {
    epitope <- rownames(hrefs_bigWig)[i];
    df <- data.frame();
    for (j in 1:ncol(hrefs_bigWig)) {
        cell_type <- gsub("\\/", "_", gsub("\ ", "_", colnames(hrefs_bigWig)[j]));
        if (!is.null(hrefs_bigWig[i,j][[1]])) {
            max_file <- hrefs_bigWig[i,j][[1]]$href[which.max(hrefs_bigWig[i,j][[1]]$file_size)];
            df <- rbind(df, data.frame(epitope=epitope, cell_type=cell_type, 
                                       file=paste("https://www.encodeproject.org", max_file, sep="/")));
        }
    }
    write.table(df, file=paste(links_dir, "/", epitope, ".csv", sep=""), 
                quote=FALSE, row.names=FALSE, sep=",");
}

###########################################################################################################

links_dir <- paste("file_links", type, sep="/");
load(file=paste("experiments_", type, ".RData", sep=""))

assay <- "DNase-seq";
sel <- which(experiments$assay_term_name == assay);
hrefs_bigWig <- tapply(sel, list(experiments$biosample_term_name[sel]), 
                       function(x) {
                           print(experiments$biosample_term_name[x[1]])
                           get_hrefs(x, filetype="bigWig");
                       });

save(hrefs_bigWig, file=paste(links_dir, "hrefs_DNase.RData", sep="/"));

###

epitope <- "DNase"
df <- data.frame();
max_files <- sapply(hrefs_bigWig, function(x) {
                    if (nrow(x) > 0) {
                        as.character(x$href[which.max(x$file_size)]);
                    }
                       });
if (length(max_files) > 0) {
    df <- data.frame(epitope=epitope, cell_type=gsub("\\/", "_", gsub("\ ", "_", names(max_files)[!sapply(max_files, is.null)])),
                     file=paste("https://www.encodeproject.org", max_files[!sapply(max_files, is.null)], sep="/"));
    write.table(df, file=paste(links_dir, "/", epitope, ".csv", sep=""), 
                quote=FALSE, row.names=FALSE, sep=",");
}



