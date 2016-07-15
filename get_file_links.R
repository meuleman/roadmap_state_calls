library(jsonlite)
library(httr)
library(gplots)
marks.only=FALSE

args=(commandArgs(TRUE))
if (length(args)==0) {
    stop("No arguments supplied.")
} else {        
    print(args)
    for (i in 1:length(args)){
        eval(parse(text=args[[i]])) # parse argument: type 
    }
}

if (marks.only){
    epitopes=c("Control","H3K4me1","H3K4me3","H3K27me3","H3K36me3","H3K9me3","H3K27ac")
    lbl = '_marks'
} else { lbl = '' }

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

get_hrefs <- function(sel, filetype="bam") {
    hrefs <- data.frame();
    for (dataset in experiments[["@id"]][sel]) {
        url_files <- paste("https://www.encodeproject.org/search/?type=file&dataset=", dataset, "&file_format=", filetype, "&frame=object&format=json&limit=all", sep="");
        call <- GET(url_files, config(c("Authorization" = paste("Basic", secret))))
        obj_files <- fromJSON(rawToChar(call$content))

        files <- obj_files[["@graph"]]
        hrefs <- rbind(hrefs, data.frame(href=files$href, file_size=files$file_size));
    }
    invisible(hrefs);
}

assay <- "ChIP-seq";
sel <- which(experiments$assay_term_name == assay);
exp.info <- experiments[sel,'target']
if (marks.only){ 
    sel <- sel[exp.info$label %in% epitopes]
}

hrefs_bam <- tapply(sel, list(experiments$target$label[sel], experiments$biosample_term_name[sel]), 
                    function(x) {
                        print(paste(experiments$target$label[x[1]], experiments$biosample_term_name[x[1]], sep=" - "))
                        get_hrefs(x, filetype="bam");
                    });

num_cell_types <- colSums(apply(hrefs_bam, 1, sapply, function(x) sum(!is.null(x))))
ord <- order(num_cell_types, decreasing=TRUE);
hrefs_bam <- hrefs_bam[ord,];
save(hrefs_bam, file=paste0(links_dir, "/hrefs",lbl,".RData"));

# Make availability matrix:
mat <- apply(hrefs_bam,1:2,function(x){length(unlist(x))/2})


for (i in 1:nrow(hrefs_bam)) {
    epitope <- rownames(hrefs_bam)[i];
    df <- data.frame();
    for (j in 1:ncol(hrefs_bam)) {
        cell_type <- gsub("\\/", "_", gsub("\ ", "_", colnames(hrefs_bam)[j]));
        if (!is.null(hrefs_bam[i,j][[1]])) {
            max_file <- hrefs_bam[i,j][[1]]$href[which.max(hrefs_bam[i,j][[1]]$file_size)];
            df <- rbind(df, data.frame(epitope=epitope, cell_type=cell_type, 
                                       file=paste0("https://www.encodeproject.org", max_file)));
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
hrefs_bam <- tapply(sel, list(experiments$biosample_term_name[sel]), 
                    function(x) {
                        print(experiments$biosample_term_name[x[1]])
                        get_hrefs(x, filetype="bam");
                    });

save(hrefs_bam, file=paste(links_dir, "hrefs_DNase.RData", sep="/"));

###
epitope <- "DNase"
df <- data.frame();
max_files <- sapply(hrefs_bam, function(x) {
                    if (nrow(x) > 0) {
                        as.character(x$href[which.max(x$file_size)]);
                    }
                    });

if (length(max_files) > 0) {
    df <- data.frame(epitope=epitope, cell_type=gsub("\\/", "_", gsub("\ ", "_", names(max_files)[!sapply(max_files, is.null)])),
                     file=paste0("https://www.encodeproject.org", max_files[!sapply(max_files, is.null)]));
    write.table(df, file=paste(links_dir, "/", epitope, ".csv", sep=""), 
                quote=FALSE, row.names=FALSE, sep=",");
}


# Merge matDN with mat (ChIP-seq data), make availability table: 
repDN <- unlist(lapply(hrefs_bam,length))
DNase <- rep(0,ncol(mat))
names(DNase) <- colnames(mat)
namDN <- colnames(mat)[colnames(mat) %in% names(repDN)]
DNase[namDN] <- repDN[namDN]
mat <-  rbind(mat,DNase)

write.table(mat, file=paste0(links_dir, "/replicates",lbl,".tsv"),quote=F,row.names=TRUE,col.names=TRUE,sep='\t')

# Find cell types with full coverage:
prod <- which(apply(mat[-8,],2,prod) != 0) # DNase coverage not necessary
avail <- names(prod)
write.table(avail, file=paste0(links_dir, "/available",lbl,".tsv"),quote=F,row.names=F,col.names=F);

# Plot available types
pmat <- mat
pmat[pmat > 0] <- 1
ord <- order(-colSums(pmat),colnames(pmat),decreasing=TRUE)
pmat <- pmat[,ord]

png(paste0(links_dir,'/plot_availability',lbl,'.png'),res=450,units='in',width=10,height=15)
par(mar = c(2,10,1.5,1))
image(pmat, axes=FALSE, col=c('white','darkblue'),main=paste('Availability for',type))
grid(nx=nrow(pmat), ny=ncol(pmat),col='grey',lty='solid',lwd=.25)
text(y=seq(0,1,length.out=ncol(pmat)), x=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=colnames(pmat), srt=0, adj=1, xpd=TRUE,cex=.6)
text(x=seq(0,1,length.out=nrow(pmat)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(pmat), srt=0, adj=.5, xpd=TRUE,cex=.75)
dev.off()
###


