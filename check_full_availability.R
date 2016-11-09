# Check for full data availability (all + submitted)
library(jsonlite)
library(httr)
library(gplots)
library(RColorBrewer)
colrwb <- colorRampPalette(c('darkred','white','royalblue'))(200)
colramp <- colorRampPalette(c('white','royalblue'))(100)
setwd('~/cHMM/db')

marks.only=FALSE
img <- '~/cHMM/img/'
type='all_submitted_released'

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

assay <- "ChIP-seq";
sel <- which(experiments$assay_term_name == assay);
exp.info <- experiments[sel,'target']
sample.info <- experiments[sel,'biosample_term_name']
exp.df <- experiments[sel,c(1:7,10,11)]



assays <- c('DNase-seq','ATAC-seq')
sel <- which(experiments$assay_term_name %in% assays);
dhs.df <- experiments[sel,c(1:7,10,11)]

all.df <- rbind(dhs.df,exp.df)
all.df$mark <- all.df$target$label                                          
all.df$mark[all.df$assay_title == 'DNase-seq'] <- 'DHS'
all.df$mark[all.df$assay_title == 'ATAC-seq'] <- 'ATAC'
# FLAG and eGFP
all.df$construct <- 'None'
fid <- grep('FLAG',all.df$mark)
eid <- grep('eGFP',all.df$mark)
all.df$construct[fid] <- 'FLAG'
all.df$construct[eid] <- 'eGFP'
all.df$mark[fid] <- sub('FLAG-','',all.df$mark[fid])
all.df$mark[eid] <- sub('eGFP-','',all.df$mark[eid])

# TALEs in K562?
cells <- sort(unique(all.df$biosample_term_name))
marks <- sort(unique(all.df$mark))

mat <- matrix(0,nrow=length(marks),ncol=length(cells),dimnames=list(marks,cells))
ids <- data.frame(mark=all.df$mark,cell=all.df$biosample_term_name,count=1)
idcount <- aggregate(count ~ mark + cell,ids,sum)
idm <- as.matrix(idcount[,1:2])
mat[idm] <- idmat[,3]

# Reorder matrix:
mat <- mat[,order(colSums(mat),decreasing=TRUE)]
mat <- mat[order(rowSums(mat),decreasing=FALSE),]

# ===================
# All available data:
# ===================
rows <- (seq(0,nrow(mat),20) + .5)/ (nrow(mat) -1)
skip <- ncol(mat) %/% 20
cols <- (seq(0,ncol(mat),skip) + .5)/ (ncol(mat) -1)

png(paste0(img,type,'_full.png'),res=450,width=17,height=11,units='in')
# pdf(paste0(img,type,'_full.pdf'),width=17,height=11)
par(mar = c(5,1,3,10))
par(xaxs="i")
par(cex.main=2)
image(log(mat),axes=F)
title(main="ENCODE available ChIP-seq + ATAC/DNase-seq", xlab="Dataset (all released in hg19 or GRCh38 and submitted)")
mtext('Cell Type', 4, line = 7)
grid(nx = nrow(mat), ny = ncol(mat),col='grey',lty='solid',lwd=.15)
grid(nx = nrow(mat), ny = ncol(mat),col='grey',lty='solid',lwd=.15)
text(x=seq(0,1,length.out=nrow(mat)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(mat), srt=90, adj=1, xpd=TRUE,cex=.25)
text(y=seq(0,1,length.out=ncol(mat)), x=par()$usr[4]+0.01*(par()$usr[4]-par()$usr[3]), labels=colnames(mat), srt=0, adj=0, xpd=TRUE,cex=.25)
abline(v=rows,col='black',lty='solid',lwd=.5)
abline(h=cols,col='black',lty='dashed',lwd=.5)
dev.off()

# ==========
# Core marks
# ==========
core <- rev(scan('core_marks_with_atac','c'))
core <- sub('\\.','F',core)
cmat <- mat[core,]
cmat <- cmat[,which(colSums(cmat[-nrow(cmat),]) > 0)]
cols <- (seq(0,ncol(cmat),skip) + .5)/ (ncol(cmat) -1)
cmat <- cmat[,order(colSums(cmat),decreasing=TRUE)]

pdf(paste0(img,type,'_core.pdf'),width=8.5,height=11)
# png(paste0(img,type,'_core.png'),res=450,width=8.5,height=11,units='in')
par(mar = c(5,1,3,8))
par(xaxs="i")
par(cex.main=1)
image(log(cmat),axes=F)
title(main="ENCODE available ChIP-seq + ATAC/DNase-seq", xlab="Core Dataset")
mtext(paste0('Cell Type (total = ',ncol(cmat),')'), 4, line = 5)
grid(nx = nrow(cmat), ny = ncol(cmat),col='grey',lty='solid',lwd=.15)
grid(nx = nrow(cmat), ny = ncol(cmat),col='grey',lty='solid',lwd=.15)
text(x=seq(0,1,length.out=nrow(cmat)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(cmat), srt=90, adj=1, xpd=TRUE,cex=.5)
text(y=seq(0,1,length.out=ncol(cmat)), x=par()$usr[4]+0.05*(par()$usr[4]-par()$usr[3]), labels=colnames(cmat), srt=0, adj=0, xpd=TRUE,cex=.25)
abline(h=cols,col='black',lty='dashed',lwd=.5)
dev.off()


# ========================================
# Without control or dhs/atac (Just marks)
# ========================================
jm <- marks[!(marks %in% c('Control','DHS','ATAC'))]
jmat <- mat[jm,]
jmat <- jmat[,which(colSums(jmat) > 0)]
jmat <- jmat[,order(colSums(jmat),decreasing=TRUE)]
jmat <- jmat[order(rowSums(jmat),decreasing=FALSE),]
cols <- (seq(0,ncol(jmat),skip) + .5)/ (ncol(jmat) -1)

# pdf(paste0(img,type,'_marks.pdf'),width=17,height=11)
png(paste0(img,type,'_marks.png'),res=450,width=17,height=11,units='in')
par(mar = c(5,1,3,10))
par(xaxs="i")
par(cex.main=1)
image(log(jmat),axes=F)
title(main="ENCODE available ChIP-seq + ATAC/DNase-seq", xlab="ChIP-seq Dataset")
mtext(paste0('Cell Type (total = ',ncol(jmat),')'), 4, line = 6)
grid(nx = nrow(jmat), ny = ncol(jmat),col='grey',lty='solid',lwd=.15)
grid(nx = nrow(jmat), ny = ncol(jmat),col='grey',lty='solid',lwd=.15)
text(x=seq(0,1,length.out=nrow(jmat)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(jmat), srt=90, adj=1, xpd=TRUE,cex=.25)
text(y=seq(0,1,length.out=ncol(jmat)), x=par()$usr[4]+0.01*(par()$usr[4]-par()$usr[3]), labels=colnames(jmat), srt=0, adj=0, xpd=TRUE,cex=.25)
abline(h=cols,col='black',lty='dashed',lwd=.5)
dev.off()

# ======================
# submitted vs. released
# ======================
rmat <- matrix(0,nrow=length(marks),ncol=length(cells),dimnames=list(marks,cells))
rids <- data.frame(mark=all.df$mark,cell=all.df$biosample_term_name,released=all.df$status=='released')
ridcount <- aggregate(released ~ mark + cell,rids,sum)
total <- merge(idcount,ridcount)
total$summary <- 2
total$summary[total$released == 0] <- 1
idm <- as.matrix(total[,1:2])
rmat[idm] <- total$summary

# Without control or dhs/atac (Just marks)
rmat <- rmat[,which(colSums(rmat) > 0)]
rmat <- rmat[,order(colSums(rmat),decreasing=TRUE)]
rmat <- rmat[order(rowSums(rmat),decreasing=FALSE),]
rows <- (seq(0,nrow(rmat),50) + .5)/ (nrow(rmat) -1)
cols <- (seq(0,ncol(rmat),skip) + .5)/ (ncol(rmat) -1)

pdf(paste0(img,type,'_rs.pdf'),width=17,height=11)
# png(paste0(img,type,'_rs.png'),res=450,width=17,height=11,units='in')
par(mar = c(5,1,3,10))
par(xaxs="i")
par(cex.main=1)
image(rmat,axes=F,col=c('white','blue','red'))
title(main="ENCODE available ChIP-seq + ATAC/DNase-seq", xlab="Has Released (red) vs Only Submitted (blue) Data")
mtext(paste0('Cell Type (total = ',ncol(rmat),')'), 4, line = 6)
grid(nx = nrow(rmat), ny = ncol(rmat),col='grey',lty='solid',lwd=.15)
grid(nx = nrow(rmat), ny = ncol(rmat),col='grey',lty='solid',lwd=.15)
text(x=seq(0,1,length.out=nrow(rmat)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(rmat), srt=90, adj=1, xpd=TRUE,cex=.25)
text(y=seq(0,1,length.out=ncol(rmat)), x=par()$usr[4]+0.01*(par()$usr[4]-par()$usr[3]), labels=colnames(rmat), srt=0, adj=0, xpd=TRUE,cex=.25)
abline(v=rows,col='black',lty='solid',lwd=.5)
abline(h=cols,col='black',lty='dashed',lwd=.5)
dev.off()



# ===========================
# submitted vs. release core:
# ===========================
crmat <- rmat[core,]
crmat <- crmat[,which(colSums(crmat[-nrow(crmat),]) > 0)]
cols <- (seq(0,ncol(crmat),skip) + .5)/ (ncol(crmat) -1)
crmat <- crmat[,order(colSums(crmat),decreasing=TRUE)]

pdf(paste0(img,type,'_core_rs.pdf'),width=8.5,height=11)
# png(paste0(img,type,'_core_rs.png'),res=450,width=8.5,height=11,units='in')
par(mar = c(5,1,3,8))
par(xaxs="i")
par(cex.main=1)
image(crmat,axes=F,col=c('white','royalblue','indianred3'))
title(main="ENCODE available ChIP-seq + ATAC/DNase-seq", xlab="Core Released (red) vs. Submitted Only (blue)")
mtext(paste0('Cell Type (total = ',ncol(crmat),')'), 4, line = 5)
grid(nx = nrow(crmat), ny = ncol(crmat),col='grey',lty='solid',lwd=.15)
grid(nx = nrow(crmat), ny = ncol(crmat),col='grey',lty='solid',lwd=.15)
text(x=seq(0,1,length.out=nrow(crmat)), y=par()$usr[3]-0.01*(par()$usr[4]-par()$usr[3]), labels=rownames(crmat), srt=90, adj=1, xpd=TRUE,cex=.5)
text(y=seq(0,1,length.out=ncol(crmat)), x=par()$usr[4]+0.05*(par()$usr[4]-par()$usr[3]), labels=colnames(crmat), srt=0, adj=0, xpd=TRUE,cex=.25)
abline(h=cols,col='black',lty='dashed',lwd=.5)
dev.off()



