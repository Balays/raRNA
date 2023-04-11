
##### Worklfow PART 2.)

### Evaluate mapping and make plots
#source('bam.prep.all.R')

######  2. Cluster reads   ####
### Cluster the reads on exact matching in the dataframe of alignments (bam.all)
source('cluster.reads.R')

## Merge 'Transcripts' with 'Exons'
all.data <- merge(tr.gt, ex.sp[,c("EX_ID", by, "strand")],
                  by.x=c("start",  "end", "strand"),
                  by.y=c("start",  "end", "strand"), all=T)
## Merge with reads
if (is.lortia) {
  all.data <- merge(all.data, unique.data.frame(bam.filt[,c('qname', 'sample', 'correct_tes', 'correct_tss')]), by='qname', all=T)
} else {
  all.data <- merge(all.data, unique.data.frame(bam.filt[,c('qname', 'sample')]), by='qname', all=T)
}
all.data <- all.data[!is.na(all.data$qname), ]
## Merge with metadata
all.data <- merge(all.data, metafilt[,metacols], by.x='sample', by.y=metacols[1], all.x=T)
stopifnot(nrow(all.data) == nrow(bam.filt))
#### ####


######  3. Feature annotation ####
all.data$tr.ORF <- NA
TR.EX    <- unique.data.frame(all.data[,c("EX_ID", "TR_ID",by,"strand")])
if (annotate.transcripts.w.features) {
  source('feature.annotation.R')
}
stopifnot(nrow(all.data) == nrow(bam.filt))

######  4. Detection of Leaders and Trailers, Differentiating of sub-genomic RNAs and genomic RNAs  ####
## !!! This section is not finished! Many chunks are from SARS-CoV anaylsis and are not required!
## !!! Have to correct!
## Parameters for leader and trailer identification
leader.thresh  <- c(85,55) ## Leaders will be considered as mapped regions ('exons') that end in this range
trailer.thresh <- 29749 ## Trailers will be considered if the mapped region ends after this position.

source('leaders.and.trailers.R')
stopifnot(nrow(all.data) == nrow(bam.filt))
#### ####
##


######  6. Sum data ####
####
source('sum.data.R')

prime5s <- tr.all %>% group_by() %>% summarise(count=n())

#### ####
##


######  5. Associating reads to TSS clusters  ####
##
if (is.lortia) {
  TSS.all   <- unique.data.frame(all.data[,c("qname", "correct_tes", "correct_tss")])
  TSS.all   <- inner_join(TSS.all, tr.sp[,1:9], by='qname')
  TSS.all   <- inner_join(tr.uni[,1:4], TSS.all, by=c('TR_ID'))
} else {
  TSS.all   <- inner_join(tr.uni[,1:4], tr.sp[,1:9], by=c('TR_ID'))
}
colnames(TSS.all)[2:3] <- c('start', 'end')
TSS.all <- add.primes(TSS.all)

stopifnot((nrow(bam.filt) -length(dup(bam.filt$qname))) == nrow(TSS.all))


source('tss.assoc.R')

#### ####
##

## NOTE::
## This contains the transcript annotation and feature association results



######  7. Generate 5-prime end plots ####
####
if(make.plots) {
  combine.groups <- F
  source('cov.plots.R')
}

#### ####
##



######  7. Write outputs ####
####
write_tsv(read.end.count, paste0(outdir,'/read.end.count.tsv'))
write_tsv(prime5.sp, paste0(outdir,'/prime5.count.tsv'))
write_tsv(prime3.sp, paste0(outdir,'/prime3.count.tsv'))
write_tsv(TSS.sum,          paste0(outdir,'/TSS.counts.tsv'))


write_tsv(tr.all,           paste0(outdir,'/Transcripts.counts.tsv'))
#write_tsv(tr.gt,            paste0(outdir,'/Transcripts.Long.tsv'))
write_tsv(tr.sum,           paste0(outdir,'/Transcripts.tsv'))
write_tsv(tr.sum.cluster,   paste0(outdir,'/Transcripts.clusters.tsv'))
write_tsv(tr.sum.leader,    paste0(outdir,'/Transcripts.leader.tsv'))
write_tsv(tr.sum.genomic,   paste0(outdir,'/Transcripts.genomic.tsv'))
write_tsv(tr.sum.subgenomic,paste0(outdir,'/Transcripts.subgenomic.tsv'))
write_tsv(tr.sum.cluster.genomic, paste0(outdir,'/Transcripts.cluster.genomic.tsv'))
write_tsv(tr.uni,           paste0(outdir,'/Uni.Transcripts.tsv'))
write_tsv(tr.sp,            paste0(outdir,'/Reads.tsv'))
write_tsv(all.data,         paste0(outdir,'/Exons.tsv'))
write_tsv(orf.sum,          paste0(outdir,'/ORFs.tsv'))
write_tsv(ex.all.sum,       paste0(outdir,'/Exons.uniq.tsv'))
write_tsv(ex.sp,            paste0(outdir,'/Exons.uniq.sp.tsv'))
write_tsv(read.counts,      paste0(outdir,'/read.counts.tsv'))
#write_tsv(genomic.tr.df,    paste0(outdir,'/genomic.tr.df.tsv'))

#### ####
##
