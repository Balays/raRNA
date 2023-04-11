
## Alignment import
seqnames.tofilt <- genome ## Filter for alignments that mapped to this contig (viral genome)NA ##
flag.tokeep     <- NA # This will not drop alignments in .bam files based on bamflag only
#c(0,16) ## this will leave primary alignments only and filter out supplementary and secondary alignments as well
mapq.filt       <- 0 ## Mapping quality threshold (in minimap2 60 is the highest)
write.filtered  <- F ## Write out filtered .bam files
rm.gaps.in.aln  <- T
force.create    <- T ## Overwrite?
filtering       <- F ## filter out supplementary alignments from bamfiles
filter.bams     <- F ## subsequent filtering: get the best alignment for each read (to make qnames unique)
flag.tokeep     <- NA
flag.tocrop     <- NA

outfilt         <- NA

##
is.lortia       <- F

##


## file name pattern for .bam files (this will be cropped from the file name to get sample name)
pattern   <- '.bam'
bamfiles <- list.files(outfilt, pattern = pattern, recursive = T, full.names = T)
bamfiles <- bamfiles[grep('.bai', bamfiles, invert = T)]

nbam <- length(bamfiles)
bamnames <- gsub('.*\\/', '', bamfiles)
bamnames <- gsub(pattern, '', bamnames)

## Import
bam.filt  <- import.bams(bamfiles[],
                         bamnames[],
                         what = c("rname", "qname", "qwidth", "flag", "pos", "mapq", "cigar", "strand", "seq", "qual"),
                         write.filtered = F, force.create=force.create,
                         rm.gaps.in.aln = rm.gaps.in.aln, filtering=F, mapq.filt = mapq.filt,
                         flag.tokeep    = flag.tokeep,    flag.tocrop = flag.tocrop, seqnames.tofilt = seqnames.tofilt,
                         is.lortia = is.lortia)

bam.filt <- bam.filt[!is.na(bam.filt$qname), ]
if(!is.na(seqnames.tofilt)) {stopifnot(all(bam.filt$seqnames == seqnames.tofilt))
  bam.filt$seqnames <- as.character(seqnames.tofilt)}


#### Read counts (filtered)
reads           <- bam.filt[,] %>% group_by(sample, qname, qwidth) %>% summarise()
reads           <- merge(reads, metafilt[,metacols], by=metacols[1])
## This will not count the unmapped reads -> they will count as 1, as they do not have seqnames (which is fine), nor qname (I'm not sure why, though)
reads           <- reads[!is.na(reads$qname), ]
## So we have to calculate unmapped reads separately:
reads.unmapped  <- bam.filt[is.na(bam.filt$seqnames),] %>% group_by(sample) %>% summarise(unmapped_read_count=n())
## summarise
read.counts     <- as.data.frame(reads %>% group_by(across(any_of(c('seqnames', metacols)))) %>% summarise(mapped_read_count=n()))
## merge
read.counts     <- merge(read.counts, reads.unmapped, by='sample', all=T)
###

##
write_tsv(read.counts, paste0(outdir, '/read_counts.filt.tsv'))

read.counts.filt <- read.counts
read.counts.filt <- read.counts.filt[,c(1:4)]
colnames(read.counts.filt)[4] <- 'viral_nochim_read_count'
reads.filt <- reads
stopifnot(luniq(reads.filt$qname) == nrow(reads.filt))
rm(reads)
