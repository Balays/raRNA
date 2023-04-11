
## Alignment import settings
seqnames.tofilt <- NA ## genome ## Filter for alignments that mapped to this contig (viral genome)
flag.tokeep     <- NA # This will not drop alignments in .bam files based on bamflag only
#c(0,16) ## this will leave primary alignments only and filter out supplementary and secondary alignments as well
mapq.filt       <- 0 ## Mapping quality threshold (in minimap2 60 is the highest)
write.filtered  <- T ## Write out filtered .bam files
rm.gaps.in.aln  <- T
force.create    <- T ## Overwrite?
filtering       <- F ## filter out supplementary alignments from bamfiles
filter.bams     <- F ## subsequent filtering: get the best alignment for each read (to make qnames unique)
flag.tokeep     <- NA
flag.tocrop     <- NA


##
is.lortia       <- F

##

nbam <- length(bamfiles)
bamnames <- gsub('.*\\/', '', bamfiles)
bamnames <- gsub(pattern, '', bamnames)


#### Import bams
bam.all  <- import.bams(bamfiles,
                        bamnames,
                        write.filtered = F, force.create=force.create, #add.primes=F,
                        rm.gaps.in.aln = rm.gaps.in.aln, filtering= F, mapq.filt = mapq.filt,
                        flag.tokeep    = flag.tokeep,    flag.tocrop = flag.tocrop, seqnames.tofilt = NA, # genome, #
                        is.lortia = is.lortia)


all.stats <- idxstats.bams(bamfiles, bamnames = bamnames )

#### Summarise seqnames and samples
bam.sum <- bam.all %>% group_by(sample, seqnames) %>% summarise(sub_alignment_count=n())

for (i in unique(all.stats$sample)) {
  stats <- all.stats[all.stats$sample == i, ]
  try({ all.stats$unmapped[all.stats$sample == i] <- stats$unmapped[stats$seqnames == '*'] })
}

all.stats <- all.stats[all.stats$seqnames != '*', ]
all.stats$ratio <- all.stats$mapped / all.stats$unmapped
##

gc()
#### ####
if (write.filtered) {
  bam.filt <- import.bams(bamfiles[],
                          bamnames[],
                          ## for import.bams
                          write.filtered = T, force.create=force.create,
                          rm.gaps.in.aln = rm.gaps.in.aln, filtering='all.reads.w.supp.ali',
                          ## for ov.from.bam2
                          mapq.filt = mapq.filt,
                          flag.tokeep    = flag.tokeep,    flag.tocrop = flag.tocrop, seqnames.tofilt = seqnames.tofilt,
                          is.lortia = is.lortia)
  table(bam.filt[,c('seqnames', 'flag')])
  ##
  target <- gsub('\\.bam', '.filtered.bam', bamfiles)
  target <- c(target, paste0(target, '.bai'))
  dir.create(outfilt)
  destination     <- gsub(paste0(bamdir, '/'), paste0(outfilt, '/'), target)
  destination     <- gsub('.filtered', '', destination)
  file.rename(target, destination)

  bamfiles <- destination[grep('.bai', destination, invert = T)] #paste0(stringi::stri_replace_last_regex(bamfiles, '.bam', ''), '.filtered.bam')

  bam.filt <- bam.filt[!is.na(bam.filt$qname), ]
  if(!is.na(seqnames.tofilt)) {stopifnot(all(bam.filt$seqnames == seqnames.tofilt)); bam.filt$seqnames <- as.character(bam.filt$seqnames)}



  #### Summarise seqnames and samples
  filt.stats <- idxstats.bams(bamfiles, bamnames = bamnames )

  for (i in unique(filt.stats$sample)) {
    stats <- filt.stats[filt.stats$sample == i, ]
    try({ filt.stats$unmapped[filt.stats$sample == i] <- stats$unmapped[stats$seqnames == '*'] })
  }

  filt.stats <- filt.stats[filt.stats$seqnames != '*', ]
  filt.stats$ratio <- filt.stats$mapped / all.stats$unmapped
  ##
  bam.filt.sum <- bam.filt %>% group_by(sample, seqnames) %>% summarise(sub_alignment_count=n())
  ##
}




#### Read counts (all)
reads           <- bam.all[,] %>% group_by(sample, qname, qwidth) %>% summarise()
reads           <- merge(reads, metafilt[,metacols], by=metacols[1])
##
## This will not count the unmapped reads -> they will count as 1, as they do not have seqnames (which is fine), nor seqnames (I'm not sure why, though)
reads           <- reads[!is.na(reads$qname), ]
## So we have to calculate unmapped reads separately:
reads.unmapped  <- bam.all[is.na(bam.all$seqnames), ] %>% group_by(sample) %>% summarise(unmapped_read_count=n())
## summarise
read.counts     <- as.data.frame(reads %>% group_by(across(any_of(c('seqnames', metacols)))) %>% summarise(mapped_read_count=n()))
## merge
read.counts     <- merge(read.counts, reads.unmapped, by='sample', all=T)
## OK!
read.counts.all <- read.counts

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


## merge filt w all data
read.counts <- merge(read.counts.all, read.counts.filt, by=1:3, all=T)
write_tsv(read.counts, paste0(outdir, '/read_counts.tsv'))


##
write_tsv(read.counts, paste0(outdir, '/read_counts.all.tsv'))


