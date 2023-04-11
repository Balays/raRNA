

### Calculate coverage per bp for mapped regions
bam.cov         <- merge(bam.filt, metafilt, by='sample', all=F)
#bam.cov$sample  <- bam.cov$hpi
#bam.cov         <- bam.cov[bam.cov$seqnames == genome, ]

mapped.cov      <- cov.from.bam(bam.cov, samples = unique(bam.cov$sample), fasta.ref = fasta.ref)

### Summarise coverage per window
win.cov.sum <- window.cov(mapped.cov, fasta.ref,
                          window_size = 100, # window_size,
                          window_step = 100  # window_step
                          )

### Total coverage sum per window
win.cov.all.sum <- win.cov.sum %>% group_by(seqnames, window_id, start.window, end.window, strand) %>% summarise(all.sum=sum(window_sum.cov))

### Coverage summary
cov.sum     <- mapped.cov %>% group_by(sample) %>% summarise(mean=mean(count), median=median(count), sd=sd(count), sum=sum(count))

gc()

if (whole.aln) {
  ### Import bamfiles (rm.gaps.in.aln = F !!!)
  bam.reads <- import.bams(bamfiles, bamnames, write.filtered=write.filtered, force.create=force.create,
                           rm.gaps.in.aln = F, filtering=filtering, mapq.filt = mapq.filt,
                           flag.tokeep    = flag.tokeep,    flag.tocrop = flag.tocrop, seqnames.tofilt = seqnames.tofilt,
                           is.lortia = is.lortia)


  ### Calculate coverage per bp for whole alignments
  aln.cov <- cov.from.bam(bam.reads, samples = unique(bam.reads$sample), fasta.ref =fasta.ref)

}

