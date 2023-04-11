##### Import libraries and functions
library(hrbrthemes)
library(ggsci)
library(seqinr)
library(Rsamtools)
library(tidyverse)
library(stringi)
library(ggpubr)
library(tidygenomics)
library(tidyr)
library(dplyr)
library(fuzzyjoin)

## download the functions in  from https://github.com/Balays/Rlyeh and set its directory below
rlyeh_dir <- '../Rlyeh-main'


source(paste0(rlyeh_dir, '/', 'cov.from.bam.R'))
source(paste0(rlyeh_dir, '/', 'cigar.sum.R'))
source(paste0(rlyeh_dir, '/', 'misc.functions.R'))
source(paste0(rlyeh_dir, '/', 'ov.from.bam2.R'))
source(paste0(rlyeh_dir, '/', 'get.best.aln.R'))
source(paste0(rlyeh_dir, '/', 'filter.and.import.bams.R'))
source(paste0(rlyeh_dir, '/', 'feature.OV.from.polyC.TR.R'))
source(paste0(rlyeh_dir, '/', 'idxstats.bams.R'))
source(paste0(rlyeh_dir, '/', 'mapping.eval.R'))
source(paste0(rlyeh_dir, '/', 'genome_plot.data.R'))
source(paste0(rlyeh_dir, '/', 'plot.genome.region.R'))


bam.flags <- read.delim(paste0(rlyeh_dir, '/', 'bam.flags.tsv'))


##### Settings for .bam files  ####
## Directory of .bam files
bamdir          <- 'bam'
## Directory of filtered .bam files
outfilt         <- 'bam.filt'; dir.create(outfilt)
pattern   <- '.bam'
bamfiles <- list.files(bamdir, pattern = pattern, recursive = T, full.names = T)
bamfiles <- bamfiles[grep('.bai', bamfiles, invert = T)]

##### Metadata ####

metadata         <- data.frame(sample   =gsub('.*/', '', gsub(pattern, '', bamfiles)))
metadata$virus   <- gsub('_.*', '', metadata$sample)
metadata$rep     <- gsub('_', '', stringi::stri_extract_first_regex(metadata$sample, '_._'))
metadata$date    <- gsub('^[ABC]_', '', stringi::stri_extract_first_regex(metadata$sample, '[ABC]_.*_'))
metadata$date    <- gsub('_', '', metadata$date)
metadata$barcode <- gsub('.*_', '', metadata$sample)

metafilt     <- metadata

## which columns of the metadata table to be included in the downstream analysis?
## The first element of this vector should be sample ID, which should be the same as the .bam file's names.
metacols     <- colnames(metadata) # c('sample', 'sample_ID')

#### ####
##

#### Settings ####
## Directory for output tables and figures
outdir <- paste0('WF.out')
dir.create(outdir)


### Miscallenaous settings
by  <- c('seqnames', 'start', 'end')
palette <- pal_npg()(10)
type <- 'start'
maxgap <- 10
minoverlap <- 1

window_size <- 1
window_step <- 1
whole.aln   <- F

make.plots <- T
combine.groups <- F
#### ####


######  Import .bam files      ####
## load already saved data?
load <- F
## filtered or unfiltered bam files?
filtered <- F
##
if (filtered) {
  save.data <- 'CAGE.filt.RData'
  if (load) {message('loading ', save.data, '...'); load(save.data)} else {
    source('import.filt.bams.R')
  }
} else {
  save.data <- 'CAGE.all.RData'
  if (load) {message('loading ', save.data, '...'); load(save.data)} else {
    source('import.all.bams.R')
  }
}
#length(unique(bam.filt$qname));length(unique(bam.prep.filt$qname))


### START Workflow part 1  #####

genomes <- as.character(na.omit(unique(bam.all$seqnames))); genomes <- genomes[order(genomes)]
viruses <- c('KSHV', 'EHV-1')
ori.regions <- read.delim('Ori_regions.tsv')

#### Add mapping width
bam.all$width <- abs(bam.all$start - bam.all$end)+1
bam.all$seqnames  <- as.character(bam.all$seqnames)

bam.filt$width    <- abs(bam.filt$start - bam.filt$end)+1
bam.filt$seqnames <- as.character(bam.filt$seqnames)


#### Filter to viral reads, mapping evaluation
bam.prep <- bam.all[!is.na(bam.all$seqnames), ]
stopifnot(nrow(bam.prep) == nrow(bam.filt))

bam.prep$org[bam.prep$seqnames == genomes[1]] <- viruses[1]
bam.prep$org[bam.prep$seqnames == genomes[2]] <- viruses[2]
bam.prep$org[bam.prep$seqnames == genomes[3]] <- viruses[3]
plyr::count(bam.prep$org)

source('bam.prep.R')

source('bam.prep.all.R')

#### Subsequent process .bam files?
filter.bams <- T
if (filter.bams) {
  ## give alignment IDs _1 and _2 for secondary alignments!!!
  #bam.seq   <- bam.prep[bam.prep$seqnames == genome, ]
  bam.filt  <- get.best.aln(bam.prep,#[bam.prep$sample == 'dRNA',],
                            bam.flags, best.mapq = F, rm.supplementary = F, rm.secondary = F, keep.chim = F,
                            give.qname.ID.to.secondary = T)

  ## switch aln_ID and qname columns
  cnames <- colnames(bam.filt)
  bam.filt <- data.frame(dplyr::select(bam.filt, -c('qname', 'aln_ID')),
                         qname=bam.filt$aln_ID,
                         aln_ID=bam.filt$qname
                         )
  bam.filt <- bam.filt[,cnames]

} else {
  bam.filt  <- bam.prep
}
## NOTE::
## in bam.filt, the alignments now have unique identifiers.
## bam.prep contains the original alignments and IDs
bam.filt$aln_count[is.na(bam.filt$aln_count)] <- 1

write_tsv(bam.filt, 'all.viral.reads.tsv')

rm(bam.all)
gc()

save.image('CAGE.filt.RData')
##### End of WF part1 #####


##### WF part 2 for each virus  #####
rm(list = ls(all=T))
##

### KSHV
load('CAGE.filt.RData')
bam.all <- bam.filt
##
i <- 1
genome <- genomes[i]; virus <- viruses[i]
outdir <- 'WF.out'
outdir <- paste0(virus, '_', outdir); dir.create(outdir)
fasta.ref <- paste0(genome, '.fasta')
fasta <- read.fasta(fasta.ref); l_genome <- length(fasta[[1]])
bam.filt <- bam.all[bam.all$org == viruses[i], ]

source('coverage.R')
write_tsv(cov.sum, paste0(outdir, '/KSHV.coverage.stats.tsv'))
kshv.mapped.cov <- mapped.cov
kshv.cov.sum <- cov.sum
kshv.win.cov.sum <- win.cov.sum

## annotation
gff          <- as.data.frame(rtracklayer::import.gff(paste0(genome,'.gff3')))
CDS.df       <- gff[gff$type == 'CDS', c("seqnames", "start", "end", "strand", "type", "Name", 'product')]
CDS.df$ID    <- CDS.df$Name
CDS.df$gene  <- CDS.df$product

# generate feature tables
feature.df <- CDS.df
feature.colname <- 'ID'
feature.df <- feature.df[, c(feature.colname, "strand", by, 'gene')]

for (i in unique(dup(feature.df[,feature.colname]))) {
  feature.df[feature.df[,feature.colname] == i, feature.colname] <-
    paste0(i, '_', 1:nrow(feature.df[feature.df[,feature.colname] == i, ]))
}
stopifnot(nrow(feature.df) == luniq(feature.df[,feature.colname]))

#write.table(feature.df, paste0(genome, '_feature.df.tsv'), sep = '\t', row.names = F, quote = F)
## ORi added manually !
feature.df <- read.delim(paste0(genome, '_feature.df.tsv'))

## 5-prime end data
TSS.all <- bam.filt[bam.filt$org == virus, ] %>% group_by(seqnames, qname, strand) %>% summarise(start=min(start), end=max(end))
TSS.all <- add.primes(TSS.all)

source('cov.plots.R')

save.image(paste0(virus, '_CAGE.RData'))
rm(list = ls(all=T))
gc()



#### EHV
load('CAGE.filt.RData')
bam.all <- bam.filt

i <- 2
genome <- genomes[i]; virus <- viruses[i]
outdir <- 'WF.out'
outdir <- paste0(virus, '_', outdir); dir.create(outdir)
fasta.ref <- paste0(genome, '.fasta')
fasta <- read.fasta(fasta.ref); l_genome <- length(fasta[[1]])
bam.filt <- bam.all[bam.all$org == viruses[i], ]

source('coverage.R')
write_tsv(cov.sum, paste0(outdir, '/EHV.coverage.stats.tsv'))
ehv.mapped.cov <- mapped.cov
ehv.cov.sum <- cov.sum
ehv.win.cov.sum <- win.cov.sum

## annotation
gff          <- as.data.frame(rtracklayer::import.gff(paste0(genome,'.gff3')))
CDS.df       <- gff[gff$type == 'CDS', c("seqnames", "start", "end", "strand", "type", "gene", 'Name')]
CDS.df$ID    <- CDS.df$Name


## generate feature tables
feature.df <- CDS.df
feature.colname <- 'ID'
feature.df <- feature.df[, c(feature.colname, "strand", by, 'gene')]

for (i in unique(dup(feature.df[,feature.colname]))) {
  feature.df[feature.df[,feature.colname] == i, feature.colname] <-
    paste0(i, '_', 1:nrow(feature.df[feature.df[,feature.colname] == i, ]))
}
stopifnot(nrow(feature.df) == luniq(feature.df[,feature.colname]))

#write.table(feature.df, paste0(genome, '_feature.df.tsv'), sep = '\t', row.names = F, quote = F)
## Ori regions added manually
feature.df <- read.delim(paste0(genome, '_feature.df.tsv'))

## 5-prime end data
TSS.all <- bam.filt[bam.filt$org == virus, ] %>% group_by(seqnames, qname, strand) %>% summarise(start=min(start), end=max(end))
TSS.all <- add.primes(TSS.all)

source('cov.plots.R')

save.image(paste0(virus, '_CAGE.RData'))
rm(list = ls(all=T))
gc()




### Figure 5
load('CAGE.filt.RData')
load(paste0(viruses[1], '_CAGE.RData'))
load(paste0(viruses[2], '_CAGE.RData'))


Fig5     <- cowplot::plot_grid(Fig5_upper ,
                               Fig5_low ,
                               nrow = 2, rel_heights = c(1,1))

ggsave(paste0(fig.dir, '/', 'Fig5.pdf'), plot = Fig5, width=14, height=9, limitsize = FALSE)

