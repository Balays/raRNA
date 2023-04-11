##### Import libraries and functions
library(cowplot)
library(grid)
library(ggplotify)
library(hrbrthemes, quietly = T)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(gggenes)
library(tidygenomics)
library(seqinr)
library(Rsamtools)
library(tidyverse)
library(stringi)
library(ggpubr)
library(tidygenomics)
library(ggsci)
library(tidyr)
library(dplyr)
library(fuzzyjoin)

#### Settings
by <- c("seqnames", "start",  "end")

## Important !!
fig.dir <- outdir #'figures.all.cov'
dir.create(fig.dir)
##

## colors and sizes
palette <-  c(pal_npg()(10), pal_aaas()(10))[c(1:10,14,18,15,13:11,16,17,20)]
alpha <- data.frame(cov_geom        = 1,
                    gene_geom       = 0.9,
                    unstranded_geom = 0.7)

vis.add <- 0 # c(250,250)

sizes <- 10
gene.sizes <- data.frame(gene_arrowhead_height=8,
                         gene_arrow_body_height=6,
                         gene_label_height=5,
                         gene_feature_height=6.0,
                         gene_feature_width=6.0,
                         gene_feature_label_height=5.0,
                         gene_feature_label_text=6.0,
                         unstranded_rect_height=0.5,
                         unstranded_feature_height=6,
                         unstranded_feature_width=6,
                         unstranded_feature_label_height=5,
                         unstranded_feature_label_text=6,
                         genome_feature_arrowhead_height=2.5,
                         genome_feature_arrow_body_height=2.5
                           )
## scale of annotation to coverage
genomplot.scale <- 6

#### Genes
genome.plotdata <- make.genome.plot.data(feature.df, 6, 'gene')
genome.plotdata$strand <- factor(genome.plotdata$strand, levels = c('+', '-', '*'))
gene.plotdata <- genome.plotdata
##
make.fromto <- function(l=l_genome, n=4) {

  fromto    <- data.frame(from=round(seq(1, l, l/n),0))
  fromto$to <- c(fromto$from[-1], l)
  return(fromto)
}

fromto <- make.fromto(n=4)
##
add.unstranded  <- T ### Correct!!
##

#### CAGEFightR results (by Adam Fulop)
cagefr       <- as.data.frame(rtracklayer::import.gff(paste0(virus, '_CAGE.gff3')))
cagefr.clust <- cagefr[cagefr$type == 'Cluster', c("seqnames", "start", "end", "strand", "type", "ID")]
cagefr.clust <- make.genome.plot.data(cagefr.clust, 6, 'gene')
cagefr.clust$strand <- factor(cagefr.clust$strand, levels = c('+', '-', '*'))

#### 5-prime end data
plot.data        <- TSS.all
plot.data$strand <- factor(plot.data$strand, levels = c('+', '-', '*'))
##
crop.FALSE      <- T
##
if (!is.lortia) {
  plot.data$correct_tes <- T
  plot.data$correct_tss <- T
}
##
if (combine.groups != F) {
  plot.data$hpi <- plot.data[,combine.groups]
} else {
  plot.data$hpi <- genome
}
##

#### Whole genome plots ####
breakseq <- 5000
visfrom  <- 1; visto=l_genome
## Barplot
message('Generating barplot for the whole genome ...')


###################### Plots for article

## EHV-1
if(virus == 'EHV-1') {
  gene.sizes$gene_label_height <- 3
  message('Generating plots for the Ori regions of ', virus, '...')

  fromto    <- ori.regions[1:2,]
  fig.names <- c('Figure 5A - EHV-1 OriL', 'Figure 5B - EHV-1 OriS internal repeat')

  colnames(fromto)[2:3] <- c('from', 'to')
  #for (i in 1:nrow(fromto) ) {

  ##
  i <- 1
  visfrom  <- fromto$from[i]
  visto    <- fromto$to[i]
  fig.name <- fig.names[i]
  ## when (i==1) {gene.sizes[,-c(4:6)] <- gene.sizes[,-c(4:6)] * 1.5 }
  message('Generating plot from ', visfrom, ' to ', visto)

  Fig5A <- plot.genome.region(visfrom=visfrom, visto=visto,
                                   plot.data=plot.data,
                                   gene.plotdata,
                                   prime = 'prime5',
                                   gene.label=F, add.unstranded=T, add.feature=F,
                                   add.cageTSS=T, cagefr.clust=cagefr.clust,
                                   crop.FALSE = crop.FALSE, scales = 'free_y',
                                   #geom   = geom_area(aes(y=count, x=prime5, fill=strand, color=strand), data=plot.sum, alpha=0.6),
                                   geom   = geom_density(aes(y=after_stat(density)), adjust = 0.025, n=512, color='black', trim=F, alpha=0.7), #abs(visto-visfrom)
                                   #samples='hpi12', #gggenome=gggenome,
                                   palette=palette, alpha=alpha, genomplot.scale=5, #genomplot.scale,
                                   gene.sizes = gene.sizes, sizes = sizes/2,
                                   y.multip = 2,
                                   force.gene.y = F, cagefr.clust.dist=3,
                                   gene.label.col='black',
                                   legend.position='none', plot.title = 'A',
                                   breakseq=breakseq, #ylim=c(0, 0.002),
                                   margins=unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
                                   gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene) )

  ggsave(paste0(fig.dir, '/', fig.name, '.jpg'), plot = Fig5A, width=6, height=9, limitsize = FALSE)

  ##
  i <- 2
  visfrom  <- fromto$from[i]
  visto    <- fromto$to[i]
  fig.name <- fig.names[i]
  ## when (i==1) {gene.sizes[,-c(4:6)] <- gene.sizes[,-c(4:6)] * 1.5 }
  message('Generating plot from ', visfrom, ' to ', visto)

  Fig5B <- plot.genome.region(visfrom=visfrom, visto=visto,
                                   plot.data=plot.data,
                                   gene.plotdata,
                                   prime = 'prime5',
                                   gene.label=F, add.unstranded=T, add.feature=F,
                                   add.cageTSS=T, cagefr.clust=cagefr.clust,
                                   crop.FALSE = crop.FALSE, scales = 'free_y',
                                   #geom   = geom_area(aes(y=count, x=prime5, fill=strand, color=strand), data=plot.sum, alpha=0.6),
                                   geom   = geom_density(aes(y=after_stat(density)), adjust = 0.025, n=512, color='black', trim=F, alpha=0.7), #abs(visto-visfrom)
                                   #samples='hpi12', #gggenome=gggenome,
                                   palette=palette, alpha=alpha, genomplot.scale=5, #genomplot.scale,
                                   gene.sizes = gene.sizes, sizes = sizes/2,
                                   y.multip = 2,
                                   force.gene.y = F, cagefr.clust.dist=3,
                                   gene.label.col='black',
                                   legend.position='none', plot.title = 'B',
                                   breakseq=breakseq, #ylim=c(0, 0.005),
                                   margins=unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
                                   gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene) )

  ggsave(paste0(fig.dir, '/', fig.name, '.jpg'), plot = Fig5B, width=6, height=9, limitsize = FALSE)

  cccccccccc <- cowplot::plot_grid(Fig5A ,
                                   Fig5B ,
                                   nrow = 1, rel_widths = c(1,1), align = 'h')

  ggsave(paste0(fig.dir, '/', 'Fig5_upper.jpg'), plot = Fig5_upper, width=14, height=9, limitsize = FALSE)

  #}

}
##

## KSHV
if(virus == 'KSHV') {
  gene.sizes$gene_label_height <- 3
  message('Generating plots for the Ori regions of ', virus, '...')

  fromto    <- ori.regions[1:2,]
  fig.names <- c('Figure 5C - KHSV OriLyt-R', 'Figure 5D - KHSV OriLyt-L', 'Figure 5E - KSHV OriLyt-L - 0.005 limit')

  colnames(fromto)[2:3] <- c('from', 'to')


  ## Fig5C
  i <- 2
  visfrom  <- fromto$from[i]
  visto    <- fromto$to[i]

  fig.name <- fig.names[1]

  message('Generating plot from ', visfrom, ' to ', visto)

  ##
  Fig5C <- plot.genome.region(visfrom=visfrom, visto=visto,
                                   plot.data=plot.data,
                                   gene.plotdata,
                                   prime = 'prime5',
                                   gene.label=F, add.unstranded=T, add.feature=F,
                                   add.cageTSS=T, cagefr.clust=cagefr.clust,
                                   crop.FALSE = crop.FALSE, scales = 'free_y',
                                   #geom   = geom_area(aes(y=count, x=prime5, fill=strand, color=strand), data=plot.sum, alpha=0.6),
                                   geom   = geom_density(aes(y=after_stat(density)), adjust = 0.025, n=512, color='black', trim=F, alpha=0.7), #abs(visto-visfrom)
                                   #samples='hpi12', #gggenome=gggenome,
                                   palette=palette, alpha=alpha, genomplot.scale=5, #genomplot.scale,
                                   gene.sizes = gene.sizes, sizes = sizes/2,
                                   y.multip = 2,
                                   force.gene.y = T, cagefr.clust.dist=3,
                                   gene.label.col='black',
                                   legend.position='none', plot.title = 'C',
                                   breakseq=breakseq, #ylim=c(0, 0.002),
                                   margins=unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
                                   gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene) )

  ggsave(paste0(fig.dir, '/', fig.name, '.jpg'), plot = Fig5C, width=6, height=9, limitsize = FALSE)

  ## Fig5D
  i <- 1
  visfrom  <- fromto$from[i]
  visto    <- fromto$to[i]

  fig.name <- fig.names[2]

  message('Generating plot from ', visfrom, ' to ', visto)

  ##
  Fig5D <- plot.genome.region(visfrom=visfrom, visto=visto,
                                   plot.data=plot.data,
                                   gene.plotdata,
                                   prime = 'prime5',
                                   gene.label=F, add.unstranded=T, add.feature=F,
                                   add.cageTSS=T, cagefr.clust=cagefr.clust,
                                   crop.FALSE = crop.FALSE, scales = 'free_y',
                                   #geom   = geom_area(aes(y=count, x=prime5, fill=strand, color=strand), data=plot.sum, alpha=0.6),
                                   geom   = geom_density(aes(y=after_stat(density)), adjust = 0.025, n=512, color='black', trim=F, alpha=0.7), #abs(visto-visfrom)
                                   #samples='hpi12', #gggenome=gggenome,
                                   palette=palette, alpha=alpha, genomplot.scale=5, #genomplot.scale,
                                   gene.sizes = gene.sizes, sizes = sizes/2,
                                   y.multip = 2,
                                   force.gene.y = T, cagefr.clust.dist=3,
                                   gene.label.col='black',
                                   legend.position='none', plot.title = 'D',
                                   breakseq=breakseq, #ylim=c(0, 0.002),
                                   margins=unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
                                   gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene) )

  ggsave(paste0(fig.dir, '/', fig.name, '.jpg'), plot = Fig5D, width=6, height=9, limitsize = FALSE)

  ## Fig5E
  i <- 1
  visfrom  <- fromto$from[i]
  visto    <- fromto$to[i]

  fig.name <- fig.names[3]

  Fig5E <- plot.genome.region(visfrom=visfrom, visto=visto,
                                   plot.data=plot.data,
                                   gene.plotdata,
                                   prime = 'prime5',
                                   gene.label=F, add.unstranded=T, add.feature=F,
                                   add.cageTSS=T, cagefr.clust=cagefr.clust,
                                   crop.FALSE = crop.FALSE, scales = 'free_y',
                                   #geom   = geom_area(aes(y=count, x=prime5, fill=strand, color=strand), data=plot.sum, alpha=0.6),
                                   geom   = geom_density(aes(y=after_stat(density)), adjust = 0.025, n=512, color='black', trim=F, alpha=0.7), #abs(visto-visfrom)
                                   #samples='hpi12', #gggenome=gggenome,
                                   palette=palette, alpha=alpha, genomplot.scale=5, #genomplot.scale,
                                   gene.sizes = gene.sizes, sizes = sizes/2,
                                   y.multip = 2,
                                   force.gene.y = T, cagefr.clust.dist=3,
                                   gene.label.col='black',
                                   legend.position='none', plot.title = 'E',
                                   breakseq=breakseq, ylim=c(0, 0.005),
                                   margins=unit(c(0.1, 0.5, 0.1, 0.5), "cm"),
                                   gene.aes = aes(xmin = start, xmax = end, y = ymin, fill = strand, forward = orientation, label = gene) )

  ggsave(paste0(fig.dir, '/', fig.name, '_ylim_0.002.jpg'), plot = Fig5E, width=6, height=9, limitsize = FALSE)

  Fig5_low <- cowplot::plot_grid(Fig5C ,
                                 Fig5D ,
                                 Fig5E ,
                                 nrow = 1, rel_widths = c(2,1,1), align = 'h')

  ggsave(paste0(fig.dir, '/', 'Fig5_low.jpg'), plot = Fig5_low, width=14, height=9, limitsize = FALSE)



}

##

