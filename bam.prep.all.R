
#### evaluate mapping in all samples
map.n_host <- spread(unique.data.frame(
  bam.prep %>% select(sample, qname, has.secondary, has.supplementary, aln.cat, sum.aln.width, width.ratio, n_seq_org, org)),
  org, n_seq_org, fill=0)
map.n_host$is.chim <- apply(map.n_host[,c(viruses)], 1, function(x) as.logical(length(x[x != 0])-1 ))

for (i in 1:length(viruses)) {
  map.n_host$org[map.n_host[,viruses[i]] == 1 ] <- viruses[i]
}

map.n_host$org[map.n_host$is.chim == T] <- 'chimaera'
plyr::count(map.n_host$org)


ggh <- gghistogram(map.n_host,
                   'width.ratio', fill='aln.cat', bins=100) +
  coord_cartesian() +
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10)) +
  facet_grid(rows = vars(org), cols=vars(aln.cat), scales = 'free_y') +
  ggtitle('Mapping ratio according to alignment category')

if(make.plots) {
  ggsave(paste0(outdir, '/mapping.ratio_per.aln.cat_density_all.reads_all_samples.jpg'), ggh, width = 12, height = 10)
}
gc()
##


