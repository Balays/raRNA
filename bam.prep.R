

nrow_bam <- nrow(bam.prep)

bam.prep <- calc.alignment.ratio(bam.prep)
if (nrow(bam.prep) != nrow_bam) { bam.prep <- unique.data.frame(bam.prep) } else NULL

bam.prep <- calc.n_seq(bam.prep)
if (nrow(bam.prep) != nrow_bam) { bam.prep <- unique.data.frame(bam.prep) } else NULL

bam.prep <- calc.n_org(bam.prep)
if (nrow(bam.prep) != nrow_bam) { bam.prep <- unique.data.frame(bam.prep) } else NULL

bam.prep <- map.eval(bam.prep, bam.flags)
if (nrow(bam.prep) != nrow_bam) { bam.prep <- unique.data.frame(bam.prep) } else NULL

bam.prep$aln.cat <- 'primary_only'
bam.prep$aln.cat[bam.prep$has.secondary == T & bam.prep$has.supplementary == T ] <- 'sec_and_supp'
bam.prep$aln.cat[bam.prep$has.secondary == T & bam.prep$has.supplementary == F ] <- 'sec'
bam.prep$aln.cat[bam.prep$has.secondary == F & bam.prep$has.supplementary == F ] <- 'primary_only'
bam.prep$aln.cat[bam.prep$has.secondary == F & bam.prep$has.supplementary == T ] <- 'supp'





