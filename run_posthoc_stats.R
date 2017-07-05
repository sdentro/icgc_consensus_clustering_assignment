#' Produce some post-hoc stats

q = 0.05
outdir = "output_wm"

source("~/repo/moritz_mut_assignment/MutationTime.R")
source("~/repo/moritz_mut_assignment/util.R")
source("~/repo/dpclust3p/R/interconvertMutationBurdens.R")

all_stats = data.frame()
for (infile in list.files(outdir, pattern="_assignment.RData", full.names=T)) {
  load(infile)
  qq_snv <- mean(MCN$D$pMutCNTail < q/2 | MCN$D$pMutCNTail > 1-q/2, na.rm=T)
  qq_indel <- mean(MCN_indel$D$pMutCNTail < q/2 | MCN_indel$D$pMutCNTail > 1-q/2, na.rm=T)
  qq_sv <- mean(MCN_sv$D$pMutCNTail < q/2 | MCN_sv$D$pMutCNTail > 1-q/2, na.rm=T)
  all_stats = rbind(all_stats, data.frame(samplename, qq_snv=qq_snv, qq_indel=qq_indel, qq_sv=qq_sv))
}

sel = all_stats[all_stats$qq_snv > 0.05, ]
sel = sel[with(sel, order(qq_snv, decreasing=T)),]
write.table(data.frame(samplename=sel$samplename), file="wm_stats_tail_snv.lst", col.names=F, sep="\t", row.names=F, quote=F)

sel = all_stats$qq_indel > 0.05
sel[is.na(sel)] = F
sel = all_stats[sel, ]
sel = sel[with(sel, order(qq_indel, decreasing=T)),]
write.table(data.frame(samplename=sel$samplename), file="wm_stats_tail_indel.lst", col.names=F, sep="\t", row.names=F, quote=F)
write.table(data.frame(samplename=all_stats$samplename[all_stats$qq_sv > 0.05]), file="wm_stats_tail_snv.lst", col.names=F, sep="\t", row.names=F, quote=F)

write.table(all_stats, file="wm_stats.txt", sep="\t", row.names=F, quote=F)

# for low numbers: Pbinom(n_out, n_total, 0.05, lower.tail=TRUE)



# 
# # Temp figure to see where badly explained SNVs are
# for (samplename in readr::read_tsv("wm_stats_tail_snv.lst", col_names=F)[[1]]) {
#   load(paste0("output_wm_old/", samplename, "_assignment.RData"))
#   temp_plot = all_data[all_data$type=="SNV",]
#   temp_plot$passed = ifelse(MCN$D$pMutCNTail < q/2 | MCN$D$pMutCNTail > 1-q/2, "FAIL", "PASS")
#   p = base_plot(temp_plot, "ccf", fill="passed")
#   png(file.path("wm_badly_explained_data", paste0(samplename, "_snv.png")), width=600, height=300)
#   print(p)
#   dev.off()
# }
# 
# for (samplename in head(readr::read_tsv("wm_stats_tail_indel.lst", col_names=F)[[1]], 50)) {
#   load(paste0("output_wm_old/", samplename, "_assignment.RData"))
#   temp_plot = all_data[all_data$type=="indel",,drop=F]
#   temp_plot$passed = ifelse(MCN_indel$D$pMutCNTail < q/2 | MCN_indel$D$pMutCNTail > 1-q/2, "FAIL", "PASS")
#   p = base_plot(temp_plot, "ccf", fill="passed")
#   png(file.path("wm_badly_explained_data", paste0(samplename, "_indel.png")), width=600, height=300)
#   print(p)
#   dev.off()
# }
# 
