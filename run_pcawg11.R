# Dump output in PCAWG11 format (i.e. the pcawg11_output.RData file from run.R)
# later runs do not need the samplename parsing as the variable is already included in the archive

args = commandArgs(T)
infile = args[1]
outdir = args[2]
print(infile)
load(infile)
samplename = paste(unlist(stringr::str_split(basename(infile), "_"))[1], collapse="_")
#samplename = unlist(stringr::str_split(basename(infile), "_"))[1]
write.table(final_pcawg11_output$final_clusters, file=paste0(outdir, samplename, "_subclonal_structure.txt"), row.names=F, sep="\t", quote=F)
write.table(final_pcawg11_output$snv_assignments, file=paste0(outdir, samplename,"_mutation_assignments.txt"), row.names=F, sep="\t", quote=F)
write.table(final_pcawg11_output$snv_mult, file=paste0(outdir, samplename,"_multiplicity.txt"), row.names=F, sep="\t", quote=F)
