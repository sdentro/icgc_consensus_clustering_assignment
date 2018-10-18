# ICGC PCAWG-11 consensus mutation assignment

## Overview

This repo contains code used to assign mutations (SNV, indel and SV) to consensus mutation clusters from subclonal reconstruction methods.

 
## Dependencies

Software packages used to develop the code and run the pipeline on the PCAWG dataset. Installation of these packages should normally take a few minutes via Bioconductor.

```
R (version 3.1.0)
```

Internally, the pipeline calls [MutationTimeR](https://github.com/gerstung-lab/MutationTimeR)
```
MutationTimeR
```

R libraries (all installed via Bioconductor)
```
Bioconductor (version 3.0)
BiocInstaller (version 1.16.5)
ggplot2
gridExtra
grid
```

## How to run

```
R --no-save --no-restore --vanilla -f run.R --args \
-l [path to where the repository is downloaded] \
--sam [samplename] \
-o [output directory] \
--snv [PCAWG consensus SNV VCF file] \
--cna [PCAWG consensus copy number profile] \
--struct [sample subclonal architecture] \
--pur [tumour purity] \
--summ ${summ_tab} \
--ind [PCAWG consensus SNV indel file] \
--sv [PCAWG consensus SV VCF file] \
--sv_vaf [output file from SVclone with VAF values and copy number mapping for each SV]
```

## Output - run.R

This script runs the assignments and produces the following output files per sample

|Filename | Description |
|---|---|
|[samplename]_assignment.RData | All output stored in an RData archive (requires loading of dependencies) |
|[samplename]_pcawg11_output.RData | All PCAWG-11 output (no dependencies required) |
|[samplename]_summary_table_entry.txt | A PCAWG-11 summary table entry for this sample |
|[samplename]_final_assignment.png | A figure showing the data and assignment, with comparison to regular binomial assignment |

## Output - run_final_output.R

This script produces the output that is shared PCAWG-wide

 
## Output - run_final_mutccf.R

This script produces per mutation raw CCF estimates

## Output - run_pcawg11.R

Produces PCAWG-11 formatted output for internal use

 * [samplename]_subclonal_structure.txt
 * [samplename]_mutation_assignments.txt
 * [samplename]_multiplicity.txt