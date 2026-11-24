# Phased variant and ctDNA fragment filtering tool
_Created as part of the the **Integrative genomic and fragmentomic analysis of cell-free DNA in DLBCL reveals locus- and patient-specific fragmentation patterns associated with progression-free survival** article._\
\
This tool is designed with minimal requirements in mind, so that analysis can be performed without the installation of specialised tools. The scripts utilise data from bam/sam files to determine the allelic orientation of variant pairs in close proximity, such as MD tags and sequence patterns to identify cis-occurring variant pairs. The tool gathers variants in <170bp distance - the approximate length of cfDNA fragments (call function), the determination of allelic orientation (filter function), and the creation of a bam file containing reads with phased variants (extract function).\
\
A snakemake wrapper script is also available, which makes the calling and performing all functions at once and calling multiple files easier.
## Requirements:
- samtools (1.23)
- Rscript (4.4.1)
- R packages: dplyr, readr

Options:
```
    -h, --help                show brief help
    -p, --process             'call', 'only-call', 'only-filter', 'extract', 'recall'
    -i, --input               specify input variants csv file
    -b, --bam                 specify corresponding input bam file
    -o, --outdir              specify output directory
    -t, --ct-variants         specify ctdna variants file
    -m, --mode                0, 1, 2 variant requirements for recall function
```
## Processes:
### call
Outputs two csv files:
- unfiltered: containing all potential PV pairs
- filtered: only containinc cis-allelic variant pairs confirmed from corresponding bam files
```
bash phased_variants.sh -p call -i /data/variants.csv -b /data/tissue.bam -o /data/outfolder
```
### recall
Recalls identified PV pairs from corresponding cfDNA bam file and creates an output csv file with variant pairs\
modes (default: 1):
- 0: doesn't require variants called from cfDNA
- 1: requires at least one variant called from cfDNA
- 2: requires both variants of a pair to be called from cfDNA
```
bash phased_variants.sh -p recall -i /data/filtered_pv.csv -b /data/cfdna.bam -t /data/cfdna_variants.csv -o /data/outfolder -m 1
```

### extract
Creates bam file only containing reads with PV pairs from the cfDNA sample
```
bash phased_variants.sh -p extract -i /data/variants.csv -b /data/variants.bam -o /data/outfolder
```

### only-call
Outputs a csv file with potential phased variant pairs\
(has not been checked for allelic position)

### only-filter
Outputs a filtered PV file checked for allelic positions\
(only works with one already existing unfiltered file)
