# Phased variant and ctDNA fragment filtering tool
_Created as part of the the **Integrative genomic and fragmentomic analysis of cell-free DNA in DLBCL reveals locus- and patient-specific fragmentation patterns associated with progression-free survival** article._\
\
This tool is designed with minimal requirements in mind, so that analysis can be performed without the installation of specialised tools. The scripts utilise data from bam/sam files to determine the allelic orientation of variant pairs in close proximity, such as MD tags and sequence patterns to identify cis-occurring variant pairs. The tool gathers variants in <170bp distance - the approximate length of cfDNA fragments (call function), the determination of allelic orientation (filter function), and the creation of a bam file containing reads with phased variants (extract function).\
\
A snakemake wrapper script is also available, which makes the calling and performing all functions at once easier.
## Requirements:
- samtools (1.23)
- Rscript (4.4.1)
- R packages
