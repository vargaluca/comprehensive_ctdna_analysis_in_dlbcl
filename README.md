# Phased variant and ctDNA fragment filtering tool
_Created as part of the the **Integrative genomic and fragmentomic analysis of cell-free DNA in DLBCL reveals locus- and patient-specific fragmentation patterns associated with progression-free survival** article._\
\
This tool is designed with minimal requirements in mind, so that analysis can be performed without the installation of specialised tools. The scripts utilise data from bam/sam files to determine the allelic orientation of variant pairs in close proximity, such as MD tags and sequence patterns to identify cis-occurring variant pairs. The tool gathers variants in <170bp distance - the approximate length of cfDNA fragments (call function), the determination of allelic orientation (filter function), and the creation of a bam file containing reads with phased variants (extract function).\
\
The tool works in a straightforward manner: a csv file containing already filtered variants can be used as an input, together with a bam file, from which variants were called. The scripts stitch together patterns to search from the bam file, to confirm the cis orientation of potential phased variants.
\
A snakemake wrapper script is also available, which makes the calling and performing all functions at once and calling multiple files easier.
## Software requirements:
- samtools (developed with: v1.23)
- Rscript (developed with: v4.4.1)
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

## Data processing
Example input file:

| chromosome | gene | entrez_id | position | reference_allele | tumor_allele | flanking_bps | t_depth | t_ref_count | t_alt_count |
|------------|------|-----------|----------|------------------|--------------|--------------|---------|-------------|-------------|
| chr2 | IGKJ4 | 0 | 89160428 | A | G | AAA | 105 | 90 | 15 |
| chr2 | IGKJ3 | 0 | 89160739 | T | G | ATA | 202 | 168 | 34 |
| chr2 | IGKJ2 | 0 | 89161049 | G | C | AGC | 252 | 235 | 17 |
| chr2 | IGKV1-8 | 0 | 89292393 | C | T | CCT | 48 | 44 | 4 |
| chr2 | IGKV1-9 | 0 | 89309565 | A | G | GAT | 57 | 55 | 2 |
| chr2 | IGKV3-15 | 0 | 89384698 | T | A | GTA | 15 | 8 | 7 |
| chr2 | IGKV3-15 | 0 | 89384705 | C | G | ACT | 18 | 10 | 8 |
| chr2 | IGKV3-15 | 0 | 89384715 | T | G | CTT | 24 | 12 | 12 |
| chr2 | IGKV3-15 | 0 | 89384761 | C | A | CCC | 128 | 66 | 62 |
| chr2 | IGKV3-15 | 0 | 89384831 | C | G | GCC | 53 | 30 | 23 |
| chr2 | IGKV3-15 | 0 | 89384840 | G | A | GGT | 43 | 23 | 20 |
| chr2 | IGKV1-17 | 0 | 89417262 | C | T | CCA | 109 | 105 | 4 |
| chr2 | IGKV2D-29 | 0 | 89986780 | G | A | AGA | 223 | 219 | 4 |
| chr2 | IGKV1D-43 | 0 | 90248944 | C | A | GCG | 92 | 88 | 4 |
| chr2 | MAP4K4 | 9448 | 102476297 | C | A | CCA | 85 | 81 | 4 |
| chr2 | MAP4K4 | 9448 | 102493511 | C | A | ACA | 142 | 138 | 4 |

For example, the first two IGKV1-15 variants in this table are potential phased variants, as they are <170 bp distance. The *call* process will piece together an unfiltered list of variant pairs that can be called from the respective bam file. The *filter* process then searches the patterns created by *call*: 
- a search pattern to find likely fragments (GAA....AGT)

In the case of the mentioned first two IGKV3-15 variants, the pattern will use flanking bases altered to include mutations and the number of wildcards between them.

- an MD-tag pattern (T6C)

This confirms that the fragments retrieved by the mutation patterns are results of mutation and not wild-type sequences captured by mistake 
Results of the capture can be output into a separate bam containing only ctDNA fragments and validated manually using IGV.

