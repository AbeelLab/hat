# HAT: Haplotype Assembly Tool 

HAT is a haplotype assembly tool that use NGS and TGS data along a reference genome to reconstruct haplotypes. 
HAT gets a VCF file containing SNPs, sorted bam file of NGS data alignment to the reference, and sorted bam file of TGS data alignment to the reference as input.


# Installation

    python setup.py install

## Requirements

- Python
- pysam
- Biopython
- numpy
- matplotlib
- seaborn

# Options

```
usage: HAT [-h] [-rl READ_LENGTH] [-pl PHASING_LOCATION] [-r REFERENCE_FILE] [-lf LONGREADS_FASTA]
           [-th TRUE_HAPLOTYPES] [-ma MULTIPLE_GENOME_ALIGNMENT] [-ha HAPLOTYPE_ASSEMBLY]
           chromosome_name vcf_file short_read_alignment long_read_alignment ploidy output output_dir

positional arguments:
  chromosome_name       The chromosome which is getting phased
  vcf_file              VCF file name
  short_read_alignment  short reads alignment file
  long_read_alignment   long reads alignment file
  ploidy                ploidy of the chromosome
  output                output prefix file name
  output_dir            output directory

optional arguments:
  -h, --help            show this help message and exit
  -rl READ_LENGTH, --read_length READ_LENGTH
                        short reads length
  -pl PHASING_LOCATION, --phasing_location PHASING_LOCATION
                        the location in the chromosome which is phased
  -r REFERENCE_FILE, --reference_file REFERENCE_FILE
                        reference file
  -lf LONGREADS_FASTA, --longreads_fasta LONGREADS_FASTA
                        long reads fasta file
  -th TRUE_HAPLOTYPES, --true_haplotypes TRUE_HAPLOTYPES
                        the correct haplotypes file
  -ma MULTIPLE_GENOME_ALIGNMENT, --multiple_genome_alignment MULTIPLE_GENOME_ALIGNMENT
                        Multiple genome alignment file of haplotypes to the reference
  -ha HAPLOTYPE_ASSEMBLY, --haplotype_assembly HAPLOTYPE_ASSEMBLY
                        Assembly of the haplotype sequences
```
# Example



