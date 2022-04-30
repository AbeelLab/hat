# HAT: Haplotype Assembly Tool 

HAT is a haplotype assembly tool that use NGS and TGS data along a reference genome to reconstruct haplotypes. 
HAT gets a VCF file containing SNPs, sorted bam file of NGS data alignment to the reference, and sorted bam file of TGS data alignment to the reference as input.


# Installation
```
git clone https://github.com/AbeelLab/hat
cd hat
pip install .
```
or
```
pip install HAT-phasing
```
or
```
conda install -c bioconda hat-phasing
```

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
           [-sf1 SHORTREADS_1_FASTQ] [-sf2 SHORTREADS_2_FASTQ] [-th TRUE_HAPLOTYPES]
           [-ma MULTIPLE_GENOME_ALIGNMENT] [-ha HAPLOTYPE_ASSEMBLY]
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
  -sf1 SHORTREADS_1_FASTQ, --shortreads_1_fastq SHORTREADS_1_FASTQ
                        first pair fastq file
  -sf2 SHORTREADS_2_FASTQ, --shortreads_2_fastq SHORTREADS_2_FASTQ
                        second pair fastq file
  -th TRUE_HAPLOTYPES, --true_haplotypes TRUE_HAPLOTYPES
                        the correct haplotypes file
  -ma MULTIPLE_GENOME_ALIGNMENT, --multiple_genome_alignment MULTIPLE_GENOME_ALIGNMENT
                        Multiple genome alignment file of haplotypes to the reference
  -ha HAPLOTYPE_ASSEMBLY, --haplotype_assembly HAPLOTYPE_ASSEMBLY
                        Assembly of the haplotype sequences
```

The -ha option requires minimap2, bwa, miniasm, seqkit and Pilon.  Currently, HAT looks in the PATH to find these tools. In addition, -ha option needs long reads fasta file and short reads fastq file.

# Example

To reconstruct the haplotypes, HAT needs 3 input:

- Sorted bam file from the alignment of short reads to the reference genome
- Sorted bam file from the alignment of long reads to the reference genome
- SNPs selected from a VCF file created by any variant calling tool.

In this example, we have used bwa mem and minimap2 for the alignments and FreeBayes to find variants. Later, we have used vcffilter to select the SNPs from the VCF file. The input files are:

- Example/haplosim-triploid-CP048984.1-highhetero/short_reads_alignment.sorted.bam
- Example/haplosim-triploid-CP048984.1-highhetero/long_reads_alignment.sorted.bam
- Example/haplosim-triploid-CP048984.1-highhetero/snp-var.vcf.gz

To reconstruct the haplotypes, we use the following command:

```
cd Example/haplosim-triploid-CP048984.1-highhetero/
HAT -r CP048984.1.fna CP048984.1 snp-var.vcf.gz short_reads_alignment.sorted.bam \
 long_reads_alignment.sorted.bam 3 hat_output results/
```

This command phase the CP048984.1 chromosome, and provide 4 outputs in the results directory:

- hat_output_ploidy_blocks figure
- hat_output_phase_matrix
- hat_output_phased_blocks

The ploidy blocks are the regions that have sufficient differences between the haplotypes. HAT operates in these regions and find the alleles of the haplotypes. The following figure shows the ploidy blocks HAT found in the example dataset.

![Ploidy blocks HAT finds for the example dataset. The long, black vertical lines at the bottom show
the SNPs and their positions on the chromosome found by FreeBayes. From these SNPS, HAT finds the seeds shown in short, black vertical lines in panel above the SNPs. The seeds are
placed vertically based on the number of combination of alleles they have, ranging from 1 to 6 (y axis). HAT uses these seeds to find ploidy blocks, which are shaded regions encapsulating
the seeds and the color of the region indicates the estimated ploidy level. See the legend for the colors corresponding to different ploidy levels.](Example/haplosim-triploid-CP048984.1-highhetero/results/ploidy_blocks.png "")


 The phase_matrix output has the alleles haplotypes for the SNP loci.

The phased_blocks output has the read clustering and shows the reads that are assigned to haplotypes and blocks.
