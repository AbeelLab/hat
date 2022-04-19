from Bio import SeqIO
import pysam


def ParseVCF_pysam(filepath, gene_location, chr_name):
   variations = {}
   vcf_file = pysam.VariantFile(filepath)
   if gene_location != None:
       for rec in vcf_file.fetch(chr_name, gene_location[0], gene_location[1]):
           if rec.rlen > 1:
               continue
           variations[rec.pos] = rec
   else:
       for rec in vcf_file.fetch(chr_name):
           if rec.rlen > 1:
               continue
           variations[rec.pos] = rec
   return variations

def Genome_size_from_reference(reference_fa):
   record = SeqIO.read(reference_fa, "fasta")
   return len(record.seq)

def Variation_Coverage_bam(variations, positions):
   cov = []
   if type(variations) == dict:
       vars = variations.keys()
   else:
       vars = variations
   for pos in positions:
       if pos in vars:
           cov.append(pos)
   return cov

def ReadAlignmentsShortRead(bam_file, variations, gene_location, chr_name):
   alignments = {}
   counter = 0
   bamfile = pysam.AlignmentFile(bam_file, "rb")
   if gene_location != None:
       iter = bamfile.fetch(chr_name, gene_location[0], gene_location[1])
   else:
       iter = bamfile.fetch(chr_name)
   for align in iter:
       read_name = align.qname
       if align.flag in [99, 147, 83, 163]:
           counter += 1
           if align.is_read1:
               pair = 'pair_1'
           else:
               pair = 'pair_2'
           if align.is_reverse:
               strand = "-"
           else:
               strand = "+"
           variant_coverage = Variation_Coverage_bam(variations, align.positions)
           if read_name in alignments.keys():
               alignments[read_name][pair] = [align.reference_start, align.reference_end, strand, variant_coverage, align.cigarstring, align]
           else:
               alignments[read_name] = {'pair_1': [], 'pair_2': []}
               alignments[read_name][pair] = [align.reference_start, align.reference_end, strand, variant_coverage,
                                              align.cigarstring, align]
   return alignments

def Average_Read_Length(long_reads_fa):
   tot = 0
   num = 0
   for record in SeqIO.parse(long_reads_fa, "fastq"):
       num += 1
       tot += len(record.seq)
   result = float(tot) / num
   return result

def ReadAlignmentsLongRead(cram_file, variations, gene_loc, chr_name):
   haplotype_of_long_reads = {}

   alignments = {}
   counter = 0
   cramfile = pysam.AlignmentFile(cram_file, "rb")
   if gene_loc != None:
       iter = cramfile.fetch(chr_name, gene_loc[0], gene_loc[1] )
   else:
       iter = cramfile.fetch(chr_name)
   for align in iter:
       if len(haplotype_of_long_reads.keys()) > 0:
           read_name = align.qname + "-" + str(haplotype_of_long_reads[align.qname])
           align.qname = align.qname + "-" + str(haplotype_of_long_reads[align.qname])
       else:
           read_name = align.qname
           align.qname = align.qname
       counter += 1
       if align.is_reverse:
           strand = "-"
       else:
           strand = "+"
       variant_coverage = Variation_Coverage_bam(variations, align.positions)
       if read_name in alignments.keys():
           alignments[read_name] = [align.reference_start, align.reference_end, strand, variant_coverage, align.cigarstring, align]
       else:
           alignments[read_name] = [align.reference_start, align.reference_end, strand, variant_coverage, align.cigarstring, align]
   return alignments

def Pos_to_Nuc_bam(pair, poses):
   nucs = []
   for pos in poses:
       pos = pos - 1
       if len(pair) < 6:
           return None
       for query_to_reference in pair[5].aligned_pairs:
           if query_to_reference[1] == pos:
               q_pos = query_to_reference[0]
               if q_pos == None or pair[5].seq == None or not pair[5].seq[q_pos].upper() in ["A", "C", "G", "T"]:
                   nucs.append("-")
               else:
                   nucs.append(pair[5].seq[q_pos].upper())
               break
   if len(nucs) != len(poses):
       return None
   return ",".join(nucs)

