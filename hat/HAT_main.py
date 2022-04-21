import os
from hat.ReadInput import ParseVCF_pysam, Genome_size_from_reference, ReadAlignmentsShortRead, ReadAlignmentsLongRead, Average_Read_Length
from hat.Seeds import Creating_Seeds, Filter_variations, Filter_seeds, Remove_overlaping_seeds
from hat.PloidyBlocks import Finding_ploidy_blocks_alg3, Ploidy_blocks_coverage, Get_Ploidy_Block_seeds
from hat.Visualization import Draw_seeds_on_chromosome
from hat.Phasing import Phase_from_Seeds, Assign_ShortReads_to_Blocks, Assign_LongReads_to_Blocks_Simmilarity, Connecting, Fill_block
import argparse
import csv
import numpy as np

def print_number_phased_variants(phase_matrix, blocks):
   count = 0
   for j in range(len(phase_matrix[0])):
       check = True
       for i in range(len(phase_matrix)):
           if phase_matrix[i][j] == 0:
               check = False
       if check:
           count += 1
   print("blocks = ", len(blocks), " phased_variants = ", count)
   return count, len(blocks)

def Run_HAT(vcf_file, read_length, gene_location, chromosome_name, reference_fa, short_read_alignment_file, long_read_alignment_file, long_reads_fa, assembly, ploidy, output_prefix, output_dir, longreads_fasta, shortreads_1_fastq, shortreads_2_fastq):
    variations = ParseVCF_pysam(vcf_file, gene_location, chromosome_name)
    var_locs = list(variations.keys())
    genome_size = Genome_size_from_reference(reference_fa)
    print("variations", len(variations.keys()))
    alignments = ReadAlignmentsShortRead(short_read_alignment_file, variations, gene_location, chromosome_name)
    print(len(alignments.keys()), "number of short reads")
    seeds, reads_to_seeds = Creating_Seeds(alignments)
    print(len(seeds.keys()), seeds.keys())
    seeds = Filter_variations(seeds, ploidy)
    seeds = Filter_seeds(seeds, 2)
    print(len(seeds.keys()), "removedlowsupport")
    seeds = Remove_overlaping_seeds(seeds)
    print(len(seeds.keys()), "removedoverlapping")
    if long_reads_fa != None:
        average_read_length = Average_Read_Length(long_reads_fa)
    else:
        average_read_length = 6300
    ploidy_blocks = Finding_ploidy_blocks_alg3(seeds, average_read_length)
    Ploidy_blocks_coverage(ploidy_blocks, var_locs, genome_size, ploidy)
    Draw_seeds_on_chromosome(var_locs, seeds, ploidy_blocks, genome_size, output_dir, output_prefix)
    all_blocks = []
    all_short_reads = []
    all_long_reads = []
    all_phase_matrix = []
    print("ploidy blocks", ploidy_blocks)
    for ploidy_b in ploidy_blocks.keys():
        print(ploidy_b)
        k = ploidy_blocks[ploidy_b]
        start_ploidy_b = int(ploidy_b.split(",")[0])
        end_ploidy_b = int(ploidy_b.split(",")[-1])
        b_variations = ParseVCF_pysam(vcf_file, [start_ploidy_b - 1, end_ploidy_b + 1], chromosome_name)
        b_alignments = ReadAlignmentsShortRead(short_read_alignment_file, b_variations,
                                           [start_ploidy_b - 1, end_ploidy_b + 1], chromosome_name)
        b_long_read_alignments = ReadAlignmentsLongRead(long_read_alignment_file, b_variations,
                                                             [start_ploidy_b - 1, end_ploidy_b + 1],
                                                             chromosome_name,
                                                             )
        b_seeds = Get_Ploidy_Block_seeds(start_ploidy_b, end_ploidy_b, seeds)
        b_seeds = Filter_seeds(b_seeds, k)
        phase_matrix, blocks, blocks_reads = Phase_from_Seeds(b_variations, b_seeds, k)
        print("Initialization is done")
        old_phase_matrix_count, old_block_count = print_number_phased_variants(phase_matrix, blocks)
        we_have_change = True
        iteration = 0
        while we_have_change or iteration < 3:
            print("ITERATION--------: ", iteration)
            if iteration > 1:
                iteration_max_count_parameter = 3
            else:
                iteration_max_count_parameter = 1
            if iteration == 0:
                blocks_reads = Assign_ShortReads_to_Blocks(phase_matrix, b_variations, blocks, k, short_read_alignment_file, chromosome_name)
                phase_matrix, blocks = Connecting(phase_matrix, blocks, blocks_reads, k, b_variations,
                                                  max_count_parameter=iteration_max_count_parameter)
                blocks_reads = Assign_ShortReads_to_Blocks(phase_matrix, b_variations, blocks, k, short_read_alignment_file, chromosome_name)
                phase_matrix = Fill_block(phase_matrix, b_variations, blocks, blocks_reads, b_alignments, iteration)
                print("Iteration 0 is overّ")
            else:
                blocks_long_reads = Assign_LongReads_to_Blocks_Simmilarity(phase_matrix, b_variations, blocks, k, long_read_alignment_file, chromosome_name)
                phase_matrix, blocks = Connecting(phase_matrix, blocks, blocks_long_reads, k, b_variations,
                                                  max_count_parameter=iteration_max_count_parameter)

                blocks_long_reads = Assign_LongReads_to_Blocks_Simmilarity(phase_matrix, b_variations, blocks, k,
                                                                           long_read_alignment_file, chromosome_name)
                phase_matrix = Fill_block(phase_matrix, b_variations, blocks, blocks_long_reads, b_alignments,
                                          iteration)
                print("iteration " + str(iteration) + " is overّ")

            new_phase_matrix_count, new_block_count = print_number_phased_variants(phase_matrix, blocks)
            if new_phase_matrix_count != old_phase_matrix_count or new_block_count != old_block_count:
                we_have_change = True
            else:
                we_have_change = False
            iteration += 1
            old_block_count = new_block_count
            old_phase_matrix_count = new_phase_matrix_count

        print("iterative done, next block?")
        blocks_reads = Assign_ShortReads_to_Blocks(phase_matrix, b_variations, blocks, k, short_read_alignment_file, chromosome_name)
        blocks_long_reads = Assign_LongReads_to_Blocks_Simmilarity(phase_matrix, b_variations, blocks, k,
                                                                   long_read_alignment_file, chromosome_name)
        all_phase_matrix.append(phase_matrix)
        all_short_reads.append(blocks_reads)
        all_long_reads.append(blocks_long_reads)
        all_blocks.append(blocks)
    print("Phasing finished")
    Save_output(all_phase_matrix, ploidy_blocks, all_blocks, all_short_reads, all_long_reads, assembly, output_prefix, output_dir, longreads_fasta, shortreads_1_fastq, shortreads_2_fastq)

def Save_output(all_phase_matrix, ploidy_blocks, all_blocks, all_short_reads, all_long_reads, assembly, output_prefix, output_dir, longreads_fasta, shortreads_1_fastq, shortreads_2_fastq):
    with open(output_dir + output_prefix + "_phase_matrix.csv", "w") as f:
        for b, ploidy_b in enumerate(list(ploidy_blocks.keys())):
            f.write(" PLOIDY BLOCK    ")
            f.write(str(ploidy_b))
            f.write("\n")
            f.write(str(all_blocks[b]))
            f.write("\n")
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(all_phase_matrix[b])
    with open(output_dir + output_prefix + "_phased_reads_blocks", "w") as f:
        for b, ploidy_b in enumerate(list(ploidy_blocks.keys())):
            f.write(" PLOIDY BLOCK     ")
            f.write(str(ploidy_b))
            f.write("\n")
            for b_h in range(len(all_short_reads[b])):
                f.write("----------_BLOCK-----------\n")
                for h in range(len(all_short_reads[b][b_h])):
                    f.write("h " + str(h) + "\n")
                    f.write("--------------short reads-----------\n")
                    for read in all_short_reads[b][b_h][h]:
                        if len(read) < 6:
                            continue
                        f.write(read[5].qname + "\n")
                    f.write("\n")
                    f.write("--------------long reads-----------\n")
                    for read in all_long_reads[b][b_h][h]:
                        if len(read) < 6:
                            continue
                        f.write(read[5].qname + "\n")
                    f.write("\n")


    zeros = 0
    size = 0
    for pb in range(len(all_phase_matrix)):
        for b in range(len(all_phase_matrix[pb])):
            for h in range(len(all_phase_matrix[pb][b])):
                if all_phase_matrix[pb][b][h] == 0:
                    zeros += 1
                size += 1

    print("phase_matrix zeros: ", zeros)
    print("phase_matrix total elements: ", size)
    if assembly:
        for pb, pb_name in enumerate(ploidy_blocks.keys()):
            for b in range(len(all_short_reads[pb])):
                for h in range(len(all_short_reads[pb][b])):
                    haplotype_dir = str(pb_name.replace(",", '-')) + "_" + str(all_blocks[pb][b][0]) + "-" + str(
                        all_blocks[pb][b][1]) + "/"
                    if not os.path.exists(output_dir + haplotype_dir):
                        os.mkdir(output_dir + haplotype_dir)
                    print(output_dir + haplotype_dir)
                    with open(output_dir + haplotype_dir + output_prefix + "_" + str(h) + '_short_read.list',
                              "w") as out_short_reads_haplotype_block:
                        for record in all_short_reads[pb][b][h]:
                            out_short_reads_haplotype_block.write(record[5].qname + "\n")
                    os.system("seqkit grep -f " + output_dir + haplotype_dir + output_prefix + "_" + str(
                        h) + '_short_read.list ' + shortreads_2_fastq + " > " + output_dir + haplotype_dir + output_prefix + "_" + str(
                        h) + '_short_read_1.fastq')
                    if short_2_path != None:
                        os.system("seqkit grep -f " + output_dir + haplotype_dir + output_prefix + "_" + str(
                            h) + '_short_read.list ' + shortreads_1_fastq + " > " + output_dir + haplotype_dir + output_prefix + "_" + str(
                            h) + '_short_read_2.fastq')
                    with open(output_dir + haplotype_dir + output_prefix + "_" + str(h) + '_long_read.list',
                              "w") as out_long_reads_haplotype_block:
                        for record in all_long_reads[pb][b][h]:
                            out_long_reads_haplotype_block.write(record[5].qname + "\n")
                    os.system("seqkit grep -f " + output_dir + haplotype_dir + output_prefix + "_" + str(
                        h) + '_long_read.list ' + longreads_fasta + " > " + output_dir + haplotype_dir + output_prefix + "_" + str(
                        h) + '_long_read.fasta')
                    pb_size = int(pb_name.split(",")[1]) + 12600 - int(pb_name.split(",")[0]) + 1
                    print("pb name", pb_name, pb_size)
                    long_reads_path = output_dir + haplotype_dir + output_prefix + "_" + str(h) + '_long_read.fasta'
                    short_1_path = output_dir + haplotype_dir + output_prefix + "_" + str(h) + '_short_read_1.fastq'
                    short_2_path = output_dir + haplotype_dir + output_prefix + "_" + str(h) + '_short_read_2.fastq'
                    long_reads_overlap = output_dir + haplotype_dir + str(h) + "_overlaps.paf"
                    os.system(
                        "minimap2 -x ava-ont " + long_reads_path + " " + long_reads_path + " > " + long_reads_overlap)
                    os.system(
                        "miniasm -f " + long_reads_path + " " + long_reads_overlap + " " + output_dir + haplotype_dir + str(
                            h) + "_assembly.gfa")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-rl", "--read_length", help="short reads length", type=int, default=125)
    parser.add_argument("-pl", "--phasing_location", help="the location in the chromosome which is phased", type=str)
    parser.add_argument("-r", "--reference_file", help="reference file", type=str)
    parser.add_argument("-lf", "--longreads_fasta", help="long reads fasta file", type=str)
    parser.add_argument("-sf1", "--shortreads_1_fastq", help="first pair fastq file")
    parser.add_argument("-sf2", "--shortreads_2_fastq", help="second pair fastq file")
    parser.add_argument("-th", "--true_haplotypes", help="the correct haplotypes file", type=str)
    parser.add_argument("-ma", "--multiple_genome_alignment",
                        help="Multiple genome alignment file of haplotypes to the reference", type=str)
    parser.add_argument("-ha", "--haplotype_assembly", help="Assembly of the haplotype sequences", type=bool,
                        default=False)

    parser.add_argument("chromosome_name", help="The chromosome which is getting phased", type=str)
    parser.add_argument("vcf_file", help="VCF file name", type=str)
    parser.add_argument("short_read_alignment", help="short reads alignment file", type=str)
    parser.add_argument("long_read_alignment", help="long reads alignment file", type=str)
    parser.add_argument("ploidy", help="ploidy of the chromosome", type=int)
    parser.add_argument("output", help="output prefix file name", type=str)
    parser.add_argument("output_dir", help="output directory", type=str)
    args = parser.parse_args()
    read_length = args.read_length
    if args.phasing_location != None:
        gene_location = args.phasing_location.split(",")
        gene_location[0] = int(gene_location[0])
        gene_location[1] = int(gene_location[1])
    else:
        gene_location = args.phasing_location
    chromosome_name = args.chromosome_name
    long_read_alignment_file = args.long_read_alignment
    short_read_alignment_file = args.short_read_alignment
    vcf_file = args.vcf_file
    true_haplotypes_file = args.true_haplotypes
    short_reads_1 = args.shortreads_1_fastq
    short_reads_2 = args.shortreads_2_fastq
    reference_fa = args.reference_file
    long_reads_fa = args.longreads_fasta
    ploidy = args.ploidy
    output_dir = args.output_dir
    reference_alignment_file = args.multiple_genome_alignment
    assembly = args.haplotype_assembly
    if assembly:
        if short_reads_1 == None or long_reads_fa == None:
            print("The long reads fasta file and short reads fastq file is required for assembling the haplotypes")
            return
    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    output_prefix = args.output
    if short_read_alignment_file != None and long_read_alignment_file != 'None':
        print("Both short reads and long reads are available")
        results = Run_HAT(vcf_file, read_length, gene_location, chromosome_name, reference_fa, short_read_alignment_file, long_read_alignment_file, long_reads_fa, assembly, ploidy, output_prefix, output_dir, long_reads_fa, short_reads_1, short_reads_2)
    elif long_read_alignment_file == 'None' and short_read_alignment_file != None:
        print("short reads only, under development")
    elif short_read_alignment_file == 'None' and long_read_alignment_file != None:
        print("long reads only, under development")
    else:
        print("No alignment file")
