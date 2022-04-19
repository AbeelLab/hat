import copy
import concurrent.futures
from hat.ReadInput import Pos_to_Nuc_bam, ReadAlignmentsShortRead, ReadAlignmentsLongRead

def Connecting(phase_matrix, blocks, blocks_reads, k, var, max_count_parameter):
   if len(blocks) == 1:
       return phase_matrix, blocks
   variations_loc = sorted(list(var.keys()))
   new_phase_matrix = copy.deepcopy(phase_matrix)
   new_blocks = []
   new_blocks_reads = copy.deepcopy(blocks_reads)
   error_1 = 0
   error_2 = 0
   all_connections = []
   for i, block in enumerate(blocks[:-1]):
       connections = []
       for j in range(k):
           max_count = -1
           max = -1
           for l in range(k):
               count = 0
               for read_1 in blocks_reads[i][j]:
                   for read_2 in blocks_reads[i+1][l]:
                       if read_1[5].query_name == read_2[5].query_name:
                           count += 1
               if count > max_count:
                   max_count = count
                   max = l
           if max_count > max_count_parameter:
               connections.append((j ,max))
       if len(connections) < k:
           if len(connections) == k - 1 and len(set([x[1] for x in connections])) == k - 1:
               missing_source = (set(range(k)) ^ set([x[0] for x in connections])).pop()
               missing_dest = (set(range(k)) ^ set([x[1] for x in connections])).pop()
               connections.append((missing_source, missing_dest))
           error_1 += 1
       if len(set([x[1] for x in connections])) < k:
           error_2 += 1
       all_connections.append(connections)
   paths = []
   for u in range(k):
       paths.append([u])
   start = variations_loc.index(blocks[0][0])
   end = variations_loc.index(blocks[0][1])
   this_merge_blocks = [blocks[0]]
   for i, block in enumerate(blocks[1:]):
       con = all_connections[i]
       check = len(set([x[1] for x in con]))
       if len(con) == k and len(set([x[1] for x in con])) == k:
           this_merge_blocks.append(block)
           end = variations_loc.index(block[1])
           for p in range(len(paths)):
               for c in range(len(con)):
                   if paths[p][-1] == con[c][0]:
                       paths[p].append(con[c][1])
                       break
           if i == len(blocks) - 2:
               for haplo in range(k):
                   for p in range(len(paths[0])):
                       haplo_in_next_block = paths[haplo][p]
                       for p in range(variations_loc.index(this_merge_blocks[p][0]), variations_loc.index(this_merge_blocks[p][1]) + 1):  # check the +1
                           new_phase_matrix[haplo][p] = phase_matrix[haplo_in_next_block][p]
               new_blocks.append((variations_loc[start], variations_loc[end]))
       else:
           for haplo in range(k):
               for p in range(len(paths[0])):
                   haplo_in_next_block = paths[haplo][p]
                   for p in range(variations_loc.index(this_merge_blocks[p][0]), variations_loc.index(this_merge_blocks[p][1]) + 1): #check the +1
                       new_phase_matrix[haplo][p] = phase_matrix[haplo_in_next_block][p]
           new_blocks.append((variations_loc[start], variations_loc[end]))
           this_merge_blocks = [block]
           paths = []
           for u in range(k):
               paths.append([u])
           start = variations_loc.index(block[0])
           end = variations_loc.index(block[1])
   if not (variations_loc[start], variations_loc[end]) in new_blocks:
       new_blocks.append((variations_loc[start], variations_loc[end]))
       for haplo in range(k):
           for p in range(len(paths[0])):
               haplo_in_next_block = paths[haplo][p]
               for p in range(variations_loc.index(this_merge_blocks[p][0]),
                              variations_loc.index(this_merge_blocks[p][1]) + 1):  # check the +1
                   new_phase_matrix[haplo][p] = phase_matrix[haplo_in_next_block][p]
   return new_phase_matrix, new_blocks

def Fill_block(phase_matrix, vars, blocks, blocks_reads, short_read_alignments, iteration):
   var_locs = list(vars.keys())
   for b, block in enumerate(blocks):
       block_start = block[0]
       block_end = block[1]
       for pos in range(var_locs.index(block_start), var_locs.index(block_end)):
           var_pos = var_locs[pos]
           for h in range(len(phase_matrix)):
               if phase_matrix[h][pos] == 0:
                   allele = Majority_vote_Allreads(var_pos, blocks_reads[b][h], short_read_alignments, iteration)
                   if allele == None:
                       continue
                   else:
                       phase_matrix[h][pos] = allele
   return phase_matrix

def Majority_vote_Allreads(pos, reads, short_read_alignments, iteration):
   alleles = []
   if iteration == 0:
       for read in reads:
           read_name = read[5].qname
           pair_1 = short_read_alignments[read_name]['pair_1']
           pair_2 = short_read_alignments[read_name]['pair_2']
           alleles.append(Pos_to_Nuc_bam(pair_1, [pos]))
           alleles.append(Pos_to_Nuc_bam(pair_2, [pos]))
   else:
       for read in reads:
           alleles.append(Pos_to_Nuc_bam(read, [pos]))
   allele = Select_allele_majority_vote(alleles)
   return allele

def Select_allele_majority_vote(alleles):
   if "-" in alleles:
       alleles = remove_values_from_list(alleles, "-")
   if None in alleles:
       alleles = remove_values_from_list(alleles, None)
   alleles_unique = list(set(alleles))
   counts = []
   for allele in alleles_unique:
       counts.append(alleles.count(allele))
   if len(counts) == 0:
       return None
   copy_of_alleles = counts
   largest_integer = max(counts)  # 39
   if largest_integer < 2:
       print("Majority_voting_less_than_two")
       return None
   elif counts.count(largest_integer) > 1:
       print("Majority_voting_two_are_identical")
       return None
   else:
       return alleles_unique[counts.index(largest_integer)]

def remove_values_from_list(the_list, val):
  return [value for value in the_list if value != val]

def Fine_Tuning(phase_matrix, variations, blocks, blocks_reads, blocks_long_reads):
   variations_loc = sorted(list(variations.keys()))
   for b in range(len(blocks)):
       blocks_variations = variations_loc[variations_loc.index(blocks[b][0]): variations_loc.index(blocks[b][1]) + 1]
       for h in range(len(blocks_reads[b])):
           for var in blocks_variations:
               phased_allele = phase_matrix[h][variations_loc.index(var)]
               alleles_long = []
               alleles = []
               for read in blocks_reads[b][h]:
                   if not var in read[3]:
                       continue
                   pos_nucs = Pos_to_Nuc_bam(read, [var])
                   if pos_nucs == None:
                       continue
                   alleles.append(pos_nucs)
               if "-" in alleles:
                   alleles = remove_values_from_list(alleles, "-")
               if len(alleles) == 0:
                   continue
               for read in blocks_long_reads[b][h]:
                   if not var in read[3]:
                       continue
                   pos_nucs = Pos_to_Nuc_bam(read, [var])
                   if pos_nucs == None:
                       continue
                   alleles_long.append(pos_nucs)
               if "-" in alleles_long:
                   alleles_long = remove_values_from_list(alleles_long, "-")
               if len(alleles_long) == 0:
                   continue
               counts = []
               alleles_unique = list(set(alleles))
               for allele in alleles_unique:
                   counts.append(alleles.count(allele))
               largest_integer = max(counts)
               best_short = alleles_unique[counts.index(largest_integer)]
               alleles_long_unique = list(set(alleles_long))
               counts_long = []
               for allele in alleles_long_unique:
                   counts_long.append(alleles_long.count(allele))
               largest_integer_long = max(counts_long)
               best_long = alleles_long_unique[counts_long.index(largest_integer_long)]
               if phased_allele != best_short:
                   if phased_allele == 0:
                       phase_matrix[h][variations_loc.index(var)] = best_short
                   else:
                       print("fine tunning, short allele is different")
                       print(var, h,  phased_allele,best_long, best_short, alleles, alleles_long)
                       if best_short == best_long and largest_integer > 0.9*len(blocks_reads[b][h]) and largest_integer_long > 0.6* len(blocks_long_reads[b][h]):
                           phase_matrix[h][variations_loc.index(var)] = best_short
               elif phased_allele != best_long:
                   print(var, phased_allele, best_long, best_short, alleles, alleles_long)
                   print("fine tunning, long allele is different")
   print("under development")
   return phase_matrix

def Get_Block_ShortReads(block, variations_loc, short_read_alignment_file, chromosome_name):
   alignments = ReadAlignmentsShortRead(short_read_alignment_file,variations_loc, (block[0], block[1]), chromosome_name)
   return alignments

def Get_Block_LongReads(block, variations_loc, long_read_alignment_file, chromosome_name):
   alignments = ReadAlignmentsLongRead(long_read_alignment_file, variations_loc, (block[0], block[1]), chromosome_name)
   return alignments

def Assign_ShortReads_to_Blocks(phase_matrix, variations, blocks, k, short_read_alignment_file, chromosome_name):
   variations_loc = sorted(list(variations.keys()))
   blocks_reads = [[] for x in range(len(blocks))]
   for b in range(len(blocks)):
       alignments = Get_Block_ShortReads(blocks[b], variations_loc, short_read_alignment_file, chromosome_name)
       for read in alignments:
           blocks_variations = variations_loc[variations_loc.index(blocks[b][0]): variations_loc.index(blocks[b][1]) + 1]
           for pair in alignments[read].values():
               if len(pair) < 4:
                   continue
               shared_pos = set(blocks_variations).intersection(set(pair[3]))
               if Produce_k_variations(phase_matrix, variations, shared_pos, k):
                   blocks_reads[b].append(pair)
                   break
   blocks_reads_per_haplo = Assign_Reads_to_Haplotypes(phase_matrix, blocks_reads, variations, blocks, k)
   return blocks_reads_per_haplo

def Assign_Reads_to_Haplotypes(phase_matrix, blocks_reads, variations, blocks, k):
   variations_loc = sorted(list(variations.keys()))
   blocks_reads_per_haplo = [[] for x in range(len(blocks))]
   for b in range(len(blocks)):
       blocks_variations = variations_loc[variations_loc.index(blocks[b][0]): variations_loc.index(blocks[b][1]) + 1]
       blocks_reads_per_haplo[b] = [[] for x in range(k)]
       for pair in blocks_reads[b]:
           correct_haplo = []
           sequences = []
           for haplo in range(k):
               sequence = ""
               phased_locs = []
               test_vars = sorted(list(set(blocks_variations).intersection(set(pair[3]))))
               for var_loc in sorted(list(set(blocks_variations).intersection(set(pair[3])))):
                   if phase_matrix[haplo][variations_loc.index(var_loc)] != 0:
                       sequence = sequence + "," + phase_matrix[haplo][variations_loc.index(var_loc)]
                       phased_locs.append(var_loc)
               sequences.append(sequence)
               check = Pos_to_Nuc_bam(pair, phased_locs)
               if check == None:
                   continue
               if check != None:
                   letssee= check
               if sequence[1:] == check:
                   correct_haplo.append(haplo)
           if len(correct_haplo) == 1:
               blocks_reads_per_haplo[b][correct_haplo[0]].append(pair)
   return blocks_reads_per_haplo

def Phase_from_Seeds(var, seeds, k):
    phase_matrix = [[0] * len(var.keys()) for i in range(k)]
    blocks = []
    blocks_reads = []
    print("Phasing")
    seeds_sorted_name = list(seeds.keys())
    all_var_nums = list(var.keys())
    for i in range(len(seeds_sorted_name)):
       for j in range(len(seeds_sorted_name) - i - 1):
           first_num = int(seeds_sorted_name[j].split(",")[0])
           second_num = int(seeds_sorted_name[j + 1].split(",")[0])
           if first_num > second_num:
               seeds_sorted_name[j], seeds_sorted_name[j + 1] = seeds_sorted_name[j + 1], seeds_sorted_name[j]
    for l, seed in enumerate(seeds_sorted_name):
       var_nums = []
       for j, var_num in enumerate(seed.split(",")):
           var_nums.append(int(var_num))
       if phase_matrix[0][all_var_nums.index(var_nums[0])] != 0:
           continue
       count_reads = []
       for variation in seeds[seed]:
           count_reads.append(len(seeds[seed][variation]))
       top_ones = sorted(range(len(count_reads)), key=lambda i: count_reads[i], reverse=True)[-k:]
       read_names = []
       for j in range(len(var_nums)):
           for i in range(k):
               temp = list(seeds[seed].keys())[top_ones[i]].split(",")
               temp = temp[j]
               phase_matrix[i][all_var_nums.index(var_nums[j])] = temp
       for i in range(k):
           read_names.append(seeds[seed][list(seeds[seed].keys())[top_ones[i]]])
       if l == 0 and len(seeds) > 1:
           first_var_next_seed_idx = all_var_nums.index(int(seeds_sorted_name[l + 1].split(",")[0]))
           blocks.append((all_var_nums[0], all_var_nums[first_var_next_seed_idx - 1]))
       elif l == 0 and len(seeds) == 1:
           blocks.append((all_var_nums[0], all_var_nums[- 1]))
       elif l < len(seeds) - 1:
           first_var_next_seed_idx = all_var_nums.index(int(seeds_sorted_name[l + 1].split(",")[0]))
           blocks.append((int(seed.split(",")[0]), all_var_nums[first_var_next_seed_idx - 1]))
       else:
           blocks.append((int(seed.split(",")[0]), all_var_nums[-1]))

       blocks_reads.append(read_names)
    return phase_matrix, blocks, blocks_reads

def Assign_LongReads_to_Blocks_Simmilarity(phase_matrix, variations, blocks, k, long_read_alignment_file, chromosome_name):
   variations_loc = sorted(list(variations.keys()))
   list_of_multiprocess_arguments = [(b, variations_loc, phase_matrix, k, variations, long_read_alignment_file, chromosome_name) for b in blocks]
   with concurrent.futures.ThreadPoolExecutor(10) as executor:
       results = executor.map(Block_similarity_multithread, list_of_multiprocess_arguments)
   blocks_reads_per_haplo = list(results)
   return blocks_reads_per_haplo

def Block_similarity_multithread(arguments):
   b, variations_loc, phase_matrix, k, variations, alignment_file, chromosome_name = arguments
   blocks_variations = variations_loc[variations_loc.index(b[0]): variations_loc.index(b[1]) + 1]
   alignments = Get_Block_LongReads(b, blocks_variations, alignment_file, chromosome_name)
   list_of_multiprocess_arguments = [(alignments[read], variations_loc, phase_matrix, k, blocks_variations, variations) for read in alignments]
   with concurrent.futures.ThreadPoolExecutor(100) as executor:
       results = executor.map(Read_assignment_haploblock_multithread, list_of_multiprocess_arguments)
   read_assignment_in_block = list(results)
   blocks_reads_per_haplo = [[] for x in range(k)]
   for i in range(len(read_assignment_in_block)):
       if read_assignment_in_block[i] != -1:
           blocks_reads_per_haplo[read_assignment_in_block[i]].append(alignments[list(alignments.keys())[i]])
   return blocks_reads_per_haplo

def Produce_k_variations(phase_matrix, variations, shared_var_pos, k):
   var_locs = list(variations.keys())
   haplo_sequences = [""]*k
   for pos in shared_var_pos:
       sequences = []
       for h in range(k):
           sequences.append(phase_matrix[h][var_locs.index(pos)])
       if not 0 in sequences:
           for h in range(k):
               haplo_sequences[h] = haplo_sequences[h] + phase_matrix[h][var_locs.index(pos)]
   if len(set(haplo_sequences)) >= k:
       return True
   else:
       return False

def Distance_Sequences(read_sequence, phase_sequences):
   distances = []
   for seq in phase_sequences:
       count = 0
       for i, char in enumerate(seq):
           if char != read_sequence[i]:
               count += 1
       distances.append(count)
   return distances

def Read_assignment_haploblock_multithread(arguments):
   read, variations_loc, phase_matrix, k, blocks_variations, variations = arguments
   if len(read) < 4:
       return -1
   shared_var_pos = set(blocks_variations).intersection(set(read[3]))
   if not Produce_k_variations(phase_matrix, variations, shared_var_pos, k):
       return -1
   sequences = []
   phased_locs = []
   for var_loc in sorted(list(set(blocks_variations).intersection(set(read[3])))):
       alleles = []
       for i in range(k):
           alleles.append(phase_matrix[i][variations_loc.index(var_loc)])
       if not 0 in alleles:
           phased_locs.append(var_loc)
   for haplo in range(k):
       sequence = ""
       for var_loc in phased_locs:
           sequence = sequence + "," + phase_matrix[haplo][variations_loc.index(var_loc)]
       sequences.append(sequence[1:].split(","))
   read_sequence = Pos_to_Nuc_bam(read, phased_locs)
   if read_sequence == None:
       return -1
   read_sequence = read_sequence.split(",")
   distances = Distance_Sequences(read_sequence, sequences)
   min_distance = min(distances)
   correct_haplo = distances.index(min_distance)
   if distances.count(min_distance) > 1:
       return -1
   if min_distance <= len(phased_locs)* 0.25:
       return correct_haplo
   else:
       return -1
