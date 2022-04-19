from hat.ReadInput import Pos_to_Nuc_bam

def Seeds_by_ploidy(seeds):
   seeds_by_ploidy_dict = {}
   for seed in seeds.keys():
       ploidy = len(seeds[seed].keys())
       if ploidy in seeds_by_ploidy_dict.keys():
           seeds_by_ploidy_dict[ploidy].append(seed)
       else:
           seeds_by_ploidy_dict[ploidy] = [seed]
   return seeds_by_ploidy_dict

def Sort_seeds_ploidy_loc(seeds):
   sorted_seeds = list(seeds)
   sorted_ploidies = []
   for i in range(len(sorted_seeds)):
       min_idx = 0
       max_ploidy = 0
       min_loc = 10000000000000000000000
       for j in range(i, len(sorted_seeds)):
           if len(seeds[sorted_seeds[j]]) > max_ploidy or  ( len(seeds[sorted_seeds[j]]) == max_ploidy  and int(sorted_seeds[j].split(",")[0]) < min_loc ):
               min_loc = int(sorted_seeds[j].split(",")[0])
               max_ploidy = len(seeds[sorted_seeds[j]])
               min_idx = j
       temp = sorted_seeds[i]
       sorted_seeds[i] = sorted_seeds[min_idx]
       sorted_seeds[min_idx] = temp
       sorted_ploidies.append(len(seeds[sorted_seeds[i]]))
   print(sorted_ploidies)
   print(sorted_seeds)
   return sorted_seeds

def Remove_overlaping_seeds(seeds):
   seeds_sorted_name = list(seeds.keys())
   for i in range(len(seeds_sorted_name)):
        for j in range(len(seeds_sorted_name) - i - 1):
            first_num = int(seeds_sorted_name[j].split(",")[0])
            second_num = int(seeds_sorted_name[j + 1].split(",")[0])
            if first_num > second_num:
                seeds_sorted_name[j], seeds_sorted_name[j + 1] = seeds_sorted_name[j + 1], seeds_sorted_name[j]
   delete_cand = []
   for i, seed in enumerate(seeds_sorted_name):
       if seed in delete_cand:
           continue
       for seed_2 in seeds_sorted_name:
           if seed == seed_2:
               continue
           if seed in delete_cand or seed_2 in delete_cand:
               continue
           if int(seed_2.split(",")[0]) > int(seed.split(",")[-1]):
               break
           if  int(seed_2.split(",")[-1]) < int(seed.split(",")[0]):
               continue
           seed_numbs = []
           seed_str = seed.split(',')
           for ch in seed_str:
               seed_numbs.append(int(ch))
           seed_2_numbs = []
           seed_2_str = seed_2.split(',')
           for ch in seed_2_str:
               seed_2_numbs.append(int(ch))
           if seed[0] < seed_2[0] or (seed[0] == seed_2[0] and len(seed_numbs) < len(seed_2_numbs)):
               first_seed = seed
               second_seed = seed_2
               first_seed_numbs = seed_numbs
               second_seed_numbs = seed_2_numbs
           else:
               first_seed = seed_2
               second_seed = seed
               first_seed_numbs = seed_2_numbs
               second_seed_numbs = seed_numbs
           for numb in first_seed_numbs:
               if numb in second_seed_numbs:
                   support_1 = Read_support(seeds[first_seed])
                   support_2 = Read_support(seeds[second_seed])
                   if len(seeds[second_seed].keys()) > len(seeds[first_seed].keys()) and support_2 >= 0.7* support_1:
                        delete_cand.append(first_seed)
                   else:
                        delete_cand.append(second_seed)
   for cand in delete_cand:
       try:
           del seeds[cand]
       except:
           continue
   return seeds


def Read_support(seed):
   read_sup = 0
   for variation_comb in seed.keys():
       read_sup += len(seed[variation_comb])
   return read_sup


def Filter_seeds(seeds, k):
   del_cand = []
   for seed in seeds.keys():
       if len(seeds[seed].keys()) < k:
           del_cand.append(seed)
   unique_del_cand = set(del_cand)
   for cand in unique_del_cand:
       del seeds[cand]
   return seeds

def RemoveNestedSeeds(seeds):
   delete_cand = []
   for seed in seeds.keys():
       for seed_2 in seeds.keys():
           if seed in seed_2:
               if seed in delete_cand or seed_2 in delete_cand:
                   continue
               if seed == seed_2:
                   continue
               delete_cand.append(seed)
   delete_cand = set(delete_cand)
   for cand in delete_cand:
       try:
           del seeds[cand]
       except:
           continue
   return seeds

def Filter_variations(seeds, ploidy):
   for seed in seeds.keys():
       if '809015' in seed or '809022' in seed:
           test = seeds[seed]
       del_cand = []
       numb_of_vars = len(seeds[seed].keys())
       numb_of_reads = 0
       for var in seeds[seed].keys(): numb_of_reads += len(seeds[seed][var])
       for var in seeds[seed].keys():
           if len(seeds[seed][var]) < 5:
               del_cand.append(var)
       unique_del_cand = set(del_cand)
       for cand in unique_del_cand:
           del seeds[seed][cand]
   return seeds

def Get_Ploidy_Block_seeds(pb_start, pb_end, seeds):
   pb_seeds = {}
   for seed in seeds.keys():
       seed_s = int(seed.split(",")[0])
       seed_e = int(seed.split(",")[-1])
       if seed_s >= pb_start and seed_e <= pb_end:
           pb_seeds[seed] = seeds[seed]
   return pb_seeds

def Creating_Seeds(alignments):
   seeds = {}
   reads_to_seeds = {}
   for read in alignments.keys():
       alignment = alignments[read]
       for pair in ['pair_1', 'pair_2']:
           if len(alignment[pair]) < 2:
               continue
           coverage_pair1 = alignment[pair][3]
           if len(coverage_pair1) > 1 :
                   for start_i in range(2, len(coverage_pair1)+1):
                           for i in range(start_i, len(coverage_pair1)+1):
                                   pair_1 = coverage_pair1[i - start_i: i]
                                   key = ','.join([str(e) for e in pair_1])
                                   alleles_pos_1 = [e for e in pair_1]
                                   alleles_1 = Pos_to_Nuc_bam(alignment[pair], alleles_pos_1)
                                   if alleles_1 == None:
                                       continue
                                   if not key in seeds:
                                       seeds[key] = {}
                                   inside_key = alleles_1
                                   if not inside_key in seeds[key].keys():
                                       seeds[key][inside_key] = []
                                   seeds[key][inside_key].append(read)
                                   if not read in reads_to_seeds.keys():
                                       reads_to_seeds[read] = {}
                                   reads_to_seeds[read][key] = inside_key
   print(len(reads_to_seeds.keys()), "Total number of seeds")
   return seeds, reads_to_seeds