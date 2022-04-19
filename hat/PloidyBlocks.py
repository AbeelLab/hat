import numpy as np
from hat.Seeds import Seeds_by_ploidy, Sort_seeds_ploidy_loc

def Finding_ploidy_blocks(seeds, dist_parameter):
   ploidy_blocks = {}
   b = ""
   current_ploidy_level = 0
   highest_ploidy_seeds = ""
   for seed in seeds.keys():
       seed_vars = seed.split(",")
       if b == "":
           b = [seed_vars[0], seed_vars[-1]]
           current_ploidy_level = len(seeds[seed])
           highest_ploidy_seeds = seed
       else:
           if len(seeds[seed]) > current_ploidy_level:
               print(b)
               if abs(int(seed_vars[-1]) - int(b[0])) < dist_parameter:
                   current_ploidy_level = len(seeds[seed])
                   highest_ploidy_seeds = seed
                   b[1] = seed_vars[-1]
               else:
                   ploidy_blocks[','.join(b)] = current_ploidy_level
                   current_ploidy_level = len(seeds[seed])
                   highest_ploidy_seeds = seed
                   b = [seed_vars[0], seed_vars[-1]]
           elif len(seeds[seed]) == current_ploidy_level:
               b[1] = seed_vars[-1]
               highest_ploidy_seeds = seed
           else:
               if abs(int(seed_vars[-1]) - int(highest_ploidy_seeds.split(',')[0])) < dist_parameter:
                   b[1] = seed_vars[-1]
               else:
                   ploidy_blocks[','.join(b)] = current_ploidy_level
                   current_ploidy_level = len(seeds[seed])
                   highest_ploidy_seeds = seed
                   b = [seed_vars[0], seed_vars[-1]]
   if not ','.join(b) in ploidy_blocks.keys():
       ploidy_blocks[','.join(b)] = current_ploidy_level
   return ploidy_blocks

def Finding_ploidy_blocks_alg2(seeds, dist_parameter):
   sorted_seeds = Sort_seeds_ploidy_loc(seeds)
   ploidy_blocks = {}
   blocked_seeds = []

   for i in range(len(sorted_seeds)):
       if sorted_seeds[i] in blocked_seeds:
           continue
       last_start_high_p = int(sorted_seeds[i].split(",")[0])
       first_start_high_p = int(sorted_seeds[i].split(",")[0])
       start_b = int(sorted_seeds[i].split(",")[0])
       end_b = int(sorted_seeds[i].split(",")[-1])
       block_p = len(seeds[sorted_seeds[i]])
       blocked_seeds.append(sorted_seeds[i])
       for j in range(i, len(sorted_seeds)):
           if sorted_seeds[j] in blocked_seeds:
               continue
           if abs(int(sorted_seeds[j].split(",")[-1]) - last_start_high_p) < dist_parameter or abs(int(sorted_seeds[j].split(",")[-1]) - first_start_high_p) < dist_parameter:
               blocked_seeds.append(sorted_seeds[j])
               if int(sorted_seeds[j].split(",")[-1]) > end_b:
                   end_b = int(sorted_seeds[j].split(",")[-1])
               if int(sorted_seeds[j].split(",")[0]) < start_b:
                   start_b = int(sorted_seeds[j].split(",")[0])
               if len(seeds[sorted_seeds[j]]) == block_p:
                   if int(sorted_seeds[j].split(",")[0]) > last_start_high_p:
                       last_start_high_p = int(sorted_seeds[j].split(",")[0])
                   elif int(sorted_seeds[j].split(",")[0]) < first_start_high_p:
                       first_start_high_p = int(sorted_seeds[j].split(",")[0])
       ploidy_blocks[str(start_b) + ',' + str(end_b)] = block_p
   print(ploidy_blocks)
   return ploidy_blocks

def Finding_ploidy_blocks_alg3(seeds, dist_parameter):
   seeds_by_ploidy = Seeds_by_ploidy(seeds)
   ploidy_blocks = {}
   blocked_seeds = []
   block_seeds = []
   block_cores = []
   block_cores_additions = []
   for i in sorted(seeds_by_ploidy.keys(), reverse=True):
       block_cores = []
       block_cores_additions = []
       for k in range(len(seeds_by_ploidy[i])):
           if seeds_by_ploidy[i][k] in blocked_seeds:
               continue
           high_ploidy_block_seeds = [seeds_by_ploidy[i][k]]
           blocked_seeds.append(high_ploidy_block_seeds[0])
           if len(seeds_by_ploidy[i]) > 0:
               for seed in seeds_by_ploidy[i]:
                   if seed in blocked_seeds:
                       continue
                   distance_of_seeds = Distances_seed1_to_seeds(seed, high_ploidy_block_seeds)
                   if distance_of_seeds < dist_parameter:
                       high_ploidy_block_seeds.append(seed)
                       blocked_seeds.append(seed)
           block_cores.append(high_ploidy_block_seeds)
           block_cores_additions.append([])
       for j in seeds_by_ploidy.keys():
           if j >= i:
               continue
           if len(seeds_by_ploidy[j]) == 0:
               continue
           for seed in seeds_by_ploidy[j]:
               if seed in blocked_seeds:
                   continue
               dist_to_cores = Distances_seed_to_cores(seed, block_cores)
               if min(dist_to_cores) < dist_parameter:
                   min_idx = np.argmin(dist_to_cores)
                   block_cores_additions[min_idx].append(seed)
                   blocked_seeds.append(seed)
       if len(block_cores) == 0:
           continue
       for i in range(len(block_cores)):
           block_seeds.append(block_cores[i])
           for j in range(len(block_cores_additions[i])):
               block_seeds[i].append(block_cores_additions[i][j])
       for i in range(len(block_cores)):
           ploidy_b = Ploidy_block_from_seeds(block_seeds[i], seeds)
           ploidy_blocks[str(ploidy_b[0]) + "," + str(ploidy_b[1])] = ploidy_b[2]
       block_seeds = []
   return ploidy_blocks


def Distances_seed_to_cores(seed, block_cores):
   distances = []
   for core in block_cores:
       distances.append(Distances_seed1_to_seeds(seed, core))
   return distances

def Ploidy_block_from_seeds(blocked_seeds, seeds):
   start = 1000000000
   end = 0
   max_ploidy = 0
   for seed in blocked_seeds:
       if int(seed.split(",")[0]) < start:
           start = int(seed.split(",")[0])
       if int(seed.split(",")[-1]) > end:
           end = int(seed.split(",")[-1])
       if len(seeds[seed]) > max_ploidy:
           max_ploidy = len(seeds[seed])
   return [start, end, max_ploidy]


def Distances_seed1_to_seeds(seed_1, seeds):
   min_distance = 10000000000
   for seed in seeds:
        dist = Distance_two_seeds(seed_1, seed)
        if dist < min_distance:
            min_distance = dist
   return min_distance


def Distance_two_seeds(seed_1, seed_2):
   distance = 0
   start_1 = int(seed_1.split(',')[0])
   start_2 = int(seed_2.split(",")[0])
   end_1 = int(seed_1.split(",")[-1])
   end_2 = int(seed_2.split(",")[-1])
   if start_1 < start_2:
       distance = abs(start_1 - end_2)
   else:
       distance = abs(start_2 - end_1)
   return distance

def Ploidy_blocks_coverage(ploidy_blocks, var_locs, genome_size, ploidy):
    length_pb = []
    for ploidy_b in ploidy_blocks.keys():
        pb_s = int(ploidy_b.split(",")[0])
        pb_e = int(ploidy_b.split(",")[1])
        length_pb.append(pb_e - pb_s + 1)
        if ploidy_blocks[ploidy_b] > ploidy:
            ploidy_blocks[ploidy_b] = ploidy

    count_vars_inside_ploidy_blocks = 0
    for var in var_locs:
        for ploidy_b in ploidy_blocks.keys():
            pb_s = int(ploidy_b.split(",")[0])
            pb_e = int(ploidy_b.split(",")[1])
            if var >= pb_s and var <= pb_e:
                count_vars_inside_ploidy_blocks += 1
                break
    print("Var inside ploidy blocks", count_vars_inside_ploidy_blocks, "total vars", len(var_locs), "percentage",
          float(count_vars_inside_ploidy_blocks) / len(var_locs))
    print("total length pb", sum(length_pb), "percentage_chromosome", sum(length_pb) / genome_size,
          "average length ploidy blocks", np.mean(length_pb), "total number of ploidy blocks", len(length_pb))

def Get_Ploidy_Block_seeds(pb_start, pb_end, seeds):
   pb_seeds = {}
   for seed in seeds.keys():
       seed_s = int(seed.split(",")[0])
       seed_e = int(seed.split(",")[-1])
       if seed_s >= pb_start and seed_e <= pb_end:
           pb_seeds[seed] = seeds[seed]
   return pb_seeds