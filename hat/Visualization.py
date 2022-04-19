import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

def Draw_variations_on_chromosome(var_locs, genome_size):
   print(len(var_locs))
   fig, ax = plt.subplots()
   plt.xlim(0, genome_size)
   for xc in var_locs:
       plt.axvline(x=xc, lw=0.1)
   ax.set_frame_on(False)
   ax.tick_params(top=False)
   ax.tick_params(labeltop=False)
   ax.set_yticklabels([])
   ax.set_xticklabels([])
   plt.gcf().set_size_inches(40, 3)
   plt.savefig("variations_sim_ploidy3.png", transparent=True)

def Draw_seeds_on_chromosome(var_locs, seeds, ploidy_blocks, genome_size, output_dir, output_prefix):
   print(ploidy_blocks)
   colors_dict = {5: "#0571b0", 3: "#92c5de", 4: "#f4a582", 2: "#ca0020"}
   cmap = ListedColormap(['#0571b0', '#92c5de"', '#f4a582', '#ca0020'])
   fig, ax = plt.subplots()
   plt.xlim(0, genome_size)
   len_seeds = []
   for seed in seeds.keys():
       len_seeds.append(len(seeds[seed]))
       locs = seed.split(",")
       x = list(range(int(locs[0]), int(locs[1]) + 1))
       y = []
       for i in range(len(x)):
           y.append(len(seeds[seed]))
       plt.plot(x, y, color="#2a2a2a", markersize=15, marker='|')
   for xc in var_locs:
       plt.axvline(x=xc, ymin=0, ymax= 0.1, lw=0.1, color='#2a2a2a', alpha=0.9)

   bar_x = []
   bar_y = []
   bar_width = []
   bar_colors = []
   bar_h = []
   for ploidy_b in ploidy_blocks.keys():
       color = colors_dict[ploidy_blocks[ploidy_b]]
       bar_x.append(int(ploidy_b.split(',')[0]))
       bar_y.append(ploidy_blocks[ploidy_b] + 0.2)
       bar_width.append(int(ploidy_b.split(',')[1]) - int(ploidy_b.split(',')[0]))
       bar_colors.append(color)
       bar_h.append(-0.4)
   plt.barh(bar_y, height=bar_h, width=bar_width, left=bar_x, color=bar_colors, align='edge', linewidth=3, edgecolor=bar_colors)
   ax.set_frame_on(False)
   ax.set_ylabel("Combination of alleles", fontsize=36)
   ax.set_xlabel("Position in the chromosome", fontsize=26)
   print("max y =", max(len_seeds), len_seeds)
   ax.set_yticks([0,1,2,3,4,5,6])
   ax.set_ylim(-1, 6)
   plt.yticks(fontsize=32)
   plt.xticks(fontsize=32)
   plt.gcf().set_size_inches(40, 20)
   five_patch = mpatches.Patch(color='#0571b0', label='5')
   three_patch = mpatches.Patch(color='#92c5de', label='3')
   four_patch = mpatches.Patch(color='#f4a582', label='4')
   two_patch = mpatches.Patch(color='#ca0020', label='2')
   six_patch = mpatches.Patch(color='darkslateblue', label='6')
   plt.legend(handles=[two_patch, three_patch, four_patch, five_patch], loc='upper right', prop={'size': 26})
   plt.savefig(output_dir + output_prefix + "_seeds_on_chromosome.pdf")
