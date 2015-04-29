
library("ggplot2")

################################################################################
# Testing phyloP scores between observed and simulated values
################################################################################

setwd("/home/gu/projects/mirna_host_evolution/conservation/phylop_simulation")

# inter and intra flank
phylop_exp = read.table("mirna_phylop_vert_plac_prim_inter_intra_flank.txt", head=T) 
phylop_exp = read.table("mirna_phylop_vert_plac_prim_inter_intra_flank_anc_nov.txt", head=T) 
# inter flank and intra random
phylop_exp = read.table("mirna_phylop_vert_plac_prim_inter_flank_intra_r.txt", head=T)
phylop_exp = read.table("mirna_phylop_vert_plac_prim_inter_flank_intra_r_anc_nov.txt", head=T)
# inter random flank and intra random
phylop_exp = read.table("mirna_phylop_vert_plac_prim_inter_rf_intra_r.txt", head=T)
phylop_exp = read.table("mirna_phylop_vert_plac_prim_inter_rf_intra_r_anc_nov.txt", head=T)

# Final plot 
phylop_exp = read.table("mirna_phylop_vert_plac_prim_inter_rf_intra_r_anc_nov_range2.txt", head=T)


# intragenic sense and intergenic
phylop_exp = phylop_exp[(phylop_exp$type == "intragenic" & phylop_exp$direction == "sense") | phylop_exp$type == "intergenic", ]

# plot only 7_12!
phylop_exp = phylop_exp[phylop_exp$age == "7_12", ]

# join (1 and 2-4 to compare inter vs intra)
phylop_exp = phylop_exp[phylop_exp$age == "1" | phylop_exp$age == "2_4", ]

phylop_exp$age <- factor(phylop_exp$age, 
                         c("1", "2_4", "5_6", "7_12"))
#phylop_exp$obs_sim <- factor(phylop_exp$obs_sim, 
#								  c("inter_obs", "inter_sim", "intra_obs_anc"))

ggplot(phylop_exp, aes(x=obs_sim, phyloP_score)) +
  #facet_grid(.~age, space="free", scales="free") +
  # width and width controls width of bar and space between bars respectively!
  geom_boxplot(width=0.8, position = position_dodge(width=0.8),aes(fill=obs_sim)) +
               #outlier.shape=NA) +
  scale_fill_manual(values=c("#ca0020", "#ca0020", "#0571b0", "#0571b0", 
  						     "#7b3294", "#7b3294"), 
                    labels=c("Inter_obs", "Inter_sim", "Intra_anc_obs",
                    		 "Intra_anc_sim", "Intra_nov_obs", "Intra_nov_sim"),
  				              name="") + 
  xlab("") +
  ylab("phyloP") +
  theme_bw()


inter_obs_1 = phylop_exp[phylop_exp$obs_sim == "inter_obs" & phylop_exp$age == "2_4", ]
intra_obs_1_anc = phylop_exp[phylop_exp$obs_sim == "intra_anc_obs" & phylop_exp$age == "2_4", ]

wilcox.test(inter_obs_1$phyloP_score, intra_obs_1_anc$phyloP_score)

inter_obs_2_4 = phylop_exp[phylop_exp$obs_sim == "intra_nov_obs", ]
intra_obs_2_4_anc = phylop_exp[phylop_exp$obs_sim == "intra_nov_obs", ]
intra_obs_2_4_nov = phylop_exp[phylop_exp$obs_sim == "intra_nov_sim", ]

wilcox.test(inter_obs_2_4$phyloP_score, intra_obs_2_4_anc$phyloP_score)


phylop_exp$host_age_cat <- factor(phylop_exp$host_age_cat, 
						 c("intergenic", "ancient", "novel"))
ggplot(phylop_exp, aes(phyloP_score, fill = obs_sim)) + 
	   geom_density(alpha = 0.35) + facet_grid(.~host_age_cat) + theme_bw()

inter_obs_7_12 = phylop_exp[phylop_exp$obs_sim == "inter_obs", ]
inter_sim_7_12 = phylop_exp[phylop_exp$obs_sim == "inter_sim", ]

intra_obs_anc_7_12 = phylop_exp[phylop_exp$obs_sim == "intra_anc_obs", ]
intra_sim_anc_7_12 = phylop_exp[phylop_exp$obs_sim == "intra_anc_sim", ]

intra_obs_nov_7_12 = phylop_exp[phylop_exp$obs_sim == "intra_nov_obs", ]
intra_sim_nov_7_12 = phylop_exp[phylop_exp$obs_sim == "intra_nov_sim", ]

wilcox.test(intra_obs_anc_7_12$phyloP_score, intra_sim_anc_7_12$phyloP_score)


phylop_exp = phylop_exp[phylop_exp$age.1 == "1_4", ]

inter_obs_1_4 = phylop_exp[phylop_exp$obs_sim == "inter_obs", ]
intra_obs_1_4 = phylop_exp[phylop_exp$obs_sim == "intra_anc_obs", ]

wilcox.test(inter_obs_1_4$phyloP_score, intra_obs_1_4$phyloP_score)

