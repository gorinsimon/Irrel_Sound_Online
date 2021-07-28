# ///////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#                  AUDITORY DISTRACTION IN SHORT-TERM MEMORY                    #
#               Comparison between online and in-person testing                 #
#                       |--- Sample size simulations ---|                       # 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////////#
#                      **** Code author: Simon Gorin****                        #

# Research authors: Raoul Bell, Emily M. Elliott, Simon Gorin, John E. Marsh, and Nick Robinson (alphabetical order)

# For any question related to this code, feel free to contact Simon at gorinsimon [ at ] gmail.com

here::i_am("data_processing.R") # Initialized project root

# ==============================================================================#
####--------------------------- 1 Loads libraries ---------------------------####
# ==============================================================================#


# The BFDA package is no on CRAN but the development version can be installed via Github
# To install the development version, use the command lines below:
### library(devtools)
### install_github("nicebread/BFDA", subdir="package")
library(BFDA)      # provides a framework to run simulations for Bayesian sequential designs
library(here)      # to work with relative path
library(gridExtra) # provides functions to arrange multiple plots together
library(ggplot2)   # provides easy way to creates nice graphics

# For more details regarding the parameters used in the simulations
# see https://rawgit.com/nicebread/BFDA/master/package/doc/BFDA_manual.html


# ==============================================================================#
####----------------------- 2 H1 'world' simulations ------------------------####
# ==============================================================================#

sim.H1 <- BFDA.sim(expected.ES = 0.50,                      # simulated population effect size (design prior) set to cohen's d = 0.5
                   type = "t.paired",                       # paired t-test was used as the type of design in the simulation
                   prior = list("Cauchy",                   # Cauchy distribution was used for the prior
                                list(prior.location = 0,      # location of the distribution was centered on 0
                                     prior.scale = sqrt(2)/2  # r scale of sqrt(2)/2 was used for the prior
                                )),
                   n.min = 10,                              # initial sample size set to 10
                   n.max = 300,                             # maximum sample size set to 300
                   alternative = "greater",                 # one-sided t-test
                   boundary = Inf,                          # Bayes factor where a sequential run is stopped (default is Inf)
                   B = 10000,                               # number of simulated studies set to 10000
                   design = "sequential",                   # type of design set to sequential
                   verbose = TRUE,                          # whether the progress of the simulations is printed set to TRUE
                   stepsize = 10                            # the number of observation added to the sample after each step of the sequential process is 10
)


# ==============================================================================#
####----------------------- 3 H0 'world' simulations ------------------------####
# ==============================================================================#

sim.H0 <- BFDA.sim(expected.ES = 0,                        # simulated population effect size (design prior) set to cohen's d = 0
                   type = "t.paired",                      # paired t-test was used as the type of design in the simulation
                   prior = list("Cauchy",                  # Cauchy distribution was used for the prior
                                list(prior.location = 0,     # location of the distribution was centered on 0
                                     prior.scale = sqrt(2)/2 # r scale of sqrt(2)/2 was used for the prior
                                )),
                   n.min = 10,                             # initial sample size set to 10
                   n.max = 300,                            # maximum sample size set to 300
                   alternative = "greater",                # one-sided t-test
                   boundary = Inf,                         # Bayes factor where a sequential run is stopped (default is Inf)
                   B = 10000,                              # number of simulated studies set to 10000
                   design = "sequential",                  # type of design set to sequential
                   verbose = TRUE,                         # whether the progress of the simulations is printed set to TRUE
                   stepsize = 10)                          # the number of observation added to the sample after each step of the sequential process is 10 


# ==============================================================================#
####------------------------ 3 Simulations analysis -------------------------####
# ==============================================================================#

# Creates a data frame that will contains the outcome of the analysis of the H1 world simulations
# Columns refers maximum sample sizes (from 50 to 200)
# Rows correspond to the proportion of simulated studies reaching evidence for H0 and H1
mat_sampleMax_50to200_H1 <- matrix(NA, ncol = 16, nrow = 2)
colnames(mat_sampleMax_50to200_H1) <- c("N=50","N=60","N=70","N=80","N=90","N=100",
                                        "N=110","N=120","N=130","N=140","N=150",
                                        "N=160","N=170","N=180","N=190","N=200")
rownames(mat_sampleMax_50to200_H1) <- c("prop_H0","prop_H1")

# Generate a data frame that will contains outcome from the H0 world simulations by duplicating
# the empty data frame genrated above for H1 simulations
mat_sampleMax_50to200_H0 <- mat_sampleMax_50to200_H1

# Analysis of the simulations. Minimum sample size is always 40 and maximum sample
# size is increased from 50 to 200 by steps of 10. The boundary to determine whether
# there is evidence for H1 or H0 is set to BF = 7
for (i in seq(50,200,10)) {
    # Computes the proportion of simulated studies reaching evidence for H0 in H1 world simulations
    mat_sampleMax_50to200_H1["prop_H0", match(i,seq(50,200,10))] <- BFDA.analyze(sim.H1,
                                                                                 design = "sequential",
                                                                                 n.min = 40,
                                                                                 n.max = i,
                                                                                 boundary = 7)$lower.hit.frac
    
    # Computes the proportion of simulated studies reaching evidence for H1 in H1 world simulations
    mat_sampleMax_50to200_H1["prop_H1", match(i,seq(50,200,10))] <- BFDA.analyze(sim.H1,
                                                                                 design = "sequential",
                                                                                 n.min = 40, n.max = i,
                                                                                 boundary = 7)$upper.hit.frac
    
    # Computes the proportion of simulated studies reaching evidence for H0 in H0 world simulations
    mat_sampleMax_50to200_H0["prop_H0", match(i,seq(50,200,10))] <- BFDA.analyze(sim.H0,
                                                                                 design = "sequential",
                                                                                 n.min = 40,
                                                                                 n.max = i,
                                                                                 boundary = 7)$lower.hit.frac
    
    # Computes the proportion of simulated studies reaching evidence for H1 in H0 world simulations
    mat_sampleMax_50to200_H0["prop_H1", match(i,seq(50,200,10))] <- BFDA.analyze(sim.H0,
                                                                                 design = "sequential",
                                                                                 n.min = 40,
                                                                                 n.max = i,
                                                                                 boundary = 7)$upper.hit.frac
}

# ==============================================================================#
####------------------------ 4 Plots the simulations-------------------------####
# ==============================================================================#

# Reshape the H1 simulations analysis into long format
H1sim_plot_matrix <- as.data.frame(cbind(rbind(mat_sampleMax_50to200_H0[2,],
                                               mat_sampleMax_50to200_H1[2,]),
                                         c("H0", "H1")))
H1sim_plot_matrix[,1:16] <- lapply(H1sim_plot_matrix[,1:16], as.character)
H1sim_plot_matrix[,1:16] <- lapply(H1sim_plot_matrix[,1:16], as.numeric)
colnames(H1sim_plot_matrix)[17] <- "H_world"
H1sim_plot_matrix <- reshape(H1sim_plot_matrix, direction = "long", varying = list(colnames(H1sim_plot_matrix)[1:16]), v.names = "Proportion", timevar = "SampleSize")

# Plots the outcome of H1 simulations
H1sim_plot <- ggplot(H1sim_plot_matrix, aes(SampleSize, Proportion, group = H_world)) +
    geom_line(aes(linetype = H_world)) +
    geom_point(aes(fill = H_world), colour = "black", shape = 21, size = 1.5, stroke = 1) +
    scale_fill_manual(values = c("black", "white")) +
    coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(breaks = seq(0,1,0.1)) +
    scale_x_continuous(breaks = seq(1,16,1), labels = seq(50,200,10)) +
    labs(x = "Maximum sample size", y = "Proportion of studies\n terminating at a boundary", subtitle = "Cohen's d = 0.5 / Min. sample size: N = 40") +
    theme_classic() +
    theme_light() +
    theme(panel.border = element_rect(fill = NA)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 12),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          legend.box.background = element_rect(),
          legend.box.margin = margin(1,1,1,1),
          legend.title=element_blank(),
          legend.text = element_text(size = 10),
          legend.key.width = unit(1, "cm"),
          legend.position = c(0.9,0.5))

# Reshape the H0 simulations analysis into long format
H0sim_plot_matrix <- as.data.frame(cbind(rbind(mat_sampleMax_50to200_H0[1,],mat_sampleMax_50to200_H1[1,]),c("H0", "H1")))
H0sim_plot_matrix[,1:16] <- lapply(H0sim_plot_matrix[,1:16], as.character)
H0sim_plot_matrix[,1:16] <- lapply(H0sim_plot_matrix[,1:16], as.numeric)
colnames(H0sim_plot_matrix)[17] <- "H_world"
H0sim_plot_matrix <- reshape(H0sim_plot_matrix, direction = "long", varying = list(colnames(H0sim_plot_matrix)[1:16]), v.names = "Proportion", timevar = "SampleSize")

# Plots the outcome of H0 simulations
H0sim_plot <- ggplot(H0sim_plot_matrix, aes(SampleSize, Proportion, group = H_world)) +
    geom_line(aes(linetype = H_world)) +
    geom_point(aes(fill = H_world), colour = "black", shape = 21, size = 1.5, stroke = 1) +
    scale_fill_manual(values = c("black", "white")) +
    coord_cartesian(ylim = c(0,1)) +
    scale_y_continuous(breaks = seq(0,1,0.1)) +
    scale_x_continuous(breaks = seq(1,16,1), labels = seq(50,200,10)) +
    labs(x = "Maximum sample size", y = "Proportion of studies\n terminating at a boundary", subtitle = "Cohen's d = 0 / Min. sample size: N = 40") +
    theme_classic() +
    theme_light() +
    theme(panel.border = element_rect(fill = NA)) +
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 12),
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 12),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          plot.subtitle = element_text(size = 13, hjust = 0.5, face = "bold"),
          legend.title=element_blank(),
          legend.text = element_text(size = 10),
          legend.key.width = unit(1, "cm"),
          legend.position = "none")

# Combines together and saves the two plots
final_figures <- grid.arrange(H1sim_plot,H0sim_plot, layout_matrix = rbind(1,2))
ggsave(here("2_sample_size_simulations", "ISE_oneline_sampleSize_figure.png"), plot = final_figures, width = 15, height = 14, units = "cm")
