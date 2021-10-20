# ///////////////////////////////////////\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\#
#                  AUDITORY DISTRACTION IN SHORT-TERM MEMORY                    #
#               Comparison between online and in-person testing                 #
#                        |--- Main analysis script ---|                         # 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////////#
#                      **** Code author: Simon Gorin****                        #

# Research authors: Raoul Bell, Emily M. Elliott, Simon Gorin, John E. Marsh, and Nick Robinson (alphabetical order)

# For any question related to this code and/or analysis, feel free to contact Simon at gorinsimon [ at ] gmail.com

here::i_am("data_processing.R") # Initialized project root


# ==============================================================================#
####--------------------- 1 Loads libraries and stuffs ----------------------####
# ==============================================================================#

library(splithalf)   # provides function to analyze internal consistency
library(here)        # to work with relative path
library(effectsize)  # provides functions to easily compute effect sizes (from the 'easystats' package series)
library(ggplot2)     # provides easy way to creates nice graphics
library(tidyr)       # provides a set of functions that help you get to tidy data
library(dplyr)       # provides a set of functions for data manipulation and wrangling
library(lubridate)   # provides set of function to work with dates
library(readr)       # provides a fast and friendly way to read data
library(purrr)       # provides tools to work with function and vectors (essentially the map()
library(stringr)     # provides set of functions to work with strings
library(UpSetR)      # provides functions to generates UpSet plots (technique to visualize set intersections in a matrix)
library(gridExtra)   # provides functions to arrange multiple plots together
library(grid)        # provides functions allowing to add manually title to plots
library(BayesFactor) # provides functions to compute Bayes factors


# Sourcing the file "geom_flat_violin.R" is required to generate the figures because
# it is not implemented by default in ggplot2. The file can be downloaded from the
# following link (https://gist.github.com/dgrtwo/eb7750e74997891d7c20).
# Be aware that there is an error in the original file that makes it crashes.
# The fix is indicated in the comment from the following link (https://gist.github.com/dgrtwo/eb7750e74997891d7c20#gistcomment-1974780)
# Finally, check that "geom_flat_violin.R" is located in the "codes" folder of your R project,
# the "codes" folder being located at the root of the project.
source(here::here("3_codes", "geom_flat_violin.R"))

# Sourcing the file "cousineau_morey_within_CI.R" to computes corrected confidence intervals for
# within-participants design.
# See Baguley (2012, https://doi.org/10.3758/s13428-011-0123-7) for more details
# Finally, check that "credible_interval_Baguley.R" is located in the "codes" folder of your R project
source(here::here("3_codes", "credible_interval_Baguley.R"))


# ==============================================================================#
####-------------------------- 2 Raw files loading --------------------------####
# ==============================================================================#

# Source the code containing the format of the different variables for the data loaded below.
# Check that "col_spec.R" is located in the "codes" folder of your R project.
source(here::here("3_codes", "col_spec.R"))

# Loads the LSU data set
LSU <- read_csv(file = here("4_data/data_LSU.csv"), col_names = TRUE, col_types = col_spec) %>%
    mutate(sample = "LSU")
# Loads the Prolific data set
Prolific <- read_csv(file = here("4_data/data_prolific.csv"), col_names = TRUE, col_types = col_spec) %>%
    mutate(sample = "Prolific")
# Binds the LSU and Prolific data sets
data_set <- bind_rows(LSU, Prolific)

# Creates a data frame that will contain participants performance across the three conditions,
# the condition of testing, as well as their answers to post-experiment questionnaire.
final_mat <- tibble(
    .rows = length(unique(data_set$code)),
    ID = as.character(NA), quiet = as.numeric(NA), steady = as.numeric(NA), changing = as.numeric(NA),
    testing = as.factor(NA), sample = as.factor(NA), finalAudioCheck = as.character(NA),
    help_person = as.logical(NA), external_help = as.logical(NA), aloud_rehearsal = as.logical(NA),
    sound_off = as.logical(NA), unplug_hp = as.logical(NA), external_distraction = as.logical(NA),
    device_used = as.character(NA), audio_device = as.character(NA), response_device = as.character(NA),
    motivation = as.factor(NA), concentration = as.factor(NA), location_time = as.character(NA),
    switch_screen = as.logical(NA), technical_issue = as.logical(NA), hearing = as.logical(NA)
)

# Set the levels of some factors
levels(final_mat$motivation) <- c("lowest", "low", "average", "high", "highest")
levels(final_mat$concentration) <- c("lowest", "low", "average", "high", "highest")

# Creates a data frame that will contains the same information as in "final_mat"
# but including trial performance level and having all the data in a long format.
reliabilityMat <- tibble(
    .rows = length(unique(data_set$code)) * 20 * 3,
    ID = as.character(NA),
    accuracy = as.numeric(NA),
    state = as.factor(NA),
    testing = as.factor(NA),
    finalAudioCheck = as.character(NA),
    help_person = as.logical(NA),
    external_help = as.logical(NA),
    aloud_rehearsal = as.logical(NA),
    sound_off = as.logical(NA),
    unplug_hp = as.logical(NA),
    external_distraction = as.logical(NA),
    device_used = as.character(NA),
    audio_device = as.character(NA),
    response_device = as.character(NA),
    motivation = as.factor(NA),
    concentration = as.factor(NA), 
    location_time = as.character(NA),
    switch_screen = as.logical(NA),
    technical_issue = as.logical(NA),
    hearing = as.logical(NA)
)

# Sets the levels of some factor variables
levels(reliabilityMat$state) <- c("quiet", "steady", "changing")
levels(reliabilityMat$testing) <- c("In-Person", "Online")
levels(reliabilityMat$motivation) <- c("lowest", "low", "average", "high", "highest")
levels(reliabilityMat$concentration) <- c("lowest", "low", "average", "high", "highest")


# ==============================================================================#
####--------------------------- 3 Data processing ---------------------------####
# ==============================================================================#

# Stores the participants ID in "final_mat" (from 1 to the number of participants)
final_mat$ID <- unique(data_set$code)

# Stores in the testing condition in "final_mat". # As for one of the two samples
# the testing condition was indicated as "online" with a lower case we need to recode
# the variable to ensure consistancy.
final_mat$testing <- recode(data_set$testing[match(final_mat$ID, data_set$code)], "online" = "Online")

# Stores in the sample type in "final_mat"
final_mat$sample <- data_set$sample[match(final_mat$ID, data_set$code)]

# Loops across all the participant IDs in "final_mat"
for (dl in 1:nrow(final_mat)) {
    
    # Creates from "final_mat" a sub-sample named "dataSub" that contains data from
    # the "dl" participant in "final_mat".
    dataSub <- filter(data_set, code == final_mat$ID[dl])
    
    # The lines below store information regarding the audio check and answers
    # to the post-experiment questionnaire.
    final_mat$finalAudioCheck[dl] <- as.character(dataSub$answer3[dataSub$sender == "Final audio check"])
    final_mat$help_person[dl] <- recode(dataSub$x1_did_you_have_any_help_from_another_person_when_remembering_the_digits[dataSub$sender == "Post-experiment Questions1"], "N" = FALSE, "Y" = TRUE)
    final_mat$external_help[dl] <- recode(dataSub$x2_did_you_use_any_external_help_e_g_paper_and_pencil_to_remember_the_digits[dataSub$sender == "Post-experiment Questions1"], "N" = FALSE, "Y" = TRUE)
    final_mat$aloud_rehearsal[dl] <- recode(dataSub$x3_did_you_say_the_digits_aloud_when_trying_to_remember_them[dataSub$sender == "Post-experiment Questions1"], "N" = FALSE, "Y" = TRUE)
    final_mat$sound_off[dl] <- recode(dataSub$x4_did_you_turn_off_the_volume_on_your_headphones_during_the_task[dataSub$sender == "Post-experiment Questions1"], "N" = FALSE, "Y" = TRUE)
    final_mat$unplug_hp[dl] <- recode(dataSub$x5_did_you_remove_or_unplug_your_headphones_during_the_task[dataSub$sender == "Post-experiment Questions1"], "N" = FALSE, "Y" = TRUE)
    final_mat$external_distraction[dl] <- recode(dataSub$x7_while_you_were_completing_the_study_were_there_any_external_sources_of_visual_or_auditory_distraction_e_g_other_people_speaking_in_the_same_room_a_running_video_a_song_playing_in_the_background_etc[dataSub$sender == "Post-experiment Questions1"], "N" = FALSE, "Y" = TRUE)
    final_mat$device_used[dl] <- as.character(dataSub$x8_what_equipment_did_you_use_to_do_the_experiment[dataSub$sender == "Post-experiment Questions1"])
    final_mat$audio_device[dl] <- as.character(dataSub$x9_what_type_of_headphones_did_you_use_to_play_the_sounds[dataSub$sender == "Post-experiment Questions2"])
    final_mat$response_device[dl] <- as.character(dataSub$x10_what_device_did_you_use_to_record_your_responses[dataSub$sender == "Post-experiment Questions2"])
    final_mat$motivation[dl] <- factor(c("lowest", "low", "average", "high", "highest"))[as.numeric(dataSub$x11_how_motivated_were_you_to_obtain_the_best_test_score_possible[dataSub$sender == "Post-experiment Questions2"])]
    final_mat$concentration[dl] <- factor(c("lowest", "low", "average", "high", "highest"))[as.numeric(dataSub$x12_how_concentrated_were_you_on_the_task[dataSub$sender == "Post-experiment Questions2"])]
    final_mat$location_time[dl] <- as.character(dataSub$x15_what_is_the_current_time_at_your_location_please_specify_am_or_pm[dataSub$sender == "Post-experiment Questions2"])
    final_mat$switch_screen[dl] <- recode(dataSub$x13_when_performing_the_task_were_you_switching_between_different_tasks_or_browsers[dataSub$sender == "Post-experiment Questions2"], "N" = FALSE, "Y" = TRUE)
    final_mat$technical_issue[dl] <- recode(dataSub$x14_did_you_experience_any_technical_difficulties_during_the_study_e_g_problems_with_the_internet_connection_delays_in_presentation_etc[dataSub$sender == "Post-experiment Questions2"], "N" = FALSE, "Y" = TRUE)
    final_mat$hearing[dl] <- as.character(dataSub$x16_if_you_reported_that_you_have_hearing_loss_at_the_start_of_the_study_please_tell_us_more_about_this_now[dataSub$sender == "Post-experiment Questions2"])
    
    # Filters the data to keep only the lines with a recall event while ignoring
    # training trials (experimental trials all have the variable "sender_id" starting
    # with the value indicated in the in the "sender_id" indicated in the line with
    # "sender" = "Outer Loop 20").
    dataSub <- dataSub %>%
        filter(
            sender == "Recall",
            str_sub(dataSub$sender_id, 1, 1) == dataSub$sender_id[dataSub$sender == "Outer Loop 20"]
        )
    
    # A matrix with the 60 target lists from the experiment
    spanMat <- map_dfr(
        map(dataSub$span, ~ as.numeric(str_split(.x, ",")[[1]])),
        ~ as.data.frame(rbind(.x))
    )
    # A matrix with the 60 lists recalled by the participant
    respMat <- map_dfr(
        map(dataSub$responses, ~ as.numeric(str_split(.x, ",")[[1]])),
        ~ as.data.frame(rbind(.x))
    )
    
    # Compares the two matrices to have a matrix of 0 (incorrect) and 1 (correct).
    # Note that the *1 is important to transform the logical format induced by the
    # comparison of "spanMat" and "respMat" into numeric.
    accuracyMat <- (spanMat == respMat) * 1
    
    # Computes the mean performance for the quiet condition (state_sim = "Q")
    final_mat$quiet[dl] <- mean(rowMeans(accuracyMat[dataSub$state_stim == "Q", ]))
    # Computes the mean performance for the steady-state condition (state_sim = "a")
    final_mat$steady[dl] <- mean(rowMeans(accuracyMat[dataSub$state_stim == "a", ]))
    # Computes the mean performance for the changing-state condition (state_sim = "Random")
    final_mat$changing[dl] <- mean(rowMeans(accuracyMat[dataSub$state_stim == "Random", ]))
    
    # The code below creates a data frame in which each row represents the response
    # from a participant to a single trial.
    
    # The index of the 60 lines for the ith participant
    dlR <- (((dl - 1) * 60) + 1):(((dl - 1) * 60) + 60)
    
    # Computes the mean recall accuracy for each trial
    reliabilityMat$accuracy[dlR] <- as.vector(rowMeans(accuracyMat))
    # Retrieves and stores the changing-state condition of each trial
    reliabilityMat$state[dlR] <- c("quiet", "steady", "changing")[match(dataSub$state_stim, c("Q", "a", "Random"))]
    
    # Retrieves and stores answers to the audio check and post experiment questionnaire
    reliabilityMat$finalAudioCheck[dlR] <- final_mat$finalAudioCheck[dl]
    reliabilityMat$help_person[dlR] <- final_mat$help_person[dl]
    reliabilityMat$external_help[dlR] <- final_mat$external_help[dl]
    reliabilityMat$aloud_rehearsal[dlR] <- final_mat$aloud_rehearsal[dl]
    reliabilityMat$sound_off[dlR] <- final_mat$sound_off[dl]
    reliabilityMat$unplug_hp[dlR] <- final_mat$unplug_hp[dl]
    reliabilityMat$external_distraction[dlR] <- final_mat$external_distraction[dl]
    reliabilityMat$switch_screen[dlR] <- final_mat$switch_screen[dl]
    reliabilityMat$technical_issue[dlR] <- final_mat$technical_issue[dl]
    reliabilityMat$testing[dl] <- final_mat$testing[dl]
    reliabilityMat$ID[dlR] <- final_mat$ID[dl]
    reliabilityMat$device_used[dlR] <- final_mat$device_used[dl]
    reliabilityMat$audio_device[dlR] <- final_mat$audio_device[dl]
    reliabilityMat$response_device[dlR] <- final_mat$response_device[dl]
    reliabilityMat$motivation[dlR] <- final_mat$motivation[dl]
    reliabilityMat$concentration[dlR] <- final_mat$concentration[dl]
    reliabilityMat$location_time[dlR] <- final_mat$location_time[dl]
}

# Computes the different distraction effects (CSE, ISE and ISE)
#|__ CSE = changing-state effect (steady - changing conditions)
#|__ SSE = steady-state effect (quiet - steady conditions)
#|__ ISE = irrelevant-sound effect (quiet - changing conditions)
final_mat <- final_mat %>%
    mutate(
        cse = steady - changing,
        sse = quiet - steady,
        ise = quiet - changing
    )


# ==============================================================================#
####---------------- 4 Filtering and writing processed data -----------------####
# ==============================================================================#

# Creating a data frame with participants from the LSU online sample and all
# 3xclusion criteria are applied (see the paper for more details).
final_mat_LSU_online <- final_mat %>%
    filter(
        testing == "Online",
        sample == "LSU",
        aloud_rehearsal == FALSE,
        finalAudioCheck == "correct",
        help_person == FALSE,
        external_help == FALSE,
        sound_off == FALSE,
        unplug_hp == FALSE,
        external_distraction == FALSE,
        switch_screen == FALSE,
        technical_issue == FALSE,
        !(device_used %in% c("Tablet", "Smartphone"))
    ) %>%
    filter(quiet >= (mean(quiet) - (3 * sd(quiet))))

# Creating a data frame with participants from the LSU in-person sample and all
# exclusion criteria are applied (see the paper for more details).
final_mat_LSU_inperson <- final_mat %>%
    filter(
        testing == "In-Person",
        sample == "LSU",
        aloud_rehearsal == FALSE,
        finalAudioCheck == "correct",
        help_person == FALSE,
        external_help == FALSE,
        sound_off == FALSE,
        unplug_hp == FALSE,
        external_distraction == FALSE,
        switch_screen == FALSE,
        technical_issue == FALSE,
        !(device_used %in% c("Tablet", "Smartphone"))
    ) %>%
    filter(quiet >= (mean(quiet) - (3 * sd(quiet))))

# Creating a data frame with participants from the Prolific sample and all
# exclusion criteria are applied (see the paper for more details).
final_mat_Prolific <- final_mat %>%
    filter(
        testing == "Online",
        sample == "Prolific",
        aloud_rehearsal == FALSE,
        finalAudioCheck == "correct",
        help_person == FALSE,
        external_help == FALSE,
        sound_off == FALSE,
        unplug_hp == FALSE,
        external_distraction == FALSE,
        switch_screen == FALSE,
        technical_issue == FALSE,
        !(device_used %in% c("Tablet", "Smartphone"))
    )


# The code below is required to incrementally check for outliers and replace them
# until we have 40 valid participants in the Prolific data set.

n40 <- 0 # status whether we have (1) or not (0) 40 valid participants after checking for outliers
id_out <- NULL # a vector containing the ID of the outliers added to the vector after each check

while (n40 == 0) {
    data_temp <- final_mat_Prolific %>%
        filter(!(ID %in% id_out)) %>%
        slice(1:40)
    
    outlier_vec <- data_temp$quiet >= (mean(data_temp$quiet) - (3 * sd(data_temp$quiet)))
    
    if (sum(outlier_vec) == 40) {
        n40 <- 1
        final_mat_Prolific <- data_temp
    } else {
        id_out <- c(id_out, data_temp$ID[!outlier_vec])
    }
}

# After filtering, binds back the  different samples for the main analysis together
final_mat_filtered <- bind_rows(final_mat_LSU_inperson, final_mat_LSU_online, final_mat_Prolific)

# Writes as a csv the filtered final matrix
write_csv(final_mat_filtered, here("5_processed_data", "final_mat_global.csv"))

# Writes as csv files paired combinations of data set (used to compare auditory
# distraction effects across samples)
write_csv(bind_rows(final_mat_LSU_inperson, final_mat_LSU_online), here("5_processed_data", "LSU_online_LSU_inperson.csv"))
write_csv(bind_rows(final_mat_LSU_online, final_mat_Prolific), here("5_processed_data", "LSU_online_Prolific.csv"))
write_csv(bind_rows(final_mat_LSU_inperson, final_mat_Prolific), here("5_processed_data", "LSU_inperson_Prolific.csv"))


# ==============================================================================#
# Below are created sub-samples to compare the effect of different filtering
# strategies on the effect size of auditory distraction in online samples
# ==============================================================================#

# Prolific sample

# Filters applied: none
final_mat_Prolific_full <- final_mat %>%
    filter(
        testing == "Online",
        sample == "Prolific"
    )

# Filters applied: no wrong device or technical issue
final_mat_Prolific_full_no_issue <- final_mat %>%
    filter(
        testing == "Online",
        sample == "Prolific",
        technical_issue == FALSE,
        !(device_used %in% c("Tablet", "Smartphone"))
    )

# Filters applied: same as above + missing audio check, external help, disabling sound
final_mat_Prolific_full_no_cheating <- final_mat %>%
    filter(
        testing == "Online",
        sample == "Prolific",
        finalAudioCheck == "correct",
        help_person == FALSE,
        external_help == FALSE,
        sound_off == FALSE,
        unplug_hp == FALSE,
        technical_issue == FALSE,
        !(device_used %in% c("Tablet", "Smartphone"))
    )

# Filters applied: same as above + switching screen and external distraction
final_mat_Prolific_full_no_distractionSwitch <- final_mat %>%
    filter(
        testing == "Online",
        sample == "Prolific",
        finalAudioCheck == "correct",
        help_person == FALSE,
        external_help == FALSE,
        sound_off == FALSE,
        unplug_hp == FALSE,
        external_distraction == FALSE,
        switch_screen == FALSE,
        technical_issue == FALSE,
        !(device_used %in% c("Tablet", "Smartphone"))
    )

# Filters applied: same as above + aloud rehearsal
final_mat_Prolific_full_no_rehearsal <- final_mat %>%
    filter(
        testing == "Online",
        sample == "Prolific",
        aloud_rehearsal == FALSE,
        finalAudioCheck == "correct",
        help_person == FALSE,
        external_help == FALSE,
        sound_off == FALSE,
        unplug_hp == FALSE,
        external_distraction == FALSE,
        switch_screen == FALSE,
        technical_issue == FALSE,
        !(device_used %in% c("Tablet", "Smartphone"))
    )

# LSU Online sample

# Filters applied: none
final_mat_LSU_full <- final_mat %>%
    filter(
        testing == "Online",
        sample == "LSU"
    )

# Filters applied: no wrong device or technical issue
final_mat_LSU_full_no_issue <- final_mat %>%
    filter(
        testing == "Online",
        sample == "LSU",
        technical_issue == FALSE,
        !(device_used %in% c("Tablet", "Smartphone"))
    )

# Filters applied: same as above + missing audio check, external help, disabling sound
final_mat_LSU_full_no_cheating <- final_mat %>%
    filter(
        testing == "Online",
        sample == "LSU",
        finalAudioCheck == "correct",
        help_person == FALSE,
        external_help == FALSE,
        sound_off == FALSE,
        unplug_hp == FALSE,
        technical_issue == FALSE,
        !(device_used %in% c("Tablet", "Smartphone"))
    )

# Filters applied: same as above + switching screen and external distraction
final_mat_LSU_full_no_distractionSwitch <- final_mat %>%
    filter(
        testing == "Online",
        sample == "LSU",
        finalAudioCheck == "correct",
        help_person == FALSE,
        external_help == FALSE,
        sound_off == FALSE,
        unplug_hp == FALSE,
        external_distraction == FALSE,
        switch_screen == FALSE,
        technical_issue == FALSE,
        !(device_used %in% c("Tablet", "Smartphone"))
    )


# ==============================================================================#
####---------------------------- 5 Effect sizes -----------------------------####
# ==============================================================================#

#### |__ 5.1 Auditory distractions effects ####

#### |_____ 5.1.1 Prolific sample ####

es_prolific_CSE <- cohens_d(
    x = final_mat_Prolific$steady,
    y = final_mat_Prolific$changing,
    paired = TRUE,
    ci = 0.95
)

es_prolific_SSE <- cohens_d(
    x = final_mat_Prolific$quiet,
    y = final_mat_Prolific$steady,
    paired = TRUE,
    ci = 0.95
)

es_prolific_ISE <- cohens_d(
    x = final_mat_Prolific$quiet,
    y = final_mat_Prolific$changing,
    paired = TRUE,
    ci = 0.95
)

#### |_____ 5.1.2 LSU online sample ####

es_LSU_online_CSE <- cohens_d(
    x = final_mat_LSU_online$steady,
    y = final_mat_LSU_online$changing,
    paired = TRUE,
    ci = 0.95
)

es_LSU_online_SSE <- cohens_d(
    x = final_mat_LSU_online$quiet,
    y = final_mat_LSU_online$steady,
    paired = TRUE,
    ci = 0.95
)

es_LSU_online_ISE <- cohens_d(
    x = final_mat_LSU_online$quiet,
    y = final_mat_LSU_online$changing,
    paired = TRUE,
    ci = 0.95
)

#### |_____ 5.1.3 LSU in-person sample ####

es_LSU_inperson_CSE <- cohens_d(
    x = final_mat_LSU_inperson$steady,
    y = final_mat_LSU_inperson$changing,
    paired = TRUE,
    ci = 0.95
)

es_LSU_inperson_SSE <- cohens_d(
    x = final_mat_LSU_inperson$quiet,
    y = final_mat_LSU_inperson$steady,
    paired = TRUE,
    ci = 0.95
)

es_LSU_inperson_ISE <- cohens_d(
    x = final_mat_LSU_inperson$quiet,
    y = final_mat_LSU_inperson$changing,
    paired = TRUE,
    ci = 0.95
)


#### |__ 5.2 Effect of filtering on auditory distraction effects ####

#### |_____ 5.2.1 Changing state effect ####

#### |________ 5.2.1.1 Prolific sample ####

# Filters applied: missing audio check, external help, disabling sound or wrong device used
es_prolific_full_CSE <- cohens_d(
    x = final_mat_Prolific_full$steady,
    y = final_mat_Prolific_full$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_prolific_full_no_issue_CSE <- cohens_d(
    x = final_mat_Prolific_full_no_issue$steady,
    y = final_mat_Prolific_full_no_issue$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_prolific_full_no_cheating_CSE <- cohens_d(
    x = final_mat_Prolific_full_no_cheating$steady,
    y = final_mat_Prolific_full_no_cheating$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: those described in the paper but not limitation to up to 40 participants
es_prolific_no_distractionSwitch_CSE <- cohens_d(
    x = final_mat_Prolific_full_no_distractionSwitch$steady,
    y = final_mat_Prolific_full_no_distractionSwitch$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: those described in the paper but not limitation to up to 40 participants
es_prolific_no_rehearsal_CSE <- cohens_d(
    x = final_mat_Prolific_full_no_rehearsal$steady,
    y = final_mat_Prolific_full_no_rehearsal$changing,
    paired = TRUE,
    ci = 0.95
)

# Table summary of Prolific sample
es_prolific_samples_comparison_CSE <- tibble(
    es = c(
        es_prolific_full_CSE$Cohens_d,
        es_prolific_full_no_issue_CSE$Cohens_d,
        es_prolific_full_no_cheating_CSE$Cohens_d,
        es_prolific_full_no_cheating_CSE$Cohens_d,
        es_prolific_no_rehearsal_CSE$Cohens_d,
        es_prolific_CSE$Cohens_d
    ),
    CI_low = c(
        es_prolific_full_CSE$CI_low,
        es_prolific_full_no_issue_CSE$CI_low,
        es_prolific_full_no_cheating_CSE$CI_low,
        es_prolific_no_distractionSwitch_CSE$CI_low,
        es_prolific_no_rehearsal_CSE$CI_low,
        es_prolific_CSE$CI_low
    ),
    CI_high = c(
        es_prolific_full_CSE$CI_high,
        es_prolific_full_no_issue_CSE$CI_high,
        es_prolific_full_no_cheating_CSE$CI_high,
        es_prolific_no_distractionSwitch_CSE$CI_high,
        es_prolific_no_rehearsal_CSE$CI_high,
        es_prolific_CSE$CI_high
    ),
    sample_type = c(
        "full",
        "full_no_issue",
        "full_no_cheating",
        "no distraction/switch screen",
        "no_rehearsal", "analyzed"
    ),
    sample_size = c(
        nrow(final_mat_Prolific_full),
        nrow(final_mat_Prolific_full_no_issue),
        nrow(final_mat_Prolific_full_no_cheating),
        nrow(final_mat_Prolific_full_no_distractionSwitch),
        nrow(final_mat_Prolific_full_no_rehearsal),
        nrow(final_mat_Prolific)
    ),
    BF10 = c(
        extractBF(ttestBF(
            x = final_mat_Prolific_full$steady,
            y = final_mat_Prolific_full$changing,
            paired = TRUE,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_issue$steady,
            y = final_mat_Prolific_full_no_issue$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_cheating$steady,
            y = final_mat_Prolific_full_no_cheating$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_distractionSwitch$steady,
            y = final_mat_Prolific_full_no_distractionSwitch$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_rehearsal$steady,
            y = final_mat_Prolific_full_no_rehearsal$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific$steady,
            y = final_mat_Prolific$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1]
    )
)

#### |________ 5.2.1.2 LSU online sample ####

# Filters applied: missing audio check, external help, disabling sound or wrong device used
es_LSU_full_CSE <- cohens_d(
    x = final_mat_LSU_full$steady,
    y = final_mat_LSU_full$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_LSU_full_no_issue_CSE <- cohens_d(
    x = final_mat_LSU_full_no_issue$steady,
    y = final_mat_LSU_full_no_issue$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_LSU_full_no_cheating_CSE <- cohens_d(
    x = final_mat_LSU_full_no_cheating$steady,
    y = final_mat_LSU_full_no_cheating$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: those described in the paper but not limitation to up to 40 participants
es_LSU_no_distractionSwitch_CSE <- cohens_d(
    x = final_mat_LSU_full_no_distractionSwitch$steady,
    y = final_mat_LSU_full_no_distractionSwitch$changing,
    paired = TRUE,
    ci = 0.95
)

# Table summary of LSU sample
es_LSU_samples_comparison_CSE <- tibble(
    es = c(
        es_LSU_full_CSE$Cohens_d,
        es_LSU_full_no_issue_CSE$Cohens_d,
        es_LSU_full_no_cheating_CSE$Cohens_d,
        es_LSU_full_no_cheating_CSE$Cohens_d,
        es_LSU_online_CSE$Cohens_d
    ),
    CI_low = c(
        es_LSU_full_CSE$CI_low,
        es_LSU_full_no_issue_CSE$CI_low,
        es_LSU_full_no_cheating_CSE$CI_low,
        es_LSU_no_distractionSwitch_CSE$CI_low,
        es_LSU_online_CSE$CI_low
    ),
    CI_high = c(
        es_LSU_full_CSE$CI_high,
        es_LSU_full_no_issue_CSE$CI_high,
        es_LSU_full_no_cheating_CSE$CI_high,
        es_LSU_no_distractionSwitch_CSE$CI_high,
        es_LSU_online_CSE$CI_high
    ),
    sample_type = c(
        "full",
        "full_no_issue",
        "full_no_cheating",
        "no distraction/switch screen",
        "analyzed"
    ),
    sample_size = c(
        nrow(final_mat_LSU_full),
        nrow(final_mat_LSU_full_no_issue),
        nrow(final_mat_LSU_full_no_cheating),
        nrow(final_mat_LSU_full_no_distractionSwitch),
        nrow(final_mat_LSU_online)
    ),
    BF10 = c(
        extractBF(ttestBF(
            x = final_mat_LSU_full$steady,
            y = final_mat_LSU_full$changing,
            paired = TRUE,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_full_no_issue$steady,
            y = final_mat_LSU_full_no_issue$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_full_no_cheating$steady,
            y = final_mat_LSU_full_no_cheating$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_full_no_distractionSwitch$steady,
            y = final_mat_LSU_full_no_distractionSwitch$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_online$steady,
            y = final_mat_LSU_online$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1]
    )
)


#### |_____ 5.2.2 Steady state effect ####

#### |________ 5.2.2.1 Prolific sample ####

# Filters applied: missing audio check, external help, disabling sound or wrong device used
es_prolific_full_SSE <- cohens_d(
    x = final_mat_Prolific_full$quiet,
    y = final_mat_Prolific_full$steady,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_prolific_full_no_issue_SSE <- cohens_d(
    x = final_mat_Prolific_full_no_issue$quiet,
    y = final_mat_Prolific_full_no_issue$steady,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_prolific_full_no_cheating_SSE <- cohens_d(
    x = final_mat_Prolific_full_no_cheating$quiet,
    y = final_mat_Prolific_full_no_cheating$steady,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: those described in the paper but not limitation to up to 40 participants
es_prolific_no_distractionSwitch_SSE <- cohens_d(
    x = final_mat_Prolific_full_no_distractionSwitch$quiet,
    y = final_mat_Prolific_full_no_distractionSwitch$steady,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: those described in the paper but not limitation to up to 40 participants
es_prolific_no_rehearsal_SSE <- cohens_d(
    x = final_mat_Prolific_full_no_rehearsal$quiet,
    y = final_mat_Prolific_full_no_rehearsal$steady,
    paired = TRUE,
    ci = 0.95
)

# Table summary of Prolific sample
es_prolific_samples_comparison_SSE <- tibble(
    es = c(
        es_prolific_full_SSE$Cohens_d,
        es_prolific_full_no_issue_SSE$Cohens_d,
        es_prolific_full_no_cheating_SSE$Cohens_d,
        es_prolific_full_no_cheating_SSE$Cohens_d,
        es_prolific_no_rehearsal_SSE$Cohens_d,
        es_prolific_SSE$Cohens_d
    ),
    CI_low = c(
        es_prolific_full_SSE$CI_low,
        es_prolific_full_no_issue_SSE$CI_low,
        es_prolific_full_no_cheating_SSE$CI_low,
        es_prolific_no_distractionSwitch_SSE$CI_low,
        es_prolific_no_rehearsal_SSE$CI_low,
        es_prolific_SSE$CI_low
    ),
    CI_high = c(
        es_prolific_full_SSE$CI_high,
        es_prolific_full_no_issue_SSE$CI_high,
        es_prolific_full_no_cheating_SSE$CI_high,
        es_prolific_no_distractionSwitch_SSE$CI_high,
        es_prolific_no_rehearsal_SSE$CI_high,
        es_prolific_SSE$CI_high
    ),
    sample_type = c(
        "full",
        "full_no_issue",
        "full_no_cheating",
        "no distraction/switch screen",
        "no_rehearsal", "analyzed"
    ),
    sample_size = c(
        nrow(final_mat_Prolific_full),
        nrow(final_mat_Prolific_full_no_issue),
        nrow(final_mat_Prolific_full_no_cheating),
        nrow(final_mat_Prolific_full_no_distractionSwitch),
        nrow(final_mat_Prolific_full_no_rehearsal),
        nrow(final_mat_Prolific)
    ),
    BF10 = c(
        extractBF(ttestBF(
            x = final_mat_Prolific_full$quiet,
            y = final_mat_Prolific_full$steady,
            paired = TRUE,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_issue$quiet,
            y = final_mat_Prolific_full_no_issue$steady,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_cheating$quiet,
            y = final_mat_Prolific_full_no_cheating$steady,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_distractionSwitch$quiet,
            y = final_mat_Prolific_full_no_distractionSwitch$steady,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_rehearsal$quiet,
            y = final_mat_Prolific_full_no_rehearsal$steady,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific$quiet,
            y = final_mat_Prolific$steady,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1]
    )
)


#### |________ 5.2.2.2 LSU online sample ####

# Filters applied: missing audio check, external help, disabling sound or wrong device used
es_LSU_full_SSE <- cohens_d(
    x = final_mat_LSU_full$quiet,
    y = final_mat_LSU_full$steady,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_LSU_full_no_issue_SSE <- cohens_d(
    x = final_mat_LSU_full_no_issue$quiet,
    y = final_mat_LSU_full_no_issue$steady,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_LSU_full_no_cheating_SSE <- cohens_d(
    x = final_mat_LSU_full_no_cheating$quiet,
    y = final_mat_LSU_full_no_cheating$steady,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: those described in the paper but not limitation to up to 40 participants
es_LSU_no_distractionSwitch_SSE <- cohens_d(
    x = final_mat_LSU_full_no_distractionSwitch$quiet,
    y = final_mat_LSU_full_no_distractionSwitch$steady,
    paired = TRUE,
    ci = 0.95
)

# Table summary of LSU sample
es_LSU_samples_comparison_SSE <- tibble(
    es = c(
        es_LSU_full_SSE$Cohens_d,
        es_LSU_full_no_issue_SSE$Cohens_d,
        es_LSU_full_no_cheating_SSE$Cohens_d,
        es_LSU_full_no_cheating_SSE$Cohens_d,
        es_LSU_online_SSE$Cohens_d
    ),
    CI_low = c(
        es_LSU_full_SSE$CI_low,
        es_LSU_full_no_issue_SSE$CI_low,
        es_LSU_full_no_cheating_SSE$CI_low,
        es_LSU_no_distractionSwitch_SSE$CI_low,
        es_LSU_online_SSE$CI_low
    ),
    CI_high = c(
        es_LSU_full_SSE$CI_high,
        es_LSU_full_no_issue_SSE$CI_high,
        es_LSU_full_no_cheating_SSE$CI_high,
        es_LSU_no_distractionSwitch_SSE$CI_high,
        es_LSU_online_SSE$CI_high
    ),
    sample_type = c(
        "full",
        "full_no_issue",
        "full_no_cheating",
        "no distraction/switch screen",
        "analyzed"
    ),
    sample_size = c(
        nrow(final_mat_LSU_full),
        nrow(final_mat_LSU_full_no_issue),
        nrow(final_mat_LSU_full_no_cheating),
        nrow(final_mat_LSU_full_no_distractionSwitch),
        nrow(final_mat_LSU_online)
    ),
    BF10 = c(
        extractBF(ttestBF(
            x = final_mat_LSU_full$quiet,
            y = final_mat_LSU_full$steady,
            paired = TRUE,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_full_no_issue$quiet,
            y = final_mat_LSU_full_no_issue$steady,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_full_no_cheating$quiet,
            y = final_mat_LSU_full_no_cheating$steady,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_full_no_distractionSwitch$quiet,
            y = final_mat_LSU_full_no_distractionSwitch$steady,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_online$quiet,
            y = final_mat_LSU_online$steady,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1]
    )
)


#### |_____ 5.2.3 Irrelevant sound effect ####

#### |________ 5.2.3.1 Prolific sample ####

# Filters applied: missing audio check, external help, disabling sound or wrong device used
es_prolific_full_ISE <- cohens_d(
    x = final_mat_Prolific_full$quiet,
    y = final_mat_Prolific_full$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_prolific_full_no_issue_ISE <- cohens_d(
    x = final_mat_Prolific_full_no_issue$quiet,
    y = final_mat_Prolific_full_no_issue$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_prolific_full_no_cheating_ISE <- cohens_d(
    x = final_mat_Prolific_full_no_cheating$quiet,
    y = final_mat_Prolific_full_no_cheating$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: those described in the paper but not limitation to up to 40 participants
es_prolific_no_distractionSwitch_ISE <- cohens_d(
    x = final_mat_Prolific_full_no_distractionSwitch$quiet,
    y = final_mat_Prolific_full_no_distractionSwitch$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: those described in the paper but not limitation to up to 40 participants
es_prolific_no_rehearsal_ISE <- cohens_d(
    x = final_mat_Prolific_full_no_rehearsal$quiet,
    y = final_mat_Prolific_full_no_rehearsal$changing,
    paired = TRUE,
    ci = 0.95
)

# Table summary of Prolific sample
es_prolific_samples_comparison_ISE <- tibble(
    es = c(
        es_prolific_full_ISE$Cohens_d,
        es_prolific_full_no_issue_ISE$Cohens_d,
        es_prolific_full_no_cheating_ISE$Cohens_d,
        es_prolific_full_no_cheating_ISE$Cohens_d,
        es_prolific_no_rehearsal_ISE$Cohens_d,
        es_prolific_ISE$Cohens_d
    ),
    CI_low = c(
        es_prolific_full_ISE$CI_low,
        es_prolific_full_no_issue_ISE$CI_low,
        es_prolific_full_no_cheating_ISE$CI_low,
        es_prolific_no_distractionSwitch_ISE$CI_low,
        es_prolific_no_rehearsal_ISE$CI_low,
        es_prolific_ISE$CI_low
    ),
    CI_high = c(
        es_prolific_full_ISE$CI_high,
        es_prolific_full_no_issue_ISE$CI_high,
        es_prolific_full_no_cheating_ISE$CI_high,
        es_prolific_no_distractionSwitch_ISE$CI_high,
        es_prolific_no_rehearsal_ISE$CI_high,
        es_prolific_ISE$CI_high
    ),
    sample_type = c(
        "full",
        "full_no_issue",
        "full_no_cheating",
        "no distraction/switch screen",
        "no_rehearsal", "analyzed"
    ),
    sample_size = c(
        nrow(final_mat_Prolific_full),
        nrow(final_mat_Prolific_full_no_issue),
        nrow(final_mat_Prolific_full_no_cheating),
        nrow(final_mat_Prolific_full_no_distractionSwitch),
        nrow(final_mat_Prolific_full_no_rehearsal),
        nrow(final_mat_Prolific)
    ),
    BF10 = c(
        extractBF(ttestBF(
            x = final_mat_Prolific_full$quiet,
            y = final_mat_Prolific_full$changing,
            paired = TRUE,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_issue$quiet,
            y = final_mat_Prolific_full_no_issue$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_cheating$quiet,
            y = final_mat_Prolific_full_no_cheating$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_distractionSwitch$quiet,
            y = final_mat_Prolific_full_no_distractionSwitch$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific_full_no_rehearsal$quiet,
            y = final_mat_Prolific_full_no_rehearsal$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_Prolific$quiet,
            y = final_mat_Prolific$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1]
    )
)

#### |________ 5.2.3.2 LSU online sample ####

# Filters applied: missing audio check, external help, disabling sound or wrong device used
es_LSU_full_ISE <- cohens_d(
    x = final_mat_LSU_full$quiet,
    y = final_mat_LSU_full$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_LSU_full_no_issue_ISE <- cohens_d(
    x = final_mat_LSU_full_no_issue$quiet,
    y = final_mat_LSU_full_no_issue$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: all from the previous sample as well as aloud rehearsal
es_LSU_full_no_cheating_ISE <- cohens_d(
    x = final_mat_LSU_full_no_cheating$quiet,
    y = final_mat_LSU_full_no_cheating$changing,
    paired = TRUE,
    ci = 0.95
)

# Filters applied: those described in the paper but not limitation to up to 40 participants
es_LSU_no_distractionSwitch_ISE <- cohens_d(
    x = final_mat_LSU_full_no_distractionSwitch$quiet,
    y = final_mat_LSU_full_no_distractionSwitch$changing,
    paired = TRUE,
    ci = 0.95
)

# Table summary of LSU sample
es_LSU_samples_comparison_ISE <- tibble(
    es = c(
        es_LSU_full_ISE$Cohens_d,
        es_LSU_full_no_issue_ISE$Cohens_d,
        es_LSU_full_no_cheating_ISE$Cohens_d,
        es_LSU_full_no_cheating_ISE$Cohens_d,
        es_LSU_online_ISE$Cohens_d
    ),
    CI_low = c(
        es_LSU_full_ISE$CI_low,
        es_LSU_full_no_issue_ISE$CI_low,
        es_LSU_full_no_cheating_ISE$CI_low,
        es_LSU_no_distractionSwitch_ISE$CI_low,
        es_LSU_online_ISE$CI_low
    ),
    CI_high = c(
        es_LSU_full_ISE$CI_high,
        es_LSU_full_no_issue_ISE$CI_high,
        es_LSU_full_no_cheating_ISE$CI_high,
        es_LSU_no_distractionSwitch_ISE$CI_high,
        es_LSU_online_ISE$CI_high
    ),
    sample_type = c(
        "full",
        "full_no_issue",
        "full_no_cheating",
        "no distraction/switch screen",
        "analyzed"
    ),
    sample_size = c(
        nrow(final_mat_LSU_full),
        nrow(final_mat_LSU_full_no_issue),
        nrow(final_mat_LSU_full_no_cheating),
        nrow(final_mat_LSU_full_no_distractionSwitch),
        nrow(final_mat_LSU_online)
    ),
    BF10 = c(
        extractBF(ttestBF(
            x = final_mat_LSU_full$quiet,
            y = final_mat_LSU_full$changing,
            paired = TRUE,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_full_no_issue$quiet,
            y = final_mat_LSU_full_no_issue$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_full_no_cheating$quiet,
            y = final_mat_LSU_full_no_cheating$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_full_no_distractionSwitch$quiet,
            y = final_mat_LSU_full_no_distractionSwitch$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1],
        extractBF(ttestBF(
            x = final_mat_LSU_online$quiet,
            y = final_mat_LSU_online$changing,
            paired = TRUE,
            rscale = sqrt(2) / 2,
            nullInterval = c(0, Inf)
        ))$bf[1]
    )
)


# ==============================================================================#
####---------------- 6 Figure (auditory distraction effect) -----------------####
# ==============================================================================#

# Pivots the tibble to a longer format to have the changing-state condition scores
# in two columns (state and accuracy) instead of three, which is easier to plot
# Also creates a group variable (testing_sample) to differentiate the three groups
final_mat_filtered_long <- final_mat_filtered %>%
    pivot_longer(quiet:changing, names_to = "state", values_to = "accuracy") %>%
    unite(testing_sample, testing:sample)

# Renames the factor levels of the 'state' variable
final_mat_filtered_long$state <- factor(recode(final_mat_filtered_long$state,
                                               "quiet" = "Silence",
                                               "steady" = "Steady",
                                               "changing" = "Changing"),
                                        levels = c("Silence",
                                                   "Steady",
                                                   "Changing"))

# Renames the factor levels of the 'testing_sample' to include line breaks (looks better in the plot)
final_mat_filtered_long$testing_sample <- factor(recode(final_mat_filtered_long$testing_sample,
                                                        "In-Person_LSU" = "Psychology students\n(in-person)\n",
                                                        "Online_LSU" = "Psychology students\n(online)\n",
                                                        "Online_Prolific" = "Online panel"),
                                                 levels = c("Psychology students\n(in-person)\n",
                                                            "Psychology students\n(online)\n",
                                                            "Online panel"))

# Plots recall accuracy (y axis) as a function of changing-state condition (color) and testing sample (x axis)
# Boxes range represent the range from first and third quartiles
# Vertical lines represent 1.5*IQR (inter-quartile range or the distance between 1st and 3rd quartiles)
# Thick vertical lines in the boxes are the median
# Points are individual data points
# Vertical curves is the distribution of the data for each condition
changing_state_by_group_plot <- final_mat_filtered_long %>%
    ggplot(aes(x = testing_sample, y = accuracy)) +
    geom_flat_violin(aes(fill = state), width = 0.75, position = position_nudge(x = .3, y = 0), adjust = 1.5, trim = TRUE, alpha = .5, colour = NA) +
    geom_boxplot(aes(fill = state), outlier.shape = NA, alpha = 0.3, width = 0.5) +
    geom_point(aes(color = state),
               alpha = 0.4,
               position = position_jitterdodge(jitter.width = .125, seed = 11689, dodge.width = 0.5)
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1)) +
    theme_classic() +
    labs(
        title = "Auditory distraction effects as a function\nof sample/procedure collection",
        x = NULL,
        y = "Recall accuracy (proportion)",
        fill = "Changing-state\ncondition",
        color = "Changing-state\ncondition"
    ) +
    theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 13),
        legend.text = element_text(size = 11)
    )

# Saves above plot
ggsave(
    filename = here("6_figures", "changing_state_by_group_plot.png"),
    plot = changing_state_by_group_plot,
    width = 9, height = 7
)


# ==============================================================================#
####-------------------- 7 Internal consistency analysis --------------------####
# ==============================================================================#

# Recodes the levels of the state variables in 'reliabilityMat'
reliabilityMat$state <- recode(reliabilityMat$state,
                               "quiet" = "Quiet",
                               "steady" = "Steady",
                               "changing" = "Changing"
)

#### |__ 7.1 Prolific samples ####

#### |_____ 7.1.1 Full sample ####

# Internal consistency analysis
int_cons_Prolific_full <- reliabilityMat %>%
    filter(ID %in% final_mat_Prolific_full$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |_____ 7.1.2 No issue ####

# Internal consistency analysis
int_cons_Prolific_full_no_issue <- reliabilityMat %>%
    filter(ID %in% final_mat_Prolific_full_no_issue$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |_____ 7.1.3 No cheating ####

# Internal consistency analysis
int_cons_Prolific_full_no_cheating <- reliabilityMat %>%
    filter(ID %in% final_mat_Prolific_full_no_cheating$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |_____ 7.1.4 No distraction ####

# Internal consistency analysis
int_cons_Prolific_full_no_distractionSwitch <- reliabilityMat %>%
    filter(ID %in% final_mat_Prolific_full_no_distractionSwitch$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |_____ 7.1.5 No rehearsal ####

# Internal consistency analysis
int_cons_Prolific_full_no_rehearsal <- reliabilityMat %>%
    filter(ID %in% final_mat_Prolific_full_no_rehearsal$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |_____ 7.1.6 Analyzed ####

# Internal consistency analysis
int_cons_Prolific <- reliabilityMat %>%
    filter(ID %in% final_mat_Prolific$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |__ 7.2 LSU online samples ####

#### |_____ 7.2.1 Full sample ####

# Internal consistency analysis
int_cons_LSU_full <- reliabilityMat %>%
    filter(ID %in% final_mat_LSU_full$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |_____ 7.2.2 No issue ####

# Internal consistency analysis
int_cons_LSU_full_no_issue <- reliabilityMat %>%
    filter(ID %in% final_mat_LSU_full_no_issue$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |_____ 7.2.3 No cheating ####

# Internal consistency analysis
int_cons_LSU_full_no_cheating <- reliabilityMat %>%
    filter(ID %in% final_mat_LSU_full_no_cheating$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |_____ 7.2.4 No distraction ####

# Internal consistency analysis
int_cons_LSU_full_no_distractionSwitch <- reliabilityMat %>%
    filter(ID %in% final_mat_LSU_full_no_distractionSwitch$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |_____ 7.2.5 Analyzed ####

# Internal consistency analysis
int_cons_LSU_online <- reliabilityMat %>%
    filter(ID %in% final_mat_LSU_online$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

#### |__ 7.3 LSU in-person sample ####

#### |_____ 7.3.1 Analyzed ####

# Internal consistency analysis
int_cons_LSU_inperson <- reliabilityMat %>%
    filter(ID %in% final_mat_LSU_inperson$ID) %>%
    splithalf(
        data = .,
        outcome = "accuracy",
        score = "average",
        conditionlist = unique(reliabilityMat$state),
        halftype = "random",
        permutations = 20000,
        var.ACC = "accuracy",
        var.condition = "state",
        var.participant = "ID",
        average = "mean",
        plot = FALSE
    )

# ==============================================================================#
####-------------- 8 Exclusion criteria intersections figures ---------------####
# ==============================================================================#

#### |__ 8.1 LSU online sample ####

png(here("6_figures", "intersection_criteria_LSU_online.png"), width = 7, height = 6, units = "in", res = 300)
final_mat %>%
    filter(sample == "LSU", testing == "Online") %>%
    select(finalAudioCheck:device_used, switch_screen:technical_issue) %>%
    mutate(across(finalAudioCheck:technical_issue, ~ if_else(.x %in% c(TRUE, "incorrect", "Tablet", "Smartphone"), 1, 0))) %>%
    rename(
        `Audio check` = finalAudioCheck,
        `Help from a person` = help_person,
        `External help` = external_help,
        `Overt rehearsal` = aloud_rehearsal,
        `Turning sound off` = sound_off,
        `Headphone plugged off` = unplug_hp,
        `External distraction` = external_distraction,
        `Wrong device` = device_used,
        `Switching screen` = switch_screen,
        `Technical issue` = technical_issue
    ) %>%
    as.data.frame(.) %>%
    upset(.,
          nsets = 10,
          nintersects = NA,
          mainbar.y.label = "Exclusion criteria intersections",
          text.scale = c(1.3, 1.3, 1, 1, 1.25, 1.25),
          order.by = "freq",
          sets = colnames(.),
          keep.order = TRUE
    )
grid.text("Psychology students\n(online)",x = 0.65, y=0.95, gp=gpar(fontsize=14))
graphics.off()

#### |__ 8.2 Prolific sample ####

png(here("6_figures", "intersection_criteria_Prolific.png"), width = 7, height = 6, units = "in", res = 300)
final_mat %>%
    filter(sample == "Prolific") %>%
    select(finalAudioCheck:device_used, switch_screen:technical_issue) %>%
    mutate(across(finalAudioCheck:technical_issue, ~ if_else(.x %in% c(TRUE, "incorrect", "Tablet", "Smartphone"), 1, 0))) %>%
    rename(
        `Audio check` = finalAudioCheck,
        `Help from a person` = help_person,
        `External help` = external_help,
        `Overt rehearsal` = aloud_rehearsal,
        `Turning sound off` = sound_off,
        `Headphone plugged off` = unplug_hp,
        `External distraction` = external_distraction,
        `Wrong device` = device_used,
        `Switching screen` = switch_screen,
        `Technical issue` = technical_issue
    ) %>%
    as.data.frame(.) %>%
    upset(.,
          nsets = 10,
          nintersects = NA,
          mainbar.y.label = "Exclusion criteria intersections",
          text.scale = c(1.3, 1.3, 1, 1, 1.25, 1.25),
          order.by = "freq",
          sets = colnames(.),
          keep.order = TRUE
    )
grid.text("Online panel",x = 0.65, y=0.95, gp=gpar(fontsize=14))
graphics.off()