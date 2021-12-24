# README

This repository is associated to the study entitled *Auditory Distraction Can Be Studied Online! A Direct Comparison Between In-Person and Online Experimentation* and submitted to the *Journal of Cognitive Psychology* ([https://doi.org/10.1080/20445911.2021.2021924](https://doi.org/10.1080/20445911.2021.2021924)). The repository contains five different folders and a main script called `data_processing.R` that can be used to reproduce the analyses and figures reported in the manuscript.

The folder `1_material` contains the experimental task saved as a .json file created with lab.js ([https://lab.js.org/](https://lab.js.org/)). The folder stimuli contains the audio files named after the letter they correspond to (e.g., a.wav is the sound file palying the spoken letter A) and `sil.wav` is a silent audio file.

The folder `2_sample_size_simulations` contains the scripts used to perform the simulation used to determine the appropriate sample size for the expected effect size, as well as to generate the appropriate figures (for more details, see the manuscript).

The folder `3_data` contains the raw data, that is, data sets from the two location sites mentioned in the paper. The file `data_LSU.csv` contains data sets from the Psychology-students tested in-person and online. The file `data_Prolific.csv` contains data sets from the non-student sample tested online via Prolific Academic (see the manuscript for more details about the samples). These data are read by the main script. A data dictionary describing the variables and their format is also available in the folder.

The folder `4_codes` contains the scripts that are sourced by the main script `data_processing.R`. Here are listed the scripts sourced and what they are used for:

+ `col_spec.R` &rarr; used to format variables when reading the raw data
+ `credible_interval_Baguley.R` &rarr; used to computes credible intervals for repeated-measures using the method described in Baguley (2012, [https://doi.org/10.3758/s13428-011-0123-7](https://doi.org/10.3758/s13428-011-0123-7)))
+ `geom_flat_violin.R` &rarr; used to generates the vertically-displayed distributions in Figure 2 (see the manuscript)

The folder `5_processed_data` &rarr; contains the processed data generated when running the main script `data_processing.R` (the folder is already populated with the processed data but the files are generated each time the main script is ran). A data dictionary describing the variables and their format is also available in the folder. The folder contains the following files:

+ `final_mat_global.csv` combines the processed data from all the samples
+ `LSU_inperson_Prolific.csv` combines the processed data from the Psychology students tested in-person and the non-students tested online
+ `LSU_online_LSU_inperson.csv` combines the processed data from the Psychology students tested online and the Psychology students tested in-person
+ `LSU_online_Prolific.csv` combines the processed data from the Psychology students tested online and the non-students tested online

The folder `6_figures` contains the figures displayed in the manuscript that are generated when running the main script

The folder `7_statistical_files` contains the .jasp statistical files with the outcomes of the statistical analysis (for more details about JASP and how to install it, see [https://jasp-stats.org/](https://jasp-stats.org/)). Here are listed the files and their content (see manuscript for more details):

+ `main_analysis.jasp` &rarr; registered analysis of variance (note that the seed has been set to '11689' to make the results reproducible)
+ `LSU_inperson_Prolific.jasp` &rarr; paired-comparison between Psychology students tested in-person and non-students tested online
+ `LSU_online_LSU_inperson.jasp` &rarr; paired-comparison between Psychology students tested online and sychology students tested in-person
+ `LSU_online_Prolific.jasp` &rarr; paired comparison between Psychology students tested online and non-students tested online
