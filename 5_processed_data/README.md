# Data dictionary

The data dictionary applies to the three data sets present in the folder: `final_mat_global.csv`, `LSU_inperson_Prolific.csv`, `LSU_online_LSU_inperson.csv`, and `LSU_online_Prolific.csv`. For each data set, one line represents data from a single participant.

+ `ID`
	+ *description*: unique anoymized identifier of the participant
    + *format*: character
+ `quiet`
	+ *description*: participant's mean recall score in the silent condition (averaged across serial positions)
    + *format*: numeric
+ `steady`
	+ *description*: participant's mean recall score in the steady state condition (averaged across serial positions)
    + *format*: numeric
+ `changing`
	+ *description*:participant's mean recall score in the changing state condition (averaged across serial positions)
    + *format*: numeric
+ `testing`
	+ *description*: condition of testing
    + *format*: factor (online/in-person)
+ `sample`
	+ *description*: site of data collection
    + *format*:	factor (LSU/Prolific)
+ `finalAudioCheck`
	+ *description*: whether the participant passed the final audio check
    + *format*:	factor (correct = passed / incorrect = failed)
+ `help_person`
	+ *description*: whether the participant reported having help from another person during the experiment
    + *format*: boolean (TRUE/FALSE)
+ `external_help`
	+ *description*: whether the participant reported having external help during the experiment
    + *format*: boolean (TRUE/FALSE)
+ `aloud_rehearsal`
	+ *description*: whether the participant reported using aloud rehearsal during the experiment
    + *format*: boolean (TRUE/FALSE)
+ `sound_off`
	+ *description*: whether the participant reported turning off the sound during the experiment
    + *format*: boolean (TRUE/FALSE)
+ `unplug_hp`
	+ *description*: whether the participant reported having unplugged the audio device during the experiment
    + *format*: boolean (TRUE/FALSE)
+ `external_distraction`
	+ *description*: whether the participant reported having experienced external distraction during the experiment
    + *format*: boolean (TRUE/FALSE)
+ `device_used`
	+ *description*: type of device used by the participant to complete the experiment
    + *format*: factor (Desktop/Laptopt/Smartphone/Tablet)
+ `audio_device`
	+ *description*: type of audio device used by the participant to complete the experiment
    + *format*: factor (On-ear/Over-ear/In-ear)
+ `response_device`
	+ *description*: Type of device used to respond during the experiment
    + *format*: factor (Mouse/Trackpad/Touchscreen)
+ `motivation`
	+ *description*: degree of motivation during the experiment
    + *format*: factor (lowest/low/average/high/highest)
+ `concentration`
	+ *description*: degree of concentration during the experiment
    + *format*: factor (lowest/low/average/high/highest)
+ `location_time`
	+ *description*: local time at the moment of completing the experiment (enter as a free text input)
    + *format*: character
+ `switch_screen`
	+ *description*: whether the participant reported having switched between different during the experiment
    + *format*: boolean (TRUE/FALSE)
+ `technical_issue`
	+ *description*: whether the participant reported having experienced technical issues during the experiment
    + *format*: boolean (TRUE/FALSE)
+ `hearing`
	+ *description*: description of hearing difficulties, if any
    + *format*: character
+ `cse`
	+ *description*: the changing state effect (difference between `steady` and `changing`)
    + *format*: numeric
+ `sse`
	+ *description*: the steady state effect (difference between `quiet` and `steady`)
    + *format*: numeric
+ `ise`
	+ *description*: the irrelevant sound effect (difference between `quiet` and `changing`)
    + *format*: numeric