# Metadata for `aphid_data_final.csv`

`week` - Week of experiment. Insect abundances were measured every 7 days. `week = 0` corresponds to the week when aphids were first added to plants.

`ID` - Unique label for each experimental mesocosm. Numbers correspond to the cage number (1-30) and letters correspond to the composition of aphid color morphs (R, Y, G) and whether the cage eventually had a parasitoid added to it (P).

`R` - Number of counter 'clicks' for red aphid morphs observed on a given number of `plants` surveyed.

`Y` - Number of counter 'clicks' for yellow aphid morphs observed on a given number of `plants` surveyed.

`G` - Number of counter 'clicks' for green aphid morphs observed on a given number of `plants` surveyed.

`W` - Number of counter 'clicks' for winged (alate) adult aphids on a given number of `plants` surveyed. Color morphs cannot be reliably distinguished at this stage. This data was not used in the analyses for the manuscript "Empirically measuring eco-evolutionary stability".

`P` - Number of adult parasitoid individuals counted in the cage. This count is for the entire cage volume and not associated with the number of plants counted.

`M` - Number of aphid mummies (parasitoid pupae) counted (always to individual level) on a given number of `plants` surveyed.

`plants` - Number of 2-week old radish seedlings surveyed for aphid abundances and mummies. This number was used to extrapolate aphid abundances to their density per 16 plants (total number of seedlings in the cage).

`aphids_per_click` - Multiplier to accurately determine aphid abundances. A value of 1 means each counter 'click' corresponds to 1 aphid individual, whereas a value of 5 means each counter 'click' corresponds to 5 aphid individuals.

`fungus_present` - Presence/absence (1/0) of a entomopathogen of the aphids. Despite its presence late in the experiment, preliminary analyses suggested it was unimportant in influencing ecological or evolutionary dynamics. This data was not used in the analyses for the manuscript "Empirically measuring eco-evolutionary stability".

`plants_no_growth` - Qualitative indicator of whether plants had normal growth (0) or were unusually small (1) during aphid sampling.

`spider_contam` - Number of spiders observed in cage. These individuals were recorded and removed immediately.

`ptoid_contam` - Number of adult parasitoids observed in cages that should not have parasitoids. These individuals were removed immediately and the cage monitored multiple times a week to remove them. This data was not used in the analyses for the manuscript "Empirically measuring eco-evolutionary stability".

`mummy_contam` - Number of aphid mummies (parasitoid pupae) observed in cages that should not have parasitoids. These individuals were killed immediately and the cage monitored multiple times a week to remove them. This data was not used in the analyses for the manuscript "Empirically measuring eco-evolutionary stability".

`aphids_extinct` - Persistence/extinction (values = 0/1) of aphid population. Note that zero values of abundances had to be observed for two consecutive weeks for the previous week to be labeled as an extinction.

`ptoids_extinct` - Persistence/extinction (values = 0/1) of parasitoid population. Note that zero values of insect abundances had to be observed for two consecutive weeks for the previous week to be labeled as an extinction.

`ptoids_added` - Number of adult female parasitoids added to the cage in week 3 and 4 of the experiment.
