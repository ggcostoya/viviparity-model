# A general model for the evolution of live birth in lizards

Guillermo Garcia-Costoya<sup>1</sup>, Saúl Domínguez-Guerrero<sup>2,3</sup>, Lutz Fromhage<sup>4</sup>, Matthew E. Gifford<sup>5</sup>, Eric A. Riddell<sup>6</sup>, Michael L. Logan<sup>1</sup>.

Affiliations:

<sup>1</sup>Department of Biology and Program in Ecology, Evolution, and Conservation Biology, University of Nevada - Reno, Reno, NV, 89557, USA.

<sup>2</sup>Department of Ecology and Evolutionary Biology, Yale University, New Haven, CT, 06511, USA.

<sup>3</sup>Yale Institute for Biospheric Studies, Yale University, New Haven, CT, 06511 USA.

<sup>4</sup>Department of Biological and Environmental Science, University of Jyvaskyla, Jyvaskyla, Finland.

<sup>5</sup>University of Central Arkansas, 201 Donaghey Ave., LSC 180, Conway, AR 72035, USA.

<sup>6</sup>Department of Biology, University of North Carolina - Chapel Hill, Chapel Hill, NC, 27599, USA.

### Abstract

Life-history traits evolve to optimize fitness trade-offs over ontogeny. Potentially to mediate the trade-off between survival and fecundity, multiple taxa have independently evolved live birth (viviparity), including at least 70 transitions in lizards alone. In lizards, viviparity is thought to evolve as a mechanism to improve embryonic development in cold climates (cold climate hypothesis, or CCH), potentially at the expense of the mother’s survival. Past comparative studies often align with the CCH’s predictions, but they usually treat core features of the hypothesis (e.g., the roles played by behavioral thermoregulation and variation in life-history traits) as implicit and, most importantly, typically infer process from pattern rather than testing causal mechanisms, leaving the CCH without unequivocal support. To address this, we developed a process-based model that integrates behavior, thermal physiology, life history, and climate to predict optimal gestation length. We generated a comprehensive trait database of 89 globally distributed lizard populations that vary in parity mode, and we used ecophysiological modelling to test our model’s predictive power. Our model produced generally accurate predictions, strongly supporting the hypothesis that sub-optimal climates favor the evolution of viviparity in lizards and revealing the ecological contexts and underlying mechanisms by which this life history strategy evolves.

### Repository structure

-   `data/`: Contains all data used in the manuscript.

    -   `adult_data.csv`: Population level data on adult characterigstics. Contains columns for parity mode, species, adult thermal physiology traits ($CT_{min}$, $T_{opt}$ and $CT_{max}$), life history traits ( $Z_a$, $Z_e$, $N$), activity hours (Diurnal, Nocturnal or Cathemeral), substrate use (Arboreal, Terrestrial, Fossorial, Saxicolous etc.). All data is associated with each population's location (latitude, longitude and elevation above sea level) and that populations timing of embryonic development and nesting conditions which are also contained in `embryo_dev_nest_conditions.csv`. This data also contains columns to determine the source of the adult thermal physiology data (`source`, either `b` for Buckley et al. 2023 or `dg` for Domínguez-Guerrero et al. 2022) and for life-history trait data (`source_lh`, either `dg` for Domínguez-Guerrero et al. 2022 or `sb` for the SquamBase in Meiri et al. 2023).

    -   `eco_data.csv`: Population level data simulated via ecophysiological modelling. Contains monthly information on $T_g$, $T_{ap}$ and $T_{ep}$ at different nest depths and shade coverage levels, together wit environmental temperature. All data is associated with each population's location (latitude, longitude and elevation above sea level) and an estimate of environmental temperature at a reference height of 2m under shade (`tenv`).

    -   `embryo_thermal_phys.csv`: Species level data on embryonic thermal physiology. Contains columns for parity mode, species, and embryonic $CT_{min}$, $T_{opt}$ and $CT_{max}$.

    -   `embryo_dev_nest_conditions.csv`: Population level information on breeding phenology and nesting conditions. Contains columns for parity mode, species and population location (latitude, longitude, elevation above sea level), together with columns indicating which month of the year are embryos first seen in female oviducts (`eggs_seen`), which month of the year are hatchlings first seen (`hatchlings_seen`), at what depth in the soil (`nest_depth`, in cm) and at what degree of shade does the species typically lay its eggs in the wild (`nest_shade` as a proportion of sun exposure). This data also contains information on the reference from which breeding phenology (`breed_ref`) and nesting conditions (`nest_ref`) as numbers that can be correlated with the `embryo_dev_nest_conditions_citations.xlsx` file. Lastly, for populations where no information on breeding phenology was available, we obtained data from the species that was both phylogenetically and geographically closest to the species of interest and this species is indicated in the column `breed_other_species` .

    -   `embryo_dev_nest_conditions_citations.xlsx`: Information on the references used to build the embryonic development and nesting conditions database. The file contains two tabs, one for embryonic development and one for nesting conditions. On each tab, there's a column indicating the reference number as used in `embryo_dev_nest_conditions.csv` and a column for the citation of the reference.

    -   `model_test_data.RData`: Population and month species data to test the predictive accuracy of our model. This data contains all information also contained on `adult_data.csv` and `eco_data.csv`. In addition, this data also contains columns to determine whether or not a month falls within the period when embryonic development is known to occur (`dev_check`, `0` if the month does not fall within the period and `1` if it does), for the parameters $\alpha$ and $\gamma$ and for the model's prediction ( $d^*$, `opt_d`) and fitness of that strategy. Each row corresponds to a test of the model for a given population, in a given month, under an assumed combination of nest depth, nest shade $\alpha$ and $\gamma$. This file is stored as an `.RData` object as opposed to a `.csv` file to accommodate for nested columns.

-   `figures/` :

-   `scripts/`:

    -   `adult_thermal_phys_life_history_traits_additional_info.R` : Combines data from Buckley et al. 2022, Domínguez-Guerrero et al. 2022 and the SquamBase (Meiri et al. 2023) to produce `adult_data.csv`.

    -   `ecophysiological_modelling.R` : Uses `adult_data.csv` and `embryonic_thermal_phys.csv` to perform the ecophysiological modelling using `NicheMapR` and produce `eco_data.csv`.

    -   `embryonic_thermal_physiology.R` : Combines data from Pettersen et al. 2023 the reptile development database (Noble et al. 2018) to produce `embryo_thermal_phys.csv`

    -   `fit_model_sensitivity_analysis.R` : Uses `adult_data.csv`, `eco_data.csv` and `embryonic_thermal_phys.csv` to produce `model_test_data.RData`.

    -   `plot_figures.R` :

    -   `summaries_and_statistics.R` :
