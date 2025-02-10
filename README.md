# Sex bias in mortality risk changes over the lifespan of bottlenose dolphins

## Description of the data and file structure

This file explains all of the variables in each of the datasets that accompany: McEntee, M.H.F.,  Foroughirad, V.,  Krzyszczyk, E., Mann, J. 2023. Sex bias in mortality risk changes over the lifespan of bottlenose dolphins. [doi/10.1098/rspb.2023.0675 ](https://doi.org/10.1098/rspb.2023.0675).

For more information please contact Molly McEntee at [mhm95@georgetown.edu](mailto:mhm95@georgetown.edu).

### Files and variables

#### File: lifespan\_data.csv

**Description:** data to produce Kaplan-Meier plots, Cox proportional hazards analysis, PAM analysis, and sensitivity analysis.

##### Variables

* Dolphin.ID = dolphin ID
* Birth.Date = date of dolphin birth (YYYY-MM-DD)
* Death.Date = date of dolphin death if known, NA if the dolphin is  alive (YYYY-MM-DD)
* Sex = Female, Male, or Unknown
* total.num.sightings = Total number of sightings of that dolphin in the dataset
* gap = maximum individual sighting gap (years), NA if the individual  occurred one time or fewer in the data.
* time = time at death if dead, time at right censorship if alive (years)
* event = 1 if dead, 0 if alive
* trunc.time = age the dolphin entered the dataset (years)
* first.sight.date = first occurrence of the dolphin in the dataset (YYYY-MM-DD)
* last.sight.date = last occurrence of the dolphin in the dataset (YYYY-MM-DD)

#### File: BaSTA\_data.csv

**Description:** data used for BaSTA analysis. This is mark-recapture data for the same juvenile and adult sample as the lifespan data. 

##### Variables

* ID =  dolphin ID
* Birth = year of birth, 0 if unknown. In this analysis we blinded the model to year of birth, so all values are 0.
* Death = year of death, 0 if unknown. In this analysis we blinded the model to year of death, so all values are 0.
* 1985-2019 = 1 if the dolphin was seen in that year, 0 if the dolphin was not seen in that year
* SexFemale = 1 if the dolphin is female, 0 if the dolphin is male
* SexMale = 1 if the dolphin is male, 0 if the dolphin is female

####

## Code/software

#### File: R.code.Lifespan.R

Code to reproduce all figures in R.
