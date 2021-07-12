_your zenodo badge here_

# tebaldi-etal_2021_natclimchange

**Extreme sea levels at different global warming levels**

Claudia Tebaldi<sup>1\*</sup>, Roshanka Ranasinghe<sup>2</sup>, Michalis Vousdoukas<sup>3</sup>, D.J. Rasmussen<sup>4</sup>, Ben Vega-Westhoff<sup>5</sup>, Ebru Kirezci<sup>6</sup>, Robert E. Kopp<sup>7</sup>, Ryan Sriver<sup>5</sup>, and Lorenzo Mentaschi<sup>3,8</sup>

<sup>1 </sup> Pacific Northwest National Laboratory, College Park, MD, USA  
<sup>2 </sup> IHE Delft Institute for Water Education, Delft, Netherlands  
<sup>3 </sup> European Commission, Joint Research Centre, Ispra, Italy  
<sup>4 </sup> Princeton University, Princeton, NJ, USA  
<sup>5 </sup> University of Illinois, Urbana-Champaign, IL, USA  
<sup>6 </sup> University of Melbourne, Melbourne, Australia  
<sup>7 </sup> Rutgers University, New Brunswick, NJ, USA  
<sup>8 </sup> University of Bologna, Bologna, Italy  

\* corresponding author: claudia.tebaldi@pnnl.gov

## Abstract
The Paris agreement focused global climate mitigation policy on limiting global warming to 1.5C or 2C above pre-industrial. Consequently, projections of hazards and risk are increasingly framed in terms of global warming levels (GWLs) rather than emission scenarios. Here, for the first time, we describe changes in extreme sea levels (ESLs) driven by changes in mean sea level associated with a wide range of GWLs, from 1.5C to 5C, and for a large number of locations providing uniform coverage over most of the world's coastlines. Our novel multi-method approach estimates that by 2100 approximately 50% of the 7,000+ locations considered will experience the present-day 100-yr ESL event at least once a year, even under 1.5C of warming, and often well before the end of the century. The tropics appear more sensitive than the Northern high latitudes, where some locations do not see this frequency change even for the highest GWLs.

## Journal reference
Tebaldi, C., & CoAuthors (2021). Extreme sea levels at different global warming levels. Nature Climate Change

## Code reference
Tebaldi, C., & CoAuthors (2021). Supporting code for Tebaldi et al. 2021 - Nature Climate Change [Code]. Zenodo. TBD

## Data reference

### Input data
TBD

### Output data
TBD

## Contributing modeling software
| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| model 1 | version | link to code repository | link to DOI dataset |

## Reproduce my experiment
Fill in detailed info here or link to other documentation that is a thorough walkthrough of how to use what is in this repository to reproduce your experiment.

1. Install the software components required to conduct the experiment from [Contributing modeling software](#contributing-modeling-software)
2. Download and install the supporting input data required to conduct the experiment from [Input data](#input-data)
3. Run the following scripts in the `workflow` directory to re-create this experiment:

| Script Name | Description | How to Run |
| --- | --- | --- |
| `step_one.py` | Script to run the first part of my experiment | `python3 step_one.py -f /path/to/inputdata/file_one.csv` |
| `step_two.py` | Script to run the last part of my experiment | `python3 step_two.py -o /path/to/my/outputdir` |

4. Download and unzip the output data from my experiment [Output data](#output-data)
5. Run the following scripts in the `workflow` directory to compare my outputs to those from the publication

| Script Name | Description | How to Run |
| --- | --- | --- |
| `compare.py` | Script to compare my outputs to the original | `python3 compare.py --orig /path/to/original/data.csv --new /path/to/new/data.csv` |

## Reproduce my figures
Use the scripts found in the `figures` directory to reproduce the figures used in this publication.

| Script Name | Description | How to Run |
| --- | --- | --- |
| `generate_figures.py` | Script to generate my figures | `python3 generate_figures.py -i /path/to/inputs -o /path/to/outuptdir` |
