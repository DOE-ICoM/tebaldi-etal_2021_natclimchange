[![DOI](https://zenodo.org/badge/385398259.svg)](https://zenodo.org/badge/latestdoi/385398259)

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
The Paris agreement focused global climate mitigation policy on limiting global warming to 1.5&deg;C or 2&deg;C above pre-industrial. Consequently,  projections of hazards and risk are increasingly framed in terms of global warming levels (GWLs) rather than emission scenarios. Here, we use a multi-method approach to describe changes in extreme sea levels  (ESLs) driven by changes in mean sea level associated with a wide range of GWLs, from 1.5&deg;C to 5&deg;C, and for a large number of locations providing uniform coverage over most of the world's coastlines. 

We estimate that by 2100 approximately 50% of the 7,000+ locations considered will experience the present-day 100-yr ESL event at least once a year, even under 1.5&deg;C of warming, and often well before the end of the century. The tropics appear more sensitive than the Northern high latitudes, where some locations do not see this frequency change even for the highest GWLs.


## Journal reference
Tebaldi, C., and CoAuthors (2021). Extreme sea levels at different global warming levels. *Nature Climate Change*, 11, 746-751, https://doi.org/10.1038/s41558-021-01127-1.

## Code reference
Tebaldi, C., and CoAuthors (2021). Supporting code for Tebaldi et al. 2021 - Nature Climate Change [Code]. Zenodo. https://doi.org/10.5281/zenodo.5103345.

## Data reference
Tebaldi, C., and CoAuthors (2021). Supporting data for Tebaldi et al. 2021 - Nature Climate Change [Data set]. Zenodo. https://doi.org/10.5281/zenodo.5095675.

## Reproduce my results
1. Set up a working directory, with naming of your choice. Create three subdirectories named 'Rcode','Rdatasets' and 'pics'. Using the DOI link in the `Data reference` section above, download and unzip the input data into the working directory. You should see a directory called 'tebaldi-etal_2021_natclimchange_data' with three subdirectories named 'Kirezci', 'Vousdoukas', and 'Rasmussen'. Download the R scripts in the `workflow_and_figures` directory and store them in the directory named 'Rcode'.
2. In all but `allfunctions.r` and `convolve_alltheway_fordatapub.R`, change `filedir` to the path of your working directory if needed, unless you are running R from it already, in which case the current setting `filedir<-"./"` works.
3. Run the R scripts according to the following steps to process the data and reproduce the main figures in this publication:

| Step | Script Name | Description |
| --- | --- | --- |
| 1 | `Read_in_data.R` | Reads and restructures the CSV files into R arrays. The CSV files contain the ESL estimates from the corresponding three approaches, matched to the two alternative SLR projections, organized by the time horizon of the projection and the Global Warming Level.
| 2 | `probabilisticprojections_and_TWL100TWL1difference_fordatapub.r` | Applies the Fisher Information Matrix approach to the ESLs parameter estimates and convolves a sample from their distribution with a sample from the SLR projections; computes the difference between 100-yr and 1-yr events.
| 3 | `votingsystem_fordatapub.R` | Applies the voting system synthesis approach to the individual distribution to produce the main results of the paper, including part of the content in Table 1 and Figure 1.
| 4 | `convolve_alltheway_fordatapub.R` | Performs the full convolution as an alternative to the voting system. Produces the remaining content of Table 2, plots ED Figures 3 and 4.  Also performs analysis of timing of change in frequency, resulting in Table 2 and Supplementary Figures 12-19.

## Correction
The newest release of this metarepository addresses a bug we found after the paper and code was released. As described 
in the paper, we match two scenarios (2+ and 5+) that include the effects of ice-sheet melt and are 
available only within the projections obtained by the Rasmussen et al. 2018 method to the 2 degrees and 5 degrees 
scenarios from the other method, by Vega-Westhoff et al., 2019. A glitch in our original code, however, did not perform 
this matching correctly, and created NAs that – undetected -- affected the results from the voting system. The new code 
fixes this problem. The results of the paper are affected quantitatively but not qualitatively. In fact, the results 
from the voting system after the correction are closer to the results of the full convolution method, presented as an 
alternative approach in the paper. As of January 2022, we are in the process of posting a correction in the journal 
to document in detail all changes.