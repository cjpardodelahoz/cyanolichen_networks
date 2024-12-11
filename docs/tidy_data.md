# Data formatting and cleaning

## Regional (ABMI) datasets

We needed to clean, integrate, and summarize site and symbiont ID data from the ABMI (regional) sampling.

###**Site data**
We have (i) remote sensing data on climatic variables and forest cover provided directly by the ABMI (quadrant level), and (ii) ecosite classification and canopy cover data downloaded from the [ABMI website](https://abmi.ca/abmi-home/data-resources/data-portal-main/rawdata.html) on December 9, 2024. We integrated those datasets, summarized them at the site level, and filtered the to include the site-years included in our cyanolichen sampling.

###**Symbiont ID data**
We published the majority of the symbiont ID data (from 2,316 specimes) in a previous paper on *Nostoc* evolution. However, there are an additional 392 nrITS and 10 **rbcLX** sequences that we are reporting in this study that were part of the same sampling but were not published before. The IDs from these datasets needed to be integrated.

Both of htese datasets were prepared for analyses and publication with this script:

```sh
Rscript scripts/tidy_data/clean_abmi_data.R
```