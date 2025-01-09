# Ecological analyses

First we are going to describe the pattern

Then We are going to answer the main question directly with:

(i) interaction modeling
(ii) Variance explained by E vs residual covariance with symbionts from JSDMs (GJAM, maybe HMSC)

Then we will test predictions from the mechanisms inferred with the modeling, including:

-Module climatic overlap vs module connectivity? This seems to explain connectivity between modules C, D, and E, but not module A
-Relationship between specialization and abundance (including deviations of interaction frequencies from relative abundance)
-Species turnover vs interaction turnover along environmental gradients (using local dataset)
-Frequency of vertical and horizontal transmission--at regional and local scale
-Distribution of non-lichenized Nostoc along env gradients and relationship to molecular traits (e.g., Nostoc scytonemic)
-Regional network structure change with increasing env heterogeneity

Not sure exactly how, but we have to test how the responses to the process are similar (or different) between Nostoc and Peltigera. Maybe it's all about the size of the niche hypervolume. This could be linked to the relationship between specialization and abundance bullet up there...

Evolutionary consequences

-Sister species differentiate along environmental but not symbiotic axes. Two examples with similar potential symbionts and which share when cooccur but don't always cooccur:
    - Peltigera leucophlebia 1 and Peltigera leucophlebia 2
    - Phylogroups V and XLII
    - Keep an eye out for malacea 2 vs malacea 5


Note for env abundance: rarefy and proportions to visualize relative abundances and check if sufficient depth was acquired. Consider count normalization for differential abundance analyses.
## The regional dataset

Include here a reference to the data clean-up:

```sh

```


Plot the map of Alberta with all sampled sites. Consider highlighting the 377 where Peltigera are present

```sh
Rscript scripts/ecology/plot_site_maps.R
```

## Structure of the regional network (Fig. 1 and Table SX)

### Quantification of network structure

We quantified network structure following the three-step procedure proposed by [Felix et al. 2022](https://doi.org/10.1111/oik.09538) to distinguish simple (i.e., modular or nested) from compound (modular with low-level nestedness) network structures.

***Step 1.*** We determined whether the network was modular by calculating Barber modularity (Q) of the observed network and comparing it to the expectation from null matrices simulated under a model (*vaznull*) that preserves the connectance and degree distribution of the empirical network. Because the networks were significantly modular (i.e., z-score > 2 and p < 0.01), we proceeded to step 2.

***Step 2.*** We quantified nestedness within (*Nsm*) and outside (*Ndm*) modules with the WNODA index. "Outside" modules refers to the adjacency areas of the modules sensu [Felix et al. 2022](https://doi.org/10.1111/oik.09538).

***Step 3.*** We evaluated both *Nsm* and *Ndm* given the network's modular structure by computing z-scores with a null model that preserves the modular structure of the empirical network and the connectance within and outside of modules (equiprobable-restricted null). We also calculated z-scores with the proportional restricted null model, which preserves degree distribution in addition to the properties preserved by the equiprobable restricted null. A compound network will have *Nsm* and *Ndm* that are significantly higher than expected under the equiprobable null, and similar *Nsm* and *Ndm* to the expectation under the proportional null.

We did these calculations for the full regional dataset (including all cyanolichens), and a subset including only *Peltigera*, which represent the majority of specimens (2068 out of 2148). This will take a while (~12-14 hrs):

```sh
Rscript scripts/ecology/quantify_regional_network_structure.R
```

The script will print the results (Table SX) to `documents/tables/regional_network_metrics.csv`. It also generates `.RData` files with the interaction matrices and the full set of network structure metrics (including results not reported from proportional-restricted null and WNODA for full matrix) under `analyses/ecology/*interaction_matrix.RData` and `analyses/ecology/*network_metrics.RData`

### Plot of *Peltigere-Nostoc* regional network and module distributions

We then plotted the regional network highlighting the compound structure we found above. This will generate two alternatives: one with *Nostoc* in the rows and one with *Nostoc* in the columns. The first one is shown in Fig. 1.

```sh
Rscript scripts/ecology/plot_compound_network.R
```

The script will also generate `analyses/ecology/peltigera_module_assignments.RData`, which contains the module assigments that match the compound plot and which we used to generate the maps with the module distributions in Fig. 1:

```sh
Rscript scripts/ecology/plot_module_maps.R
```
## Pairing mechanisms and drivers of network structure

### Spatial interaction modeling

This is the approach from [Gravel et al. (2019)](https://doi.org/10.1111/ecog.04006). The idea is to formalize the probability of observing an interaction between a pair of symbionts (i,j) at site y as a joint probability of two events conditional on the environmental conditions: 

-The probability of cooccurrence given the environmental conditions P(Xiy, Xjy | Ey) 
-The probability of interaction given cooccurrence and the environmental conditions P(Lijy | Xiy, Xjy Ey)

By using different expressions of these two terms, we can test hypothesis about the spatial mechanisms driving pairing in the network. We used four variants of the cooccurrence term and three variants of the interaction term (Table SX1) to produce six different models (Table SX2) that we fitted to our dataset.

***Table SX1***

| Term        | Label | Expression Variant         | Interpretation                                                                 |
|-------------|-------|----------------------------|---------------------------------------------------------------------------------|
| Coocurrence | C0    | P(Xiy, Xjy)                | Cooccurrence independent of E, i.e., constant across landscape without assuming independence between species |
| Coocurrence | C1    | P(Xiy) P(Xjy)              | Cooccurrence independent of E and independent between species                   |
| Coocurrence | C2    | P(Xiy, Xjy | Ey)           | Cooccurrence dependent on E and without assuming independence between species    |
| Coocurrence | C3    | P(Xiy|Ey) P(Xjy|Ey)        | Cooccurrence dependent on E and independent between species                     |
| Interaction | L0    | P(Lijy)                    | Deterministic, i.e., always 1 (if they interact whenever they cooccur) or always 0 (if they never interact) |
| Interaction | L1    | P(Lijy | Xiy, Xjy)         | Probabilistic, dependence only on cooccurrence                                  |
| Interaction | L2    | P(Lijy | Xiy, Xjy, Ey)     | Probabilistic, depends on both cooccurrence and E                               |


***Table SX2***

| Label | Coocurrence term       | Interaction term                   |
|-------|------------------------|------------------------------------|
| C0_L2 | P(Xiy, Xjy)            | P(Lijy | Xiy, Xjy)                 |
| C1_L2 | P(Xiy) P(Xjy)          | P(Lijy | Xiy, Xjy)                 |
| C2_L0 | P(Xiy, Xjy | Ey)       | P(Lijy)                           |
| C2_L1 | P(Xiy, Xjy | Ey)       | P(Lijy | Xiy, Xjy)                |
| C2_L2 | P(Xiy, Xjy | Ey)       | P(Lijy | Xiy, Xjy, Ey)            |
| C3_L2 | P(Xiy|Ey) P(Xjy|Ey)    | P(Lijy | Xiy, Xjy, Ey)            |

To evaluate the fit of the models, we generated a binary table where, for each site y, we recorded the observation of each species, Xiy and Xjy, their co-occurrence, Xijy, the observation of an interaction Lijy, and environmental co-variates Ey (Table X3 - this is an example of what the table looks like). Then, for each symbiont pair, we fitted the six model combinations (Table SX2) generalized linear models with a binomial error distribution and logit link function. We used the predicted probabilities to compute the likelihood of each observation given each model and compared them using AIC. 

```sh
sh scripts/ecology/interaction_modeling.R
```

We included only symbiont pairs that cooccurred at least 10 times, and generated figures for pairs that cooccurred at least 20 times (to ensure there was enough data to rejrect the more complex models) and which interacted at least once (42 symbiotic pairs). The models were also fitted separately for each of four environmental covariates (MAT, MAP, proportion of conifer land cover, and elevation). The scripts generates plots comparing the fit for between models C2_L2 and C2_L1 across the four E covariates (`documents/plots/delta_aic_ECOVARIATE.pdf`) as well as tables with a summary of the fit across all included pairs for each E covariate (`documents/tables/model_fit_summary_ECOVARIATE.csv`).

### Estimating species' responses to environmental variation with GJAM

We used GJAM to jointly model the distribution of Peltigera and Nostoc at the regional scale. As part of the model fitting, gjam estimates beta coefficients that represent the response of each taxon to the variation of environmental predictors included in the model. This way, we obtained posterior estimates of the responses to MAT, MAP, proportion of conifer land cover, and elevation (`documents/plots/ECOVARIATE_response_nostoc|peltigerac.pdf`). We also computed a PCA to visualize the similarities in environmental responses across taxa from different network modules (`documents/plots/env_response_pca.pdf`).

```sh
Rscript scripts/ecology/gjam_abmi.R
```

The model is data-hungry, so we only included taxa found in at least 30 sites. The script also writes the gjam output of the three chains to `analyses/ecology/gjam_output.RData`.

### Module spatial overlap vs connectivity

If symbiotic pairing depends mostly on cooccurrence, then modules that cooccurr more often should be more connected. We quantified cooccurrence of core symbiotic pairs between modules and calculated pairwise module connectivity as the sum of edge weights between the two modules divided by the total edge weight of the network. We then visualized a plot of pairwise module connectivity vs. module cooccurrence `documents/plots/module_connectivity.pdf`.

```sh
Rscript scripts/ecology/plot_module_connectivity.R
```

### Plots of symbiont section trees

We plotted trees highlighting the sections within both *Nostoc* (Pardo-De la Hoz et al., 2025) and *Peltigera* [Chagnon et al 2019](https://lutzonilab.org/wp-content/uploads/Chagnon_et_al-2019-Journal_of_Ecology.pdf):

```sh
Rscript scripts/ecology/plot_section_trees.R
```

## Vertical vs. Horizontal transmission

### Plot 16S placement with detection on environment/thalli samples and comparative genomics on motility and scytonemin

The hypothesis that interactions are largely driven by the effect of environmental variation on symbiont cooccurrence assumes that horizontal transmission of symbionts is common. For horizontal transmission to occur, *Nostoc* must be able to survive in non-lichenized states. Here I plot the phylogenetic placement of 16S ASVs of lichenized *Nostoc* (which include the great majority of *Nostoc* lineages in the regional dataset as well), and their detection in thalli and environmental samples across the 15 sites. XX of these 16S ASVs are identical to the 16S from genomes we have sequenced before. Therefore, I also mapped the number of copies of genes for hormogonia development (gvp) and scytonenim (scy) biosynthesis because these traits could explain why phylogroups such as VIIa are almost never found in the environment while other common lineages are.

```sh
Rscript scripts/ecology/plot_env_detection.R
```

### Inference of photobiont transmission mode in *Peltigera* from ITS and *rbcLX* sequence data

We can use the ITS-*rbcLX* haplotype combinations in the specimens as a proxy for their clonality: if two specimens have an identical combianation, we consider them clonal; if the combination is different, we assume there was an event of horizontal photobiont transmission. We used this principle to quantify clonality across *Peltigera* species at both local and regional scales. Specifically, we calculated the maximum haplotype combination proportion (MHP) for each species. A value of 1 means that all specimens of that species have the same haplotype combinations and are likely clonal. MHP decreases to 1/n, where n is the number of specimens, when every single specimen has a different combination, likely because horizontal photobiont transmission is common. 

At the local scale, the calculation is simply max. haplotype frequency/n. However, if we do the same at the regional scale, the result could be biased by situations where a species is very common and clonal in a few sites. Therefore, I generated 1000 pseudoreplicates of the regional scale populations where I randomly sampled a single specimen of a given species from every site where it was found. Thus, all pseudoreplicates included the same number of specimens, i.e., equal to the number of sites where the species was found. Then, I calculated MHP independently for each pseudoreplicate. Finally, I plotted both the distribution of MHP at the regional scale and the local MHP for species with more than 10 specimens at a site in the local dataset.

```sh
Rscript scripts/ecology/infer_clones.R
```