# Ecological analyses

First we are going to describe the pattern

Then We are going to answer the main question directly with:

(i) interaction modeling
(ii) Variance explained by E vs residual covariance with symbionts from JSDMs (GJAM, maybe HMSC)

Then we will test predictions from the mechanisms inferred with the modeling, including:

-Relationship between specialization and abundance (including deviations of interaction frequencies from relative abundance)
-Species turnover vs interaction turnover along environmental gradients (using local dataset)
-Frequency of vertical and horizontal transmission--at regional and local scale
-Distribution of non-lichenized Nostoc along env gradients and relationship to molecular traits (e.g., Nostoc scytonemic)
-Regional network structure change with increasing env heterogeneity

Not sure exactly how, but we have to test how the responses to the process are similar (or different) between Nostoc and Peltigera. Maybe it's all about the size of the niche hypervolume. This could be linked to the relationship between specialization and abundance bullet up there...

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

***Step 1.*** We determined whether the network was modular by calculating Barber modularity (Q) of the observed network and comparing it to the expectation from null matrices simulated under a modell (*vaznull*) that preserves the connectance and degree distribution of the empirical network. Because the networks were significantly modular (i.e., z-score > 2 and p < 0.01), we proceeded to step 2.

***Step 2.*** We quantified nestedness within (*Nsm*) and outside (*Ndm*) modules with the WNODA index. "Outside" modules refers to the adjacency areas of the modules sensu Felix et al. 2022.

***Step 3.*** We evaluated both *Nsm* and *Ndm* given the network's modular structure by computing z-scores with a null model that preserves the modular structure of the empirical network and the connectance within and outside of modules (equiprobable-restricted null). We also calculated z-scores with the proportional restricted null model, which preserves degree distribution in addition to the properties preserved by the quiprobable restricted null. A compound network will have *Nsm* and *Ndm* that are significantly higher than expected under the equiprobable null, and similar *Nsm* and *Ndm* to the expectation under the proportional null.

We did these calculations for the full regional dataset (including all cyanolichens), and a subset including only *Peltigera*, which represent the majority of specimens (2068 out of 2148). This will take a while (several hours):

```sh
Rscript scripts/ecology/quantify_regional_network_structure.R
```

The script will print the results (Table SX) to `documents/tables/regional_network_metrics.csv`. It also generates `.RData` files with the interactino matrices and the full set of network structure metrics (including results not reported from proportional-restricted null and WNODA for full matrix) under `analyses/ecology/*interaction_matrix.RData` and `analyses/ecology/*network_metrics.RData`

### Plot of *Peltigere-Nostoc* regional network and module distributions

We then plotted the regional network highlighting the compound structure we found above. This will generate two alternatives: one with *Nostoc* in the rows and one with *Nostoc* in the columns. The first one is shown in Fig. 1.

```sh
Rscript scripts/ecology/plot_compound_network.R
```

The script will also generate `analyses/ecology/peltigera_module_assignments.RData`, which contains the module assigments that match the compound plot and which we used to generate the maps with the module distributions in Fig. 1:

```sh
Rscript scripts/ecology/plot_module_maps.R
```

### Plots of symbiont section trees

We plotted trees highlighting the sections within both *Nostoc* (Pardo-De la Hoz et al., 2025) and *Peltigera* [Chagnon et al 2019](https://lutzonilab.org/wp-content/uploads/Chagnon_et_al-2019-Journal_of_Ecology.pdf):

```sh
Rscript scripts/ecology/plot_section_trees.R
```

## Vertical vs. Horizontal transmission

### Inference of clonality in *Peltigera* from ITS and *rbcLX* sequence data


```sh
Rscript scripts/ecology/infer_clones.R
```



