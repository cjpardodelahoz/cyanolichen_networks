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

## Structure of the regional network (Fig. 1?)

We first quantified the structure of the regional network using... This will take a while (several hours)

```sh
Rscript scripts/ecology/quantify_regional_network_structure.R
```

Now we plot the regional network highlighting compound structure. This will generate two alternatives: one with *Nostoc* in the columns and one with *Nostoc* in the rows.

```sh
Rscript scripts/ecology/plot_compound_network.R
```

We used the module assignments from the compound network plot to map the pairs in the map of Alberta:

```sh
Rscript scripts/ecology/plot_site_maps.R
```

## Vertical vs. Horizontal transmission

### Inference of clonality in *Peltigera* from ITS and *rbcLX* sequence data


```sh
Rscript scripts/ecology/infer_clones.R
```



