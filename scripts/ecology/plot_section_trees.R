#!/usr/bin/env Rscript

##### LOAD PACKAGES #####

# Load packages
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)

##### PELTIGERA TREE #####

# Outgroup
peltigera_outgroup <- c("L220_Nephroma_plumbeum",
                        "L746_Nephroma_expallidum",
                        "L197_Nephroma_helveticum",
                        "saccata_Solorinasaccata",
                        "L915_Solorina_simensis1",
                        "L1045_Solorina_simensis2",
                        "crocea_Solorinacrocea")

# Load the Peltigera tree
peltigera_tree <- read.tree("data/trees/peltigera_chagnon2019.tree") %>%
    drop.tip(peltigera_outgroup)


# Section nodes
hydro_node <- MRCA(peltigera_tree, c("hydrothyria", "gowardii"))
phlebia_node <- MRCA(peltigera_tree, c("venosa", "venosa1"))
chloro_node <- MRCA(peltigera_tree, c("nigripunctata", "latiloba"))
peltidea_node <- MRCA(peltigera_tree, c("frippii", "aphthosa4"))
polydactylon_node <- MRCA(peltigera_tree, c("scabrosa2", "pulverulenta1"))
horizontales_node <- MRCA(peltigera_tree, c("neckeri", "elisabethae"))
retifoveata_node <- MRCA(peltigera_tree, c("retifoveata", "sp13"))
canina_node <- MRCA(peltigera_tree, c("aubertii", "monticola9"))


peltigera_tree_plot <- ggtree(peltigera_tree) %>%
    scaleClade(canina_node, .3) %>%
    scaleClade(retifoveata_node, 1.3) %>%
    scaleClade(polydactylon_node, 0.3) %>%
    scaleClade(horizontales_node, 0.6) %>%
    scaleClade(phlebia_node, 2) %>%
    collapse(node = hydro_node, 'min', fill = "steelblue", alpha = 0.5) %>%
    collapse(node = chloro_node, 'min', fill = "#707870", alpha = 0.5) %>%
    collapse(node = phlebia_node, 'max', fill = "darkred", alpha = 0.5) %>%
    collapse(node = peltidea_node, 'min', fill = "darkgreen", alpha = 0.5) %>%
    collapse(node = polydactylon_node, 'min', fill = "darkorange", alpha = 0.5) %>%
    collapse(node = horizontales_node, 'min', fill = "darkblue", alpha = 0.5) %>%
    collapse(node = retifoveata_node, 'min', fill = "darkviolet", alpha = 0.5) %>%
    collapse(node = canina_node, 'min', fill = "darkcyan", alpha = 0.5)

ggsave(peltigera_tree_plot, filename = "documents/plots/peltigera_tree.pdf", 
    width = 3, height = 8)
