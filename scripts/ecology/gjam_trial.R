#!/usr/bin/env Rscript

# Load the required libraries
library(gjam)
library(ggplot2)
library(ggpubr)
library(knitr)
library(cowplot)

# Load the data
load("birdData.rdata")

# create base for plots
base <- ggplot() +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

# Example of point 0 mass distribution
base +
  geom_histogram(aes(x = ydata$`american crow`)) +
  labs(title = "Distribution of American crow",
       y = "Frequency", x = "Counts of American crow")

# Variation of effort across sites
base +
  geom_histogram(aes(x = edata)) +
  labs(title = "Distribution of effort",
       y = "Frequency", x = "Effort (minutes)")

# Keep top 50 most abundant species - HOW TO SUBSAMPLE ALBERTA DATA?
ydata <- as.data.frame(gjamTrimY(ydata[,1:ncol(ydata)], maxCols = 50)$y)
ydata <- sapply(ydata, as.numeric)
ydata <- as.data.frame(ydata)

# In xdata, landcover should be a factor. Convert it to a factor and set the reference level—
# here, we set “barren” as the reference level. GJAM will fit coefficients for each factor 
#relative to the reference level (e.g. forest has X effect on abundance relative to barren).
xdata$landcover <- as.factor(xdata$landcover)
xdata$landcover <- relevel(xdata$landcover, 
                           "barren")

# The first is a vector of length 𝑆 holding numbers 1 through and the second is a vector of length 𝑄 holding the effort for each observation (site)
# holding the effort for each observation.
elist <- list(columns = 1:ncol(ydata),
               values = edata)

# Set up the formula for the model
form1 <- as.formula(~ elevation + summerSMOS + precip + landcover)

# Create rlist, which contains 𝑁 and 𝑟 for dimension reduction - HOW DO I SET THIS?
rlist   <- list(r = 8, N = 20)

# Set up priors

spLo <- "great blue heron"
sp <- length(spLo)
lo <- vector("list", sp)

# add names to the list
names(lo) <- paste0("landcoverwetlands_", spLo)

# add values to the list
lo[1] <- 0

spHi <- c("blue jay", "northern parula", "killdeer", "chipping sparrow", "ovenbird")
sp <- length(spHi)
hi <- vector("list", sp)

# add names to the list
names(hi) <- paste0("precip_", spHi)

# add values to the list
hi[1:length(hi)] <- 0

prior <- gjamPriorTemplate(formula = form1, 
                           xdata = xdata, ydata = ydata,
                           lo = lo, hi = hi)

# Set up prior vector with no priors???
lo <- list()
hi <- list()

prior <- gjamPriorTemplate(formula = form1, 
                           xdata = xdata, ydata = ydata,
                           lo = NULL, hi = NULL)

# Set up the model and run parameters
mlist <- list(ng=4000, burnin=2000, typeNames = 'DA', betaPrior = prior,
                  effort = elist, reductList = rlist)

# Run the model
outbirds <- gjam(form1, xdata = xdata, ydata = ydata, modelList = mlist)

gjamPlot(output = out, plotPars = list(PLOTALLY = T, SAVEPLOTS = T, outfolder = "./g"))

t <- out$parameters$sigMu %>% rowSums()

out$parameters$sigMu["hairywoodpecker", "hairywoodpecker"]
out$parameters$sigMu["redwingedblackbird", "redwingedblackbird"]
out$parameters$sigMu["mourningdove", "mourningdove"]
out$parameters$sigMu["northernflicker", "northernflicker"]
out$parameters$sigMu["wildturkey", "wildturkey"]
out$parameters$sigMu["easternwoodpewee", "easternwoodpewee"]
out$parameters$sigMu["canadagoose", "canadagoose"]
out$parameters$sigMu["brownheadedcowbird", "brownheadedcowbird"]


out$fit$yscore
out$fit$rmspeAll


# make plot for first species, Mourning Dove
md <- base +
  geom_point(aes(x = out$inputs$y[,1],
                 y = out$prediction$ypredMu[,1])) +
  labs(x = "observed", y = "predicted", title = paste0(colnames(out$inputs$y)[1], ", common"))

# make plot for a less common species, American Redstart
ar <- base +
  geom_point(aes(x = out$inputs$y[,8],
                 y = out$prediction$ypredMu[,8])) +
  labs(x = "observed", y = "predicted", title = paste0(colnames(out$inputs$y)[8], ", uncommon"))


ggarrange(md, ar, 
          ncol = 2)
