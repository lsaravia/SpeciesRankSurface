# Importance of manuscript

In this manuscript I describe a new way to analyze species abundance distributions: the species rank surface (SRS), using multifractal methods. These methods could be applied to a wide range of ecological communities where spatial data is available. I show utility of the method aplying it to spatially explicit community models. 

#

Different communities may have the same species-abundance, the simplest example is a neutral community without local dispersal, I mean one determined only by migration from metacommunity. When I add local dispersal the same SAD is observed and also with different dispersal Kernels the SAD is conserved but SAR and SRS should be different. 

- Metacommunity from periphyton and BCI (log-series)

 Modelos usar Con y Sin Metacommunity Con y sin saturacion Con y sin hierarchical

I use the two sample Kolmogorov-smirnov test 
http://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test

ks.test in R

# Notes

Chave, J., H. C. Muller-Landau, and S. A. Levin. 2002. Comparing Classical Community Models: Theoretical Consequences for Patterns of Diversity. American Naturalist 159:1-23.

termina diciendo que hacen falta otras medidas de las comunidades 



--------------

# Two approaches to constructing species-area curves

Condit et al. (1996) discussed two approaches to constructing S-N and species-area
curves. The usual approach for species-area curves (approach 1) is to count species in
nested series of expanding quadrats. The other approach (approach 2) is to add species
counts from one-size sampling quadrats randomly taken from across the whole
sampled plot. Approach 1 results in a curve showing the increase of the number of
species in a geographically expanding area, the appropriate way to compare different
sites. The species accumulation curve derived from approach 2 is not comparable across
sites because, as explained by Condit et al. (1996), that curve is dependent on the size
(area) of the plot sampled. Our computer programme for artificial sampling is in line
with approach 1.

Condit, R., S.P. Hubbell, J.V. Lafrankie, R. Sukumar, N. Manokaran, R.B. Foster & P.S. Ashton, 1996.
Species-area and species-individual relationships for tropical trees: a comparison of three 50-ha plots.
Journal of Ecology 84: 549–562.

-----------------

# discussion about why species have self-similar distributions

1. Šizling AL, Storch D (2004) Power-law species–area relationships and self-similar species distributions within finite areas. Ecol Lett 7: 60–68. doi:10.1046/j.1461-0248.2003.00549.x.

As Hubbell (2001) and others have pointed out, within the smallest scales the SAR certainly does not have a form of the power law, and the self-similarity does not hold down the level of individuals. We just show that when dealing with sufficiently large sampling plots (grid cells), the assumption of self-similarity is valid.

But why should species distributions be self-similar? One possibility is that this feature is imposed on species by the environment, i.e. that natural landscapes have self-similar properties. There is some evidence supporting this argument. Storch et al. (2002) showed that the spatial variability of biologically relevant parameters reveal spectral properties indicating self-similarity (so called 1/f spectra; see Halley 1996).
...
The power law in this case emerges because of the spatial aggregation of species which is not fully attributable to habitats. It is thus necessary to look for spatial population processes that generate self-similarity, i.e. to consider the dynamic nature of species assemblages (Adler & Lauenroth 2003).



