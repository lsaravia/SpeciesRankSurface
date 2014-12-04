# Importance of manuscript

In this manuscript I describe a new way to analyze species abundance distributions: the species rank surface (SRS), using multifractal methods. These methods could be applied to a wide range of ecological communities where spatial data is available. I show utility of the method aplying it to spatially explicit community models. 

#

Different communities may have the same species-abundance, the simplest example is a neutral community without local dispersal, I mean one determined only by migration from metacommunity. When I add local dispersal the same SAD is observed and also with different dispersal Kernels the SAD is conserved but SAR and SRS should be different. 


# Notes

Chave, J., H. C. Muller-Landau, and S. A. Levin. 2002. Comparing Classical Community Models: Theoretical Consequences for Patterns of Diversity. American Naturalist 159:1-23.

termina diciendo que hacen falta otras medidas de las comunidades 

# Interpretation of multifractal spectrum

We propose an interpretation of the multifractal spectrum in ecological terms as diversity patterns of subsets of species with a similar spatial distribution [Yakimov2008].

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


2.Instead of one process determining changes in species richness across a wide range of scales, different processes might determine plant biodiversity at different spatial scales. [Crawley2001]


3. The SAR describes the observed increase in the
number of species identified as the area sampled
increases. The SAR is one of the most studied patterns
in ecology (e.g. Connor & McCoy 1979; Rosenzweig
1995) [@White2010a]

4. The species-area relationship (SAR) is considered to be
one of a few generalities in ecology, yet a universal model of its shape
and slope has remained elusive. [@Sizling2011]
---


5. McGill BJ, Nekola JC (2010) Mechanisms in macroecology: AWOL or purloined letter? Towards a pragmatic view of mechanism. Oikos 119: 591–603. doi:10.1111/j.1600-0706.2009.17771.x.

Harte et al. have demonstrated how spatial fracti-
cality may give rise to the power law form of the species-area
relationship (Harte 2008), and that a fractal distribution of
individuals can lead to realistic species abundance distribu-
tions (Harte et al. 1999). These theories have subsequently
been extensively elaborated and tested (Green et al. 2003,
Green and Ostling 2003). Kunin and colleagues (Kunin
et al. 2000, Hartley et al. 2004) sought to use fractaility to
extrapolate from easily measured smaller scales up to impor-
tant but hard to measure larger scales. 

Similarly, several authors (Hawkins and Diniz-Filho 2002, Connolly et al.
2003, Kerr et al. 2006) have shown that random placement
alone is not enough to fully explain variations in diversity
across landscapes. However, random placement of species
with intra-specific clumping does a good job of parsimoni-
ously producing many key macroecological patterns (McGill
and Collins 2003, Harte et al. 2005, 2008).


------

1. Yakimov BN, Gelashvili DB, Solntsev LA, Iudin DI, Rozenberg GS (2014) Nonconcavity of mass exponents’ spectrum in multifractal analysis of community spatial structure: The problem and possible solutions. Ecol Complex 20: 11–22. doi:10.1016/j.ecocom.2014.07.003.

The slope of the SAR z is its main parameter, which
determines the growth rate of species richness with increasing
area. Theoretical limits for z are 0 (fixed species richness) and 1
(simple linear proportionality). Thus, in natural communities z is a
non-integer and so it corresponds to fractal dimension features.

Multifractal analysis is a generalisation of such a simple
approach, which also finds many applications in ecology. It is
applied to simulate neutral landscapes (Gamarra, 2005; Kirkpa-
trick and Weishampel, 2005) for time series analysis (reviewed by
Seuront, 2010), the description of microphytobenthos biomass
distribution (Seuront and Spilmont, 2002; Saravia et al., 2012a,b),
description of spatial structure of the tropical forest as a whole
(Sole and Manrubia, 1995a,b) and description of the individual
species distribution (Borda-de-Agua et al., 2007), as well as to
quantify the heterogeneity of species richness (Laurie and Perrier,
2010, 2011). Applications of multifractal analysis in ecology
mentioned above applied the standard version of this analysis,
which deals with the spatial distribution of a single measure (e.g.
biomass or stem density).

1. Hentschel HGE, Procaccia I (1983) The infinite number of generalized dimensions of fractals and strange attractors. Phys D 8: 435–444.

We note that not all fractals belong to this rescaling class. However, since all the generalized dimensions of interest are defined in the limit b ~ 0 (of. eqs. (1.1)-(1.4)), it is sufficient that there exists a length scale 1o below which such self similarity exists.

1. Beck C (1990) Upper and lower bounds on the Renyi dimensions and the uniformity of multifractals. Phys D Nonlinear Phenom 41: 67–78. doi:10.1016/0167-2789(90)90028-N.

These quantities describe the scaling behaviour of the regions in phase space where the measure is most concentrated (D(+infinit)) or most rarified (D(-infinit)).


1. Hurlbert AH (2004) Species–energy relationships and habitat complexity in bird communities. Ecol Lett 7: 714–720. doi:10.1111/j.1461-0248.2004.00630.x.

In this paper, my objective is twofold. First, I test the

Spatial autocorrelation is necessarily present in these geographic data, and some believe that it leads to inflated estimates of the number of degrees of freedom in significance tests. However, the statistically significant relationships presented in this study would remain significant at P < 0.01 if only 1% of the observations were statistically independent. Thus, the probability of an increased type I error rate due to spatial autocorrelation appears to be negligible for the relationships reported here.


1. Legendre P, Dale MRT, Fortin MJ, Gurevitch J, Hohn M, et al. (2002) The consequences of spatial structure for the design and analysis of ecological field surveys. Ecography (Cop) 25: 601–615.

Simulations have been performed to check the type I
error and estimate the power of the tests of significance
in the presence of different types of sampling designs
and spatial structures. Type I error occurs when the
null hypothesis is rejected while the data conform to
H0. A test of statistical significance is valid if the
rejection rate is not larger than the significance level a,
for any value of a, when the null hypothesis is true
(Edgington 1995). A test of significance should also be
able to reject the null hypothesis in most instances
when H0 is false. The ability to reject H0 in these
circumstances is referred to as the power of the test. In
the present study, power is the empirical rate of rejec-
tion of the null hypothesis when H0 is false by con-
struct. High power is a desirable property. When two
or more procedures are available (sampling designs or
tests of statistical significance), one should use the
procedure that has the highest power.


The rate of type I error is computed as the proportion
of rejection of the null hypothesis when the data con-
form to it. In our simulations, H0 is true if the environ-
mental (E) and response (R) variables are not linked by
the transfer parameter b. A test can be said to have
correct rate of type I error if, across the simulations,
the rejection rate is approximately equal to the signifi-
cance level a used to make the statistical decision.


1. Rosindell J, Cornell SJ (2007) Species–area relationships from a spatially explicit neutral model in an infinite landscape. Ecol Lett 10: 586–595. doi:10.1111/j.1461-0248.2007.01050.x.

The species–area exponent z relating to the power law
S μ Az is a key quantity for ecologists, but existing studies
give apparently inconsistent predictions of how z should
behave under the neutrality assumption. Zillio et al. (2005)
predict that z does not vary with m provided that m is small,
whereas the formulae of Durrett & Levin (1996) predict that z
should vary weakly with m. Hubbell (2001, p. 188) produces a
diagram showing how z should vary strongly with the
Ôfundamental biodiversity numberÕ h 1⁄4 2Jm for his patch
model. Chave et al. (2002) show that z varies with spatial scale
and can take realistic values, but only investigate SARs for
h 1⁄4 5. 
All these existing studies use biphasic SARs to identify
z. We hope to resolve these differences by identifying z as the
gradient in the central phase of our triphasic model. We show
how the value of z from our model can take any value between
0 and 1 and is a simple function of m alone, independent of the
other model parameters.


1. S. R. Supp, X. Xiao, S. K. M. Ernest, and E. P. White. 2012. An experimental test of the response of macroecological patterns to altered species interactions. Ecology 93: 2505-2511. doi:10.1890/12-0370.1

http://jabberwocky.weecology.org/2013/01/16/do-macroecological-patterns-respond-to-altered-species-interactions-research-summary/

Despite a long history of documenting macroecological patterns, an understanding of what determines pattern behavior, why patterns are so easily predicted by so many different models, and how we should go about addressing real ecological problems using a macroecological approach has still not been reached.

1. Rosindell J, Cornell SJ (2013) Universal scaling of species-abundance distributions across multiple scales. Oikos 122: 1101–1111. doi:10.1111/j.1600-0706.2012.20751.x.


The species abundance distribution (SAD) describes the distribution of commonness and rarity among species, but only at one spatial scale. In contrast the species area curve (SAC) spans many spatial scales, but only gives total species richness at each scale. It is increasingly recognised that SADs and SACs are not very informative of process because many different models can produce the same patterns

The combination of many SADs at different spatial scales is likely to be much more informative
(Chisholm and Lichstein 2009) because it indicates the commonness and rarity of species and spans many spatial scales.
