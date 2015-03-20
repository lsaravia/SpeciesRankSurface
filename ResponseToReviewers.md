
# Reviewers comments

Associate Editor Comments to Author:

Associate Editor
Comments to the Author:
After reading the manuscript and the comments of the two reviewers I think that this approach has high potential but that the work needs some considerable improvements as well. in addition to the comments by the reviewers I suggest to take care of the following:
#- In the beginning, I had some difficulties to follow which different metrics were actually tested against each other in the end and what different information they give us about the patterns and how we could interpret them in ecological terms. There is some text abut this in the beginning of the result section but this is much too late (and to me also an odd place to put this info). I would suggest that this is made clearer early on and maybe it would be good to add a table that lists the indices together with their ecological meaning.

- I moved the paragraph about interpretation to the methods section, and added more explanations, but as I didn't know in advance how Dq will be, I present most interpretations as a discussion. I added a table with the methods used.
 

# -It should be stated early on in which situations and for which questions this method can be used and it should be highlighted somewhere what assumptions need to be tested to apply it
- The assumptions are stated in methods and I am testing here the power of the test to compare communities, but other possible applications are suggested for Dq to characterize rare species. 

#- Finally, please make sure that color-blind (e.g. red-green) people can understand your figures.

- All figures have now a color-blind-safe palettes.


Reviewer(s)' Comments to Author:

Reviewer: 1

Comments to the Corresponding Author
The author proposes new macroecological metrics – species rank surface (SRS) and generalized dimensions spectrum DqSRS, which is the result of SRS multifractal analysis. New metrics are compared with the classical rank abundance distribution (RAD) and the spectrum of generalized dimensions DqSAD, which reflects the scaling of species abundance distribution (SAD). To demonstrate the possibilities of new metrics, author uses a simple toy model and a more realistic dynamic model of community. Author also has developed a procedure for statistical comparison of RAD, DqSRS and DqSAD for the latter model. This procedure is applied to detect differences between communities with varying degrees of neutrality as well as to demonstrate the differences in power between the metrics.
I highly appreciate the new method of analysis and consider it as having high potential. However, the current manuscript is not suitable for publication. The main drawbacks:

#1. The author offers no meaningful interpretation of dimensions DqSRS, restricting discussion to statements that dimension for large negative q reflect the distribution of rare species, and dimensions for large positive q reflect the distribution of the dominant species. Such lack of interpretation leads to that the new metric is used rather as a "black box", which can or cannot distinguish the communities.
- I added a some more interpretations, and about wath we know and what we don't know about Dq. This is why I added the simple toy models to understand more about how D_q is related to SAD and spatial patterns.   

#2. Toy examples are chosen poorly. I cannot imagine a community consisting of linear clusters of species ordered by their overall abundance.
- these are examples of very simple species gradients and don't intend to be realistic, and mimic in a primitive way an enviromental gradient.  

#3. Statistical procedure used for comparison of spectra is very unconventional and raises suspicions in its correctness. For comparison between a pair of two-dimensional curves author uses Anderson-Darling test, which is designed to compare the empirical distribution functions of two samples of one-dimensional random variables, samples must consist of independent observations. From the point of view of Anderson-Darling test the generalized dimensions spectrum Dq is a bimodal distribution with modes at D_-inf and D_+inf and relatively homogeneous distribution between them. Consequently, what is compared by Anderson-Darling test is just a position of the minimum and maximum of spectra. Obviously, this is not a procedure that is needed for spectra comparison. I did not encountered the procedure for the statistical comparison of the generalized dimensions spectra. If the author intends to develop such a procedure, its description and justification must be given with a stronger consideration.
- I added an explanation of why this test can be used with Dq

#4. The spectrum of generalized dimensions reflects the spatial scaling of some measure. The use of multifractal analysis implies self-similarity of the studied object, which manifests itself in the power-law form of dependence of partition function Zq on scale. Analysis of the relationship Zq (A) must precede any further analysis. If there is no power-law scaling, the object is not self-similar and its fractal analysis makes no sense. The power law looks like a linear function on a logarithmic scale. The author uses the coefficient of determination R^2 for the analysis of linearity. R^2 is not a measure of linearity, it reflects only the scatter of points around the regression line. The proposed rule of thumb was not justified. 
- I added an explanation about R^2

#The author shows figure with Zq(A) scaling only for toy models but not for dynamical models, whereas it is well-studied for neutral models that their species-area relations is triphasic with power law in a narrow range of scales. It is not clear whether there is power-law scaling for dynamical model in studied range of scales.
- Eliminate figure 4 with R2 add figure with fit for neutral model, and several figures in the appendix.

I believe that the current version of manuscript requires serious improvement and at this stage must be rejected and resubmitted.

More specific comments and suggestions:

#1. In the description of generalized dimensions spectra it is better to introduce two types (DqSRS and DqSAD) immediately to avoid confusion in formulas 2 and 5. It is better to use notation D1 in formula 3. Formula 5 must specify explicitly that Zq is an average over plots of area A.
- Changed in the ms.

#2. In the description of SRS details of the algorithm of assigning ranks to equally abundant species must be specified. I think that various options (random assignment of consecutive ranks or assignment of the same  average rank) will have a strong impact on the shape of DqSRS spectrum, especially for negative q. This aspect needs special consideration.
- I added the method used to calculate ranks.

#3. In the description of dynamic model (at least in Appendix) the following details must be specified: the boundary conditions, the type of lattice update (synchronous / asynchronous), the algorithm of assignment of competitive ranks in metacommunity. It is also necessary to raise the issue of the model dynamics stationarity, I propose to provide some figures of dynamics of several structural indicators (at least in the appendix).
- More details were specified in the appendix.

#4. D_0^SRS is not a fractal dimension of SRS (p. 9, l. 20). SRS is a rugged surface and its dimension is more than 2 but less than 3. What is really analyzed is some measure distributed over two-dimensional support. D_0^SRS is a dimension of this support.
- The D_0SRS is the fractal dimension of the species' spatial distribution, which in our model simulations is always 2 because the lattice is always filled. I changed it in the ms. 

#5. The degree of difference between Dq spectra is not a direct function of Delta rho. It is not correct to analyze relation of power to Delta rho. For example, the great drop in power of DqSAD in Fig. 8 is associated exclusively with the fact that the point Delta rho = 0.9 corresponds to the comparison of simulations rho = 1 vs rho = 0.9, which has very similar spectra.
- Changed to a matrix of plots.

#6. The power of methods based on spectra of indices is compared with power of methods based on single components of spectra. It is obvious that the latter will be less powerful. But those methods use different statistical procedures so its direct comparison is not fully correct.
- I added an explanation on why I calculate power for 1 dimension, the idea is to show that single components should not be considered for community comparisons due to its high type I error and low power. When you have different methods to do a test of hypothesis you allway prefer the one with more power, so it is commoon to compare power of different methods.

#7. Figures 1 and 5 lacks color scheme, making them difficult to read.
- I added the label and changed color palette

#8. Figure 3 has different scales for DqSAD и DqSRS, making it difficult to compare. As A = epsilon^2 I recommend to use A in both cases.
- I changed x axis to epsilon for both.  

#9. RAD is a strictly nonincreasing sequence, but many curves at figure 6 lack this property.
- I corrected the calculation of ranks

Reviewer: 2

Comments to the Corresponding Author
Any description of biodiversity that does not explicitly include species abundance distributions (SADs) is  likely to miss important relationships on how species interact with the environment. This is one of the rarest publications in the area of SAD modelling in which author(s) performed an elegant study investigating the scaling of SAD.
The author introduces and innovative framework,  “the species rank surfaces (SRS)”, as a new way to study the scaling of SAD.
Some comments for improvement:
#1)      In the Discussion the author should also highlight the limitations of their modelling process, in particular the fact that it not takes into account: a) asymmetry of the dispersal ability of species; b) the fact that trophic structure and species interactions in communities  may lead to more idiosyncratic processes of community assembly; c) the fact  that the truncated log-normal SAD was not tested.

- This paper is about a methodology to compare communities the models was choosen to generate a gradient with similar and different communities, but not to study the models or the limitation of the models, that would be another paper...

#2)      The following paper could be cited in line 7 of page 1: Triantis, K.A., Guilhaumon, F. & Whittaker, R.J. (2012). The island species–area relationship: biology and statistics. Journal of Biogeography, 39: 215-23.
#3)      In line 24 of page 1 it seems that one additional reference is missing afer Yakimove et al., 2008)
#4)      In page 14 line 23 “comparable” should be “comparable”
#5)      The author also should revise the English of the text.
#The bibliographic references are a Chaos! The same reference appears in the text in different formats (e.g. Borda-de-Água (2012) and the same correctly as Borda-de-Água et al. (2012). All the bibliographic references should be revised for consistency both in main text and List of Bibliography.
#6)      In the reference list for some strange reason all appear with "et al." in the end. I think that the author should learn to use the tools in Latex (see http://jabref.sourceforge.net/; http://en.wikibooks.org/wiki/LaTeX/Bibliography_Management
#For instance the reference below is a book Chapter:
#Borda-de-Água, L., Hubbell, S.P. & He, F.et al. (2007). Scaling biodiversity under neutrality. Scaling 3 biodiversity, 347–375
#and should be more correctly:
#Borda-de-Água, L., Hubbell, S. P., & He, F. (2007). Scaling biodiversity under neutrality. #In: Storch, D., Marquet, P. A., & Brown. J. H. (eds), Scaling biodiversity. Cambridge University Press, Cambridge, 347-375

- the references were corrected and the ms was read by a native  


 My best regards to the author,

Paulo A. V. Borges
Azorean Biodiversity Group
University of Azores
