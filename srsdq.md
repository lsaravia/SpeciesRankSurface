# Species rank distributions in space: a multifractal analysis

**Leonardo A. Saravia**, Ph.D.

Instituto de Ciencias Básicas

Universidad Nacional de General Sarmiento

J.M. Gutierrez 1159 (1613), Los Polvorines

Buenos Aires, Argentina.

<lsaravia@ungs.edu.ar>

\newpage

## Abstract

Species-area relationships (SAR) are one of the most studied patterns in ecology, due to its application in both theoretical and conservation issues. SARs has been generalized incorporating abundances using Hill's diversity indexes or Renyi entropies. The study of the scaling of Renyi entropies leads to the application of generalized dimensions ($D_q$) to study species diversity-area relationships. The calculation of $D_q$ is then based on the species-abundance distribution (SAD). Here I introduce a new way to study species diversity-area relationships using spatial version of RAD: the species-rank surface (SRS) that can be analyzed using $D_q$. To demonstrate the relationship of the D_q versions with species abundances and spatial patterns I perform simulations with regular and random patterns with different SAD. Finally I compare the power of the methods to detect different community patterns using a continuum of hierarchical and neutral spatially explicit models.

Keywords: multifractals, species-rank surface, species-area relationship, multi-species spatial pattern.


\newpage
	
## Introduction

The species-area relationship (SAR) is considered the one of oldest and best-documented patterns and one of a few generalization in ecology [@Sizling2011; @Crawley2001]. The SAR is often characterized through a triphasic curve, with a range of intermediate scales corresponding to power law relationship between the number of species and the area [@Preston1960; @Hubbell2001]. Although other quantitative forms could also be appropriate [@Tjorve2003; @White2010a] the power-law is more widely accepted [@Rosindell2007]. This relationship implies a self-similar or fractal structure of species distributions [@Harte1997]. 
SARs can be extended incorporating the abundances: instead of studying using only the scaling of the number of species, the Hill's generalized diversity index [@Hill1973] can be used. This index is an expression of the Renyi's generalized entropies [@Renyi1970]. The scaling of Renyi entropies is called generalized dimension and this is used to characterize multifractals. Thus extending the analysis of SAR have led us naturally to multifractals and this is the approach presented by Borda-de-Água [-@Borda-de-Agua2002]. But as SAD are represented by RAD we can use also the ranks to represent multi-species spatial distributions.

Multifractals and fractals are related techniques mainly used in physics to characterize scaling behavior of a system, the difference is that fractals look at the geometry of presence/absence patterns, while multifractals look at the arrangement of quantities like population densities or biomass [@Saravia2012a]. 

Metapopulation models in multifractal landscapes [@Gamarra2005] Analysis of natural multifractal landscapes [@Kirkpatrick2005] Multifractal search patterns in copepods [@Seuront2014]  

Multifractals were developed for nonlinear systems and they are suitable to capture multiplicative processes that can be observed in ecological systems. They characterize variability in a scale independent way within an experimental range.

Halley et al. [-@Halley2004] suggested that multifractals could be useful for ecological systems. This is because multifractals are associated with systems governed by random multiplicative process [@Stanley1988] and multiplicative process like the interaction of survival probabilities and compound growth are common in ecological systems [@McGill2003a]. Moreover, the widespread presence of multiplicative process is argued to produce the lognormal-like shape of most species-abundance distributions [@May2007a]. Also random processes with spatial correlations can generate multifractals [@Stanley1988], these kind of processes are part of neutral community models [@Houchmandzadeh2003] and are observed in natural communities [@Condit2000]. Thus there are a priori reasons to think that multifractals can be applied to much of the spatial ecological data. Indeed it has been applied to vegetal communities [@Scheuring1994], tropical forest [@Sole1995b, @Manrubia1996], microphytobentos and periphyton biomass patterns [@Seuront2002; @Saravia2012a], and as previously mentioned the characterization of species-area relationships [@Borda-de-Agua2002; @Yakimov2008; @Laurie2011]. But the application of multifractals is not widespread [@Seuront2009], one of the reasons is that is difficult to use the available software for quantitative multifractal analysis.

Rank-abundance distributions are a representation of species-abundance distributions (SAD) that are a classical description of a communities [@Mcgill2007]. These have been used to compare different communities and to compare models and data, but different mechanisms can produce nearly identical SADs [@Chave2002]. SADs are often presented using rank-abundance diagrams (RAD) where the log-abundance is plotted on the y-axis vs rank  on the x-axis [@Mcgill2007]. RADs are equivalent to cumulative  distributions [@Newman2005] and thus are a robust way to visualize the SAD without losing information[@Li2005]. Here I propose to extend the analysis of SAD attaching the rank of each species to its spatial distribution, in this way the multivariate spatial distribution of species is summarized into an univariate 2D distribution. I called this spatial distributions the species-rank surface (SRS), that can be analyzed and compared using multifractals.

Multifractals require that the object under study should be statistically self-similar, that means that a power-law could be fitted to data in a range of scales. But that does not mean that the power-law must be the best possible model. We can use as a technique to analyze the data without claiming that it is an exact multifractal. Another point of view is not to impose self-similarity and analyze multifractality as a local property [@Bez2010], both approaches are valid but should some care is needed to make comparisons. 

One of the advantages of multifractals is that they require fewer conditions on data than more classical statistics like autocorrelation and variograms. These usually require isotropy and stationarity [@Fortin2006] but multifractals can be used with anisotropic data [@Harte2001]  and are inherently non-stationary [@Laurie2010; @Bez2010]. Anisotropy and nonstationarity are often seen in spatial ecological distributions [@Plotkin2002].

Here I show that multifractals can be applied to different data sets and models to demonstrate its utility to characterize ecological data and propose the spatial analysis of rank-abundance distributions using the species-rank surface. 

Is the calculation of a spectrum of dimensions needed to differentiate communities or a single dimension with less and simpler calculations should be enough?

# Methods

## Multifractal analysis

Extensive reviews of multifractal methods applied to ecology are available [@Seuront2009] and some good introductions are also published [@Scheuring1994; @Borda-de-Agua2007]. Thus I will only give a brief overview. Multifractals analyze the scaling properties of quantities distributed in a space that we assume to be two dimensional (a plane). Multifractasl can mathematically represented using different  ways [@Harte2001], from these the most close to ecology are the generalized dimensions $D_q$ [@Grassberger1983], also called Renyi dimensions [@Renyi1970]. D_q$ has been used to characterize the probabilistic structure of attractors derived from dynamical systems [@Hentschel1983]. 

To estimate generalized dimensions I used the method of moments based on box-counting [@Evertsz1992]. The spatial distribution of quantities $\mu$ is covered with a grid, which divided it into $N(\epsilon)$ boxes of side  $\epsilon$, allowing us to calculate the quantity $\mu_i(\epsilon)$ in each of them. Then the so called partition function is computed as: 

 (@e1) $Z_{q}(\epsilon)={\displaystyle \sum_{i}^{N(\epsilon)}} \left(\mu_{i}(\epsilon)\right)^{q}$

Where $q$ is called moment order. The operation is performed for different values of $\epsilon$ and $q$, within a predetermined range. The generalized dimension is calculated as:

 (@e2) $D_{q}=\cfrac{1}{q-1}\underset{\epsilon\to0}{\lim}\cfrac{\log\left(Z_{q}(\epsilon)\right)}{\log\epsilon}$

When $q = 1$, the denominator of the first term in $D_q$ is undefined, so it must be replaced by the following expression:

 (@e3) $D_{q}=\underset{\epsilon\to0}{\lim}\cfrac{{\displaystyle \sum_{i}^{N(\epsilon)}}\mu_{i}(\epsilon)\log\left(\mu_{i}(\epsilon)\right)}{\log\epsilon}$

In practical cases as the limit can not be assessed, the dimensions are estimated as the slope of the $log(Z_q)$ versus $log(\epsilon)$ in equation (@e1) replacing by the numerator in equation (@e3). This is done for different $q$, provided that it is a real number which yields a graphs of $D_q$ in terms of $q$, called the spectrum of generalized dimensions. 

To be an approximate multifractal the relationship $log(Z_q)$ versus $log(\epsilon)$ should be well described by a linear relationship, but also a linear relationship with superimposed oscillations is accepted [@Borda-de-Agua2007]. A range of $q$ and $\epsilon$ must be established, then $D_q$ is estimated using linear regressions. Note that $D_q$ are defined in the limit ${\epsilon\to0}$ (equations @e2 and @e3) thus it is sufficient that there exists a scale below which such linear relationship applies to use the method [@Hentschel1983]. 

To analyze species-abundance-area relationships with multifractals as Borda-de-Água [-@Borda-de-Agua2002], the boxes are replaced by species. Thus at each spatial scale $\epsilon$ each species holds the quantity of interest: its own abundance. Then the partition function is defined as a sum over the species present $S(A)$ in an area $A$ and the side of the box $\epsilon$ is replaced by the area:

(@e4) $Z_{q}(A)={\displaystyle \sum_{i}^{S(A)}} \left(\mu_{i}(A)\right)^{q}$

Where $\mu_{i}(A)$ is the abundance of species $i$ in an area $A$. And $D_q$ is defined as:

(@e5) $D_{q}=\cfrac{1}{1-q}\underset{A\to\infty}{\lim}\cfrac{\log\left(Z_{q}(A)\right)}{\log A}$

I call $D_q$ calculated from species abundances $D_q^{SAD}$ and $D_q$ calculated from SRS $D_q^{SRS}$, when I mention $D_q$ I refer to both. As I already mentioned $D_q^{SAD}$ characterize the scaling of the Hill's generalized diversity index [@Hill1973]. Thus when the moment order $q=0$, $D_q^{SAD}$ represent the exponent of the SAR power-law scaling, when $q=1$, $D_q^{SAD}$ represent the scaling of Shannon diversity index and when q=2, $D_q^{SAD}$ becomes the scaling of Simpson's index. 

Theoretically $D_q$ must be a non-increasing function of $q$ [@Hentschel1983], which means that if $q_1 \geq q_2$ then $D_{q1} \leq D_{q2}$. Some studies showed small violations to this property for $D_q^{SAD}$ [@Borda-de-Agua2002; @Zhang2006]. These violations are be related to the way that $D_q^{SAD}$ is defined: the summation of equation @e4 is over species and the summation of equation @e1 is over boxes, this changes the way in which the mathematical limits are taken and also the way of computation of $D_q^{SAD}$. A partial solution has been proposed [@Yakimov2014] but the anomalies observed may be related to the mathematical assumptions needed for $D_q$ to be non-increasing, a new mathematical proof should be developed for $D_q^{SAD}$. Thus as long as the linear relationship is reasonable I take $D_q^{SAD}$ as a useful technique of analysis. 

I proposed a new way to analyze species-abundance-area using multifractals that fits more closely to the original definitions of equations @e1 - @e3: the species-rank surface [@Saravia2014]. To construct the species-rank surface (SRS) the spatial distribution of species have to be transformed assigning to each species position its rank. First I use the species abundances, at the whole plot level, to calculate the rank ordering the species from biggest to lowest and assigning a number starting with one. Then the rank is assigned to the spatial position of the individuals of each species forming a surface. This landscape have valleys formed by the most abundant species and peaks determined by the most rare species. Finally the standard multifractal analysis is applied. If sampling was made in quadrats without taking the spatial position of individuals, the sum of the ranks of the species in the smallest quadrats could be used to form the SRS. 

I use the coefficient of determination ($R^2$) as a descriptive measure of goodness of fit [@Borda-de-Agua2002]. The source code in C++ to perform multifractal analysis is available at <https://github.com/lsaravia/mfsba>. 

### $D_q$ relation with different spatial patterns and different SAD

I simulated species spatial patterns with different SAD's to demonstrate how $D_q$ is related to them. First I used a uniform SAD, in this case all species have the same densities and the same number of individuals. To add a degree of stochasticity I take the number of individuals of each species from a Poisson distribution with the same mean. I distributed them in bands over a spatial grid so they form a regular spatial pattern. Each position of the grid is occupied by one individual and I choose the number of species to exactly divide the side of the grid so all species are strips with approximately the same width (Figure 1).
I used square grids with sides of 256 and 512 sites which contain 65536 and 262144 individuals respectively, and 8, 64 and 256 species. I calculated $D_q$ for the regular pattern, then I randomized the positions of species to compare $D_q$ obtained with these two extreme cases.
The second SAD I used is a Logseries [@Fisher1943] with the same number of species and the same sides as previously. I used the R package untb [@Hankin2007] to calculate the density for each species, basically they use Poisson distribution with the expected Logseries abundances as means. I then build the regular pattern with strips of species, but as species have different abundances the widths for each species are different (Figure 1). Then I estimated $D_q$ for the regular and randomized patterns. I simulated 10 spatial patterns for each case and calculated the mean and standard deviation of $D_q$. 

### Spatially explicit model

To simulate more realistic patterns of species-abundance-area relationships I used a stochastic spatially explicit model. I developed a stochastic cellular automata [@Molofsky2004] that can switch between a neutral or hierarchical competition representing a continuum between niche and neutral communities [@Gravel2006]. In a neutral model individuals do not interact, and have all the same mortality, colonization rates and dispersal distances. In spite of these gross simplifications neutral models are capable of predict several real community patterns [@Rosindell2011]. At the other end of the continuum are niche communities represented by a hierarchical competition model [@Tilman1994]. In this case species have differences that imply a competitive hierarchy in which some species are always better than others, then they produce competitive exclusion [@Chave2002]. Basically I added to the neutral model a probability of replacement $\rho$. When $\rho=1$ more competitive species always replace less competitive and the model behaves as pure hierarchical. When $\rho=0$, there is no replacement of species and the model is completely neutral. A more thorough description of the model is given in appendix A and the C++ source code of the model is available at <https://github/lasaravia/neutral> and figshare <http://dx.doi.org/10.6084/m9.figshare.969692>. 

Following a classical neutral scheme the model has a metacommunity: a regional collection of communities, from which migration occurs at a rate $m$. Species can also disperse locally and I assume an exponential dispersal kernel with average dispersal distance $d$. Other model parameters are the mortality rate $\mu$, the number of species in the metacommunity and also the size of the community, this last is represented as the *side* of the grid used in the simulations. I use a logseries SAD for the metacommunity, defined by the maximum number of individuals (*side x side*) and the number of species [@Fisher1943].  

The values of the parameters were in the range of estimated for BCI from the existing literature [@Anand2010; @Condit2002; @Etienne2007]. I performed 50 simulations for each combination of parameters given in Table 1. To compute the power I made comparisons of communities with different levels of $\rho$, representing more neutral or hierarchical communities, the other parameters were kept constant. I also made comparisons between repetitions with the same $\rho$ to calculate the type I error. 

+------+-------------+-------+-----+--------+--------+
| Side | No. Species | $\mu$ | $d$ |  $m$   | $\rho$ |
+======+=============+=======+=====+========+========+
|  256 |          11 |   0.2 |  25 |  0.001 |      1 |
+------+-------------+-------+-----+--------+--------+
|  512 |          86 |       |     |        |    0.1 |
+------+-------------+-------+-----+--------+--------+
|      |         341 |       |     |        |   0.01 |
+------+-------------+-------+-----+--------+--------+
|      |             |       |     |        |  0.001 |
+------+-------------+-------+-----+--------+--------+
|      |             |       |     |        |      0 |
+------+-------------+-------+-----+--------+--------+

Table: Parameters values used in the simulations of the neutral-hierarchical model 



### Statistical comparison of methods

I analyzed the performance of two kinds of methods to differentiate communities. Two methods are based on a spectrum of dimensions: $D_q^{SAD}$ and $D_q^{SRS}$. The other kind are based on a single dimension or power exponent: the SAR exponent and the information dimension.  For the first kind I used a non-parametric test related to the Kolmogorov-Smirnov (KS) test: the Anderson-Darling (AD) test [@Feigelson2012]. This test measure the differences between the empirical distribution functions (EDF) of two datasets as a weighted sum of square deviations between the EDFs. In extensive simulations the AD test has proven more sensitive than the KS test [@Stephens1974]. I use the package kSamples [@Scholz2012] in the R statistical language. 

The SAR exponent is part of the $D_q^{SAD}$ spectra when q=0 [@Borda-de-Agua2002]. An equivalent single number measure from $D_q^{SRS}$ is the information dimension [@Ricotta2000;@Chappell2001], that is the $D_q^{SRS}$ when q=1. I calculate the power of this two indexes with a T-test using the standard deviation (SD) obtained from the box-counting procedure to estimate the multifractal spectrum. These SD are obtained with autocorrelated data because small squares are nested within big squares  (See Multifractal Analysis). The consequence is that the SD may be underestimated, but the slopes estimates are still unbiased [@Kutner2005]. This should result in an increased type I error rate and also in a spurious increase in power. 

### Calculation of power and type I error

I estimated the power for simulated communities with different degree of hierarchical structure given by the parameter $\rho$ of the model. The power of a test is the ability to reject the null hypothesis ($H_0$) when it is false. The significance level to reject $H_0$ was set a priori at $\alpha$ = 0.05 in all cases, and the rejection rate of each test was calculated as the proportion of P values that were less than or equal to $\alpha$. To estimate power I used independent simulations of communities (50 repetitions) with the same parameters except $\rho$, and compute the proportion of rejection.

The type I error is the probability to reject $H_0$ when is true (false positive). In our simulations, $H_0$ is true if two simulated communities have the same $\rho$ (and also are equal in the other parameters). To estimate type I error I compare independent simulations of communities with the same set of parameters (50 repetitions) and compute the proportion of rejection.  


# Results

The interpretation of $D_q$ is related to the interpretation of SAR power law exponent, when its bigger it means that when we change the scale of observation the change in the number of species is greater. $D_q$ express the change of the quantities under study when scale change but modulated by $q$. When $q$ is positive the terms of the sums (equations 2 & 5) with more abundant species have more weight and become even more important when $q$ is greater. When $q$ is negative we have the opposite pattern, less abundant species have more weight, so $D_q$ reflects the change of rare species. When $q$ is greater in absolute value of $D_q$ is driven by more and more extreme values thus $D_q$ will have a bigger variance.  Here I present most figures with a range of $q$ from -24 to 24 but for statistical comparisons I use a smaller range from -10 to 10.

I calculated two versions of $D_q$: a) the original definition due to Borda-de-Água [-@Borda-de-Agua2002] where $D_q$ measure the change in SAD as we change scale. b) the one based in SRS, where we measure the change of the abundances and the spatial distributions of species as we change scale. 

- HACER GRAFICO con 4 curvas: logseries {randomized, uniform} regular {randomized, uniform} y mostrar la variacion con el numero de especies 

For the uniform pattern $D_q$s are greater for the negative part, as the pattern is uniform a symmetric D_q around 2 should be expected. The difference is due to the fact that the negative part is analyzing numbers close to 0 and the logarithm enhances the differences between small numbers [@Laurie2011] . 


## Ecological communities



## Neutral Model



# Discussion

Multifractal spectrum $D_q$ can be used to describe spatial patterns in all the cases I analyzed. Multifractals patterns could be produced by the existence of multiplicative interaction between species but this is not the only possibility, spatially correlated random processes like dispersal and growth to adjacent areas would also produce multifractals [@Stanley1988]. Because plant and animal species are generally aggregated in space is very likely that multifractal analysis can be used in a wide range of cases.
The occurrence of multifractal patterns with the same $D_q$ in several places does not prove that the same mechanism is acting, but may provide stronger evidence for similar mechanisms than if only one fractal dimension is estimated. Usually this dimension corresponds to $D_0$, that characterizes the geometric complexity of the spatial distribution, the other dimensions in the multifractal spectrum characterize the non-uniformity of the distribution [@Harte2001]. This adds much information to the characterization of spatial patterns that is lost if only one dimension is calculated.

 


# Acknowledgments

I am grateful to the National University of General Sarmiento for financial support and to Fernando R. Momo for our great conversations about ecological theory. I also wish to thank to Diana Marco and Luis Borda-de-Água for useful suggestions to the manuscript, and Graeme Ruxton for his english revision. 

# References  