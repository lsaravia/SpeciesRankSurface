# Species rank distributions in space: a multifractal analysis

**Leonardo A. Saravia**, Ph.D.

Instituto de Ciencias Básicas

Universidad Nacional de General Sarmiento

J.M. Gutierrez 1159 (1613), Los Polvorines

Buenos Aires, Argentina.

<lsaravia@ungs.edu.ar>

\newpage

## Abstract

Species-area relationships (SAR) are one of the most studied patterns in ecology, due to its application in both theoretical and conservation issues. SARs has been generalized incorporating abundances using Hill's diversity index and multifractals to study species diversity-area relationships. Multifractals were developed for nonlinear systems and they are suitable to capture multiplicative processes that can be observed in ecological systems. They characterize variability in a scale independent way within an experimental range. Another well studied pattern in ecology is the species-abundance distribution (SAD) that characterize rarity and commonness in ecological communities. SADs can be represented and compared in a robust way using rank-abundance diagrams (RAD). Here I introduce a new way to analyze species diversity-area relationships using spatial version of RAD: the species-rank surface (SRS) that can be analyzed using multifractals. I use the spectrum of generalized dimensions $D_q$ for multifractal representation. To demonstrate the relationship of SRS with species spatial patterns I perform simulations with regular and random patterns and also using a neutral model. Finally I compare the power of the methods to detect different patterns.

Keywords: multifractals, species-rank surface, species-area relationship, multi-species spatial pattern.


\newpage
	
## Introduction

The species-area relationship (SAR) is considered the one of oldest and best-documented patterns and one of a few generalization in ecology [@Sizling2011; @Crawley2001]. The SAR is often characterized through a triphasic curve, with a range of intermediate scales corresponding to power law relationship between the number of species and the area [@Preston1960; @Hubbell2001]. Although other quantitative forms could also be appropriate [@Tjorve2003; @White2010a] the power-law is more widely accepted [@Rosindell2007]. This relationship implies a self-similar or fractal structure of species distributions [@Harte1997]. 
SARs can be extended incorporating the abundances: instead of studying using only the scaling of the number of species, the Hill's generalized diversity index [@Hill1973] can be used. This index is an expression of the Renyi's generalized entropies [@Renyi1970]. The scaling of Renyi entropies is called generalized dimension and this is used to characterize multifractals. Thus extending the analysis of SAR have led us naturally to multifractals and this is the approach presented by Borda-de-Água [-@Borda-de-Agua2002]. But as SAD are represented by RAD we can use also the ranks to represent multi-species spatial distributions.

Multifractals and fractals are related techniques mainly used in physics to characterize scaling behavior of a system, the difference is that fractals look at the geometry of presence/absence patterns, while multifractals look at the arrangement of quantities like population densities or biomass [@Saravia2012a]. 

Metapopulation models in multifractal landscapes [@Gamarra2005] Analysis of natural multifractal landscapes [@Kirkpatrick2005] Multifractal search patterns in copepods [@Seuront2014]  


Halley et al. [-@Halley2004] suggested that multifractals could be useful for ecological systems. This is because multifractals are associated with systems governed by random multiplicative process [@Stanley1988] and multiplicative process like the interaction of survival probabilities and compound growth are common in ecological systems [@McGill2003a]. Moreover, the widespread presence of multiplicative process is argued to produce the lognormal-like shape of most species-abundance distributions [@May2007a]. Also random processes with spatial correlations can generate multifractals [@Stanley1988], these kind of processes are part of neutral community models [@Houchmandzadeh2003] and are observed in natural communities [@Condit2000]. Thus there are a priori reasons to think that multifractals can be applied to much of the spatial ecological data. Indeed it has been applied to vegetal communities [@Scheuring1994], tropical forest [@Sole1995b, @Manrubia1996], microphytobentos and periphyton biomass patterns [@Seuront2002; @Saravia2012a], and as previously mentioned the characterization of species-area relationships [@Borda-de-Agua2002; @Yakimov2008; @Laurie2011]. But the application of multifractals is not widespread [@Seuront2009], one of the reasons is that is difficult to use the available software for quantitative multifractal analysis.


Rank-abundance distributions are a representation of species-abundance distributions (SAD) that are a classical description of a communities [@Mcgill2007]. These have been used to compare different communities and to compare models and data, but different mechanisms can produce nearly identical SADs [@Chave2002]. SADs are often presented using rank-abundance diagrams (RAD) where the log-abundance is plotted on the y-axis vs rank  on the x-axis [@Mcgill2007]. RADs are equivalent to cumulative  distributions [@Newman2005] and thus are a robust way to visualize the SAD without losing information[@Li2005]. Here I propose to extend the analysis of SAD attaching the rank of each species to its spatial distribution, in this way the multivariate spatial distribution of species is summarized into an univariate 2D distribution. I called this spatial distributions the species-rank surface (SRS), that can be analyzed and compared using multifractals.

Multifractals require that the object under study should be statistically self-similar, that means that a power-law could be fitted to data in a range of scales. But that does not mean that the power-law must be the best possible model. We can use as a technique to analyze the data without claiming that it is an exact multifractal. Another point of view is not to impose self-similarity and analyze multifractality as a local property [@Bez2010], both approaches are valid but should some care is needed to make comparisons. 

One of the advantages of multifractals is that they require fewer conditions on data than more classical statistics like autocorrelation and variograms. These usually require isotropy and stationarity [@Fortin2006] but multifractals can be used with anisotropic data [@Harte2001]  and are inherently non-stationary [@Laurie2010; @Bez2010]. Anisotropy and nonstationarity are often seen in spatial ecological distributions [@Plotkin2002].

Here I show that multifractals can be applied to different data sets and models to demonstrate its utility to characterize ecological data and propose the spatial analysis of rank-abundance distributions using the species-rank surface. 

# Methods

## Multifractal analysis

Extensive reviews of multifractal methods applied to ecology are available [@Seuront2009] and some good introductions are also published [@Scheuring1994; @Borda-de-Agua2007]. Thus I will only give a brief overview. Multifractals analyze the scaling properties of quantities distributed in a space that we assume to be two dimensional (a plane). Multifractasl can mathematically represented using different  ways [@Harte2001], from these the most close to ecology are the generalized dimensions $D_q$ [@Grassberger1983], also called Renyi dimensions [@Renyi1970]. $D_q$ has been used to characterize the probabilistic structure of attractors derived from dynamical systems [@Hentschel1983]. 

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

I call $D_q$ calculated from species abundances $D_q^{SAD}$ and $D_q$ calculated from SRS $D_q^{SRS}$, when I mention $D_q$ I refer to both.

Theoretically $D_q$ must be a non-increasing function of $q$ [@Hentschel1983], but several studies showed small violations to this property for $D_q^{SAD}$ [@Borda-de-Agua2002; @Zhang2006]. These violations are be related to the way that $D_q^{SAD}$ is defined: the summation of equation @e4 is over species and the summation of equation @e1 is over boxes, this changes the way in which the mathematical limits are taken and also the way of computation of $D_q^{SAD}$. A partial solution has been proposed [@Yakimov2014] but the anomalies observed may be related to the mathematical assumptions needed for $D_q$ to be non-increasing, a new mathematical proof should be developed for $D_q^{SAD}$. Thus as long as the linear relationship is reasonable I take $D_q^{SAD}$ as a useful technique of analysis. 

I proposed a new way to analyze species-abundance-area using multifractals that fits more closely to the original definitions of equations @e1 - @e3: the species-rank surface [@Saravia2014]. To construct the species-rank surface (SRS) the spatial distribution of species have to be transformed assigning to each species position its rank. First I use the species abundances, at the whole plot level, to calculate the rank ordering the species from biggest to lowest and assigning a number starting with one. Then the rank is assigned to the spatial position of the individuals of each species forming a surface. This landscape have valleys formed by the most abundant species and peaks determined by the most rare species. Finally the standard multifractal analysis is applied. If sampling was made in quadrats without taking the spatial position of individuals, the sum of the ranks of the species in the smallest quadrats could be used to form the SRS. 

I use the coefficient of determination ($R^2$) as a descriptive measure of goodness of fit [@Borda-de-Agua2002]. The source code in C++ to perform multifractal analysis is available at <https://github.com/lsaravia/mfsba>. 

### $D_q$ relation with different spatial patterns and different SAD

I simulated different spatial patterns and SAD's to demonstrate how $D_q$ is related to them. The first SAD I used is a uniform SAD, this has the same number of individuals for all species, I distributed them in bands over a spatial grid so they form a regular spatial pattern. Each position of the grid is occupied by one individual and I choose the number of species to exactly divide the side of the grid so all species are strips of the same width.
I used square grids with sides of 256 and 512 sites which contain 65536 and 262144 individuals respectively, and 8,64, and 256 species. I calculated both $D_q$ then I randomized the positions of species to compare $D_q$ obtained with these two extreme cases.
The second SAD I used is a Logseries [@Fisher1943] with the same number of species as previously and the same sides. I also fill the grids with strips of species to build a regular pattern but as species have different abundances the widths for each species are different. Then I calculate $D_q$ for the regular and randomized patterns. I simulated 10 spatial patterns for each case and calculate the mean and standard deviation of $D_q$. 

### Spatially explicit model

To simulate more realistic patterns of species-abundance-area relationships I use a stochastic spatially explicit model. I developed a stochastic cellular automata [@Molofsky2004] that can switch between a neutral or hierarchical competition representing a continuum between niche and neutral communities [@Gravel2006]. In a neutral model individuals do not interact, and have all the same mortality, colonization rates and dispersal distances. In spite of these gross simplifications neutral models are capable of predict several real community patterns [@Rosindell2011]. At the other end of the continuum are niche communities represented by a hierarchical competition model [@Tilman1994]. In this case species have differences that imply a competitive hierarchy in which some species are always better than others, then they produce competitive exclusion [@Chave2002]. Basically I added to the neutral model a probability of replacement $\rho$. When $\rho=1$ more competitive species always replace less competitive and the model behaves as pure hierarchical. When $\rho=0$, there is no replacement of species and the model is completely neutral. A more thorough description of the model is given in appendix A and the C++ source code of the model is available at <https://github/lasaravia/neutral> and figshare <http://dx.doi.org/10.6084/m9.figshare.969692>. 

The model assume a metacommunity: a regionall collection of communities, from which migration occurs at a rate $m$. Species can also disperse locally and I assume an exponential dispersal kernel with average dispersal distance $d$. Other model parameters are the mortality rate $\mu$, the number of species in the metacommunity and also the size of the community, this last is represented as the *side* of the grid used in the simulations. I use a logseries SAD for the metacommunity, defined by the maximun number of individuals (*side x side*) and the number of species [@Fisher1943].  

The values of the parameters were in the range of estimated for BCI from the existing literature [@Anand2010; @Condit2002; @Etienne2007]. I performed 50 simulations for each combination of parameters given in Table 1. To compute the power I made comparisons of communities with different levels of $\rho$, representing more neutral or hierarchichal communities, the other parameters were kept constant. I also made comparisons between repetitions with the same $\rho$ to calculate the type I error. 

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

I analyzed the comparative performance of three methods to differentiate comunities so I needed a statistical test that could be applied to the curves produced by all methods. Thus I used a non-parametric test related to the Kolmogorov-Smirnov (KS) test: the Anderson-Darling (AD) test [@Feigelson2012]. This test measure the diferences between the empirical distribution functions (EDF) of two datasets as a weighted sum of square deviations between the EDFs. In extensive simulations the AD test has proven more sensitive to the KS test [@Stephens1974]. I use the package kSamples [@Scholz2012] in the R statistical languaje. 

The SAR exponent is part of the $D_q^{SAD}$ espectra when q=0 [@Borda-de-Agua2002]. An equivalent single number measure from $D_q^{SRS}$ is the information dimension [@Ricotta2000;@Chappell2001], that is the $D_q^{SRS}$ when q=1. I calculate the power of this two indices with a T-test using the standard deviation (SD) obtained from the box-counting procedure to estimate the multifractal spectrum. These SD are obtained with autocorrelated data because small squares are nested within big squares  (See Multifractal Analysis). The consequence is that the SD may be underestimated, but the slopes estimates are still unbiased [@Kutner2005]. Thus should result in an increased type I error rate and also an increase in rejection rate. 

I used another method to statistical compare $D_q^{SAD}$ and $D_q^{SRS}$  a permutation test which gives a global p-level (Function compareGrowthCurves from R package statmod [@Smyth2011]). This test calculates the average T statistic over all points of the curve, then the distribution of this mean T is obtained via permutations. 

### Calculation of power and type I error

Using the neutral model simulations I estimate the power and type I error for communities with different degree of hierarchical structure. To estimate the power of the test of significance I use simulations generated with different $\rho$, the power of a test is the avility to reject $H_0$ when $H_0$ is false. The significance level to reject the null hypothesis ($H_0$) was set a priori at $\alpha$ = 0.05 in all cases, and the rejection rate of each test was calculated as the proportion of P values that were less than or equal to $\alpha$. 

The type I error is the probability to reject $H_0$ when is true (false positive). In our simulations, $H_0$ is true if two simulated communities have the same $\rho$ (and also are equal in the other parameters). To estimate type I error I compare independent simulations of communities with the same parameters (repetitions) and compute the proportion of rejection. 


# Results

The two indexes derived from the multifractal spectrum capture the change in $D_q$. Information dimension $D_1$ diminish and $\Delta D_i$ increases when $p_1$ raises (Table 1). 
I also show in Table 1 the SD calculated using bootstrap and the SD calculated from the regressions used to estimate $D_q$. The last ones are estimated from the spatial distribution, and the previous needs a number of repetitions of the process. As expected the SD from regressions is very small compared with the bootstrapped SD. 

These indexes capture the gross behavior of the $D_q$ but the shape of the spectrum can give us much more information. In the case of the spectra generated by p-model the left and the right part of the spectrum are more like inverted specular images (Figure 1c) but this is not always the case. Different shapes of the $D_q$ can give us clues about different processes that should be acting. 


I compared two similar p-models with parameter $p_1$=0.5 and 0.51 (Figure 2). The images and the $D_q$ spectra were very similar, the bootstrapped CI overlaps mainly in the region with higher $q$ absolute value. But the overlaps are very small suggesting significant differences at 0.05 level [@Cumming2009]. To confirm these results I use the mentioned permutation test [@Smyth2011]. This test calculates the average T statistic over all points of the curve, then the distribution of this mean T is obtained via permutations. Because the differences in the whole $D_q$ spectra change sign from $q$ negative to positive the differences in the T statistic cancel out, thus I tested separately the negative and positive parts of $D_q$ (Table 2). The results coincide with the visual inspection of the bootstrapped CI, for both $q$ positive and $q$ negative the $D_q$ are different. I also compared the $p_1$=0.5 images with another set generated with the same $p_1$ and repeated the procedure for $p_1$=0.3 and 0.7 (Table 2). For $p_1$=0.5 the results confirm the visual inspection of the CI there are differences at both sides of the $D_q$ spectrum and there is no significant differences when the comparison is for realizations of the model with the same $p_1$. With $p_1$=0.3 the results are similar but for $p_1$=0.7 the $D_q$ positive shows no differences for distinct $p_1$. This could be explained because at this level of $p_1$ there is a greater spatial variability and also to a bigger variability between realizations of p-model. The regions of high values are dominated by extreme events that are reflected in the positive $D_q$ which make more difficult the detection of differences. The regions with lower values are described by negative $D_q$ and they do not have such extreme events thus are highly significant.  


## Ecological communities

The $D_q$ calculated from spatial biomass distribution also does not show differences in time. The range of $R^2$ starts at 0.76 so the fit is not so good (Figure S2). Anyway deviations are expected even in true multifractals if the chosen box sizes $\epsilon$ does not match the multiplicative process [@Chhabra1989]. Here I don't try to demonstrate that the underlying process is multifractal but if the fits were poor the multifractal analysis (MFA) will not be valid. Here more than 50% of $R^2$ are greater than 0.95 thus they are enough to consider MFA valid.

I present the results of MFA for three aquaria (named D E & I) with early and late periphyton development. For periphyton communities with 3 days of development we observe that $D_q$ for $q<0$ are nearly flat and for $q>0$ the curves extends to lower values, this means that these images are composed with a uniform surface of relatively low values (the negative part) and some isolated spots of relatively high biomass (Figure 4b & d). For late periphyton communities (80 days of development) the pattern is inverted $D_q$ is flatter at $q>0$, that means a spatial distribution more uniform at higher biomass. For $q<0$ the curve $D_q$ raises and the patches of relatively low biomass became important. There are differences between two eighty days images, $D_q$ for E and I aquaria are higher than D aquaria, this means that the low biomass spots are more profound in extension and intensity for E and I (Figure 4 b & c).
From previous studies I already know that $D_q$ from periphyton have a very good fit, the same happens here for this new dataset with $R^2$ greater than 0.98 (Figure S2). 


## Neutral Model

I chose to simulate the model with a set of parameters that result in a similar richness found in BCI. The output of the model is the spatial distribution of species thus I use SRS to calculate $D_q$. The fit of all $D_q$ were very good, the range of $R^2$ was 0.947-0.999.

I calculated SD using bootstrap ($SD_b$) with ten runs of the neutral model and also the SD from the regressions used to estimate $D_q$ ($SD_r$). Contrary to the results observed with p-model, the $SD_b$ and consequently the CI for the neutral model were greater than the $SD_r$ (Table S3). This means that we cannot always assume that $SD_r$ are smaller than $SD_b$. 

I analyze two stages of the model: an early stage where species are filling the available space thus the community is in a transient stage: the number of species (S) and the total density are changing (Figure 6 c, Table S2), the Shannon diversity index (H) does not have much variation. The convergence to a steady state can be exponentially fast in this kind of models [@Durrett1994a] and we can observe this in the model. When all sites are occupied, the zero-sum assumption of the neutral model is at work and a steady state is reached, approximately after time 20. I analyzed a second stage from time 30 to 60 (Figure 6d, Table S2). 
In the first stage there are big differences between $D_q$ from successive time intervals and these differences were diminishing (Figure 6a) but are still significant. These results are confirmed with the permutation procedure (Table S4). 
In the second stage there seems to be no differences between times (Figure 6b) and the permutations test confirm this (Table S5). Thus the neutral model reaches a steady state very fast and $D_q$ can be used to detect the differences between these stages.
The $D_q$ behavior for BCI SRS (Figure 3b), with overlapping CI for all years, is more similar to $D_q$ of the neutral model at a steady state (Figure 6d). This reinforces the previous results about transient states in BCI.

I also analyzed the effect of two dispersal kernels in the model's spatial pattern. I use the same mean dispersal distance for an exponential and an inverse power dispersal distributions, the other parameters were the same as in the previous simulations. In Figure 7 we see that the spatial patterns are very different and the $D_q$ have significant differences taking into account the CI and the permutation procedure. The exponential kernel produces a spatial pattern with well defined boundaries between patches (Figure 7a). Instead the power law kernel produces intermingled species patterns with no easily defined patches. The power kernel produces a spatial pattern more similar to BCI thus the $D_q$'s are closer, a better fit could be achieved adjusting the dispersal distance and increasing the grid resolution to match the one of the BCI data. 


# Discussion

Multifractal spectrum $D_q$ can be used to describe spatial patterns in all the cases I analyzed. Multifractals patterns could be produced by the existence of multiplicative interaction between species but this is not the only possibility, spatially correlated random processes like dispersal and growth to adjacent areas would also produce multifractals [@Stanley1988]. Because plant and animal species are generally aggregated in space is very likely that multifractal analysis can be used in a wide range of cases.
The occurrence of multifractal patterns with the same $D_q$ in several places does not prove that the same mechanism is acting, but may provide stronger evidence for similar mechanisms than if only one fractal dimension is estimated. Usually this dimension corresponds to $D_0$, that characterizes the geometric complexity of the spatial distribution, the other dimensions in the multifractal spectrum characterize the non-uniformity of the distribution [@Harte2001]. This adds much information to the characterization of spatial patterns that is lost if only one dimension is calculated.

I calculated $D_q$ for BCI using the spatial rank surface (SRS) and biomass, besides there seems to be no differences between censuses. Then, what is the interpretation of the spectra? $D_q$ from SRS is more asymmetric, the negative part ($q<0$) is more important, this means that the fluctuations of low values are more important, and low values represent more abundant species. This is directly related to the shape of the non spatial rank-abundance curve, if the community were more even $D_q$ will be more symmetrical. But the interrelation of the species' spatial distributions also influences $D_q$, the patchiness of the less abundant species with respect to the most abundant will change the curvature and the range of $D_q$ in the negative part. The positive part ($q>0$) reflects the rare species and also is influenced by the patchiness. The relationship of $D_q$ with patchiness is revealed with the different dispersal kernels used with the neutral model. The exponential kernel produced well defined patches and $D_q$ is more flat, the inverse power kernel produced more intermingled distributions and $D_q$ have a greater range. Multifractal spectra calculated from biomass is more symmetrical and that is related to the frecuency distribution of high and low biomass values and its spatial distribution.

The analysis of SRS using $D_q$ adds a new dimension to the comparison of species composition, because at the same time is taking into account the spatial distribution and species-abundances. Comparisons could be done for different sample sizes, taking into account that they have a relatively good fit and that the process that act are comparable. For time series of biomass or densities $D_q$ could be also applied. For multispecies time series the Borda-de-Agua [-@Borda-de-Agua2002] should be used.

I relate the variability of the spatial distribution of biomass with richness (S) and Shannon diversity (H) and it is evident that this is a factor that sets an upper limit but there are other important factors involved. Previous studies found a negative relationship between $D_1$ and H [@Alados2007, @Saravia2012], but here I found evidence of a positive relationship. That means that H can increase when the spatial distribution of species is more homogeneous. In contrast S can increase with spatial variability. This can be explained as follows: as H responds more to the variation of abundant species than S, the dominance of a few species produces spatial homogeneity and spatial variability permits the appearance of more species. But also there are other possible explanations, an analysis of $D_q$ in relation with the full species abundance distribution (not only H and S) should be worthwhile and functional diversity should also be considered [@Weigelt2008].


For the two models I calculate SD with two methods: a) from the regressions used to estimate $D_q$, this means that I used the spatial distributions to calculate it ($SD_r$). b) Using bootstrap with a set of $D_q$ obtained with different realization or simulations of the model ($SD_b$), these are called ensemble statistics. Studying natural processes where ensemble statistics can not be obtained, is usually assumed than spatially calculated statistics approximate ensemble statistics, but this assumption is not generally valid for multifractals [@Marshak1997]. In an ecological context this is called a space-for-time substitution [@Blois2013], I observed differences in SD not in averages but the issue of multifractals spatial distributions and the space for time substitutions deserves more investigation.
When the model is an exact multifractal like p-model $SD_r$ is three times smaller than $SD_b$, this is the expected result because the estimation method for $D_q$ uses autocorrelated data. When I do the same comparison for the neutral model $SD_b$ is ten times smaller than $SD_r$. This implied that the variability between repeated simulations of the model is much lower than the estimated using spatial information. Then as we usually have only the spatial information this can lead us to the erroneous conclusion that we are observing equal $D_q$ while in reality they were distinct. 


# Acknowledgments

I am grateful to the National University of General Sarmiento for financial support and to Fernando R. Momo for our great conversations about ecological theory. I also wish to thank to Diana Marco and Luis Borda-de-Água for useful suggestions to the manuscript, and Graeme Ruxton for his english revision. 

# References  