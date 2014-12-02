# Species rank distributions in space: a multifractal analysis

# Combining species abundances and spatial patterns to compare communities

**Leonardo A. Saravia**, Ph.D.

Instituto de Ciencias Básicas

Universidad Nacional de General Sarmiento

J.M. Gutierrez 1159 (1613), Los Polvorines

Buenos Aires, Argentina.

<lsaravia@ungs.edu.ar>

\newpage

## Abstract

Species-area relationships (SAR) and species abundance distributions (SAD) are between the most studied patterns in ecology, due to its application in both theoretical and conservation issues. One problem with these general patterns is that different theories can generate the same predictions, for this reason these can not be used to detect different mechanisms. To circumvent this, other more sensitive patterns  should be developed. One candidate for such pattern is the generalization of SAR. Instead of calculating the relationship of number of species with area, the whole species abundance distribution is used. This leads to the application of generalized dimensions ($D_q$) to study species diversity-area relationships. The calculation of $D_q$ is then based on the species-abundance distribution (SAD). Here I introduce a new way to study species diversity-area relationships: a robust way to express SAD are rank-abundance distribution (RAD), I use spatial version of RAD: the species-rank surface (SRS) that can be analyzed using $D_q$. To demonstrate the relationship of the different $D_q$ versions with species abundances and spatial patterns I perform simulations of spatial patterns with different SAD. Finally I compare the power of the $D_q$'s, SAD and SAR to detect different community patterns using a continuum of hierarchical and neutral spatially explicit models.

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

I call $D_q$ calculated from species abundances $D_q^{SAD}$ and $D_q$ calculated from species rank surface $D_q^{SRS}$, when I mention $D_q$ without superscript I refer to both. $D_q^{SAD}$ represent the scaling of the Hill's generalized diversity index [@Hill1973], when the moment order is $q=0$,then $D_q^{SAD}$ becomes the exponent of the SAR power-law scaling, when $q=1$, $D_q^{SAD}$ represent the scaling of Shannon diversity index and when q=2, $D_q^{SAD}$ becomes the scaling of Simpson's index. This is why $D_q$ can characterize diversity-area relationships.

Theoretically $D_q$ must be a non-increasing function of $q$ [@Hentschel1983], which means that if $q_1 \geq q_2$ then $D_{q1} \leq D_{q2}$. Some studies showed small violations to this property for $D_q^{SAD}$ [@Borda-de-Agua2002; @Zhang2006]. These violations are be related to the way that $D_q^{SAD}$ is defined: the summation of equation @e4 is over species and the summation of equation @e1 is over boxes, this changes the way in which the mathematical limits are taken and also the way of computation of $D_q^{SAD}$. A partial solution has been proposed [@Yakimov2014], but the anomalies observed may be related to the mathematical assumptions needed for $D_q$ to be non-increasing, a new mathematical proof should be developed for $D_q^{SAD}$. Thus as long as the linear relationship is reasonable I take $D_q^{SAD}$ as a useful technique of analysis. 

I proposed a new way to analyze species-abundance-area using multifractals that fits more closely to the original definitions of equations @e1 - @e3: the species-rank surface [@Saravia2014]. To construct the species-rank surface (SRS) the spatial distribution of species have to be transformed assigning to each species position its rank. First I use the species abundances, at the whole plot level, to calculate the rank ordering the species from biggest to lowest and assigning a number starting with one. Then the rank is assigned to the spatial position of the individuals of each species forming a surface. This landscape have valleys formed by the most abundant species and peaks determined by the most rare species. Finally the standard multifractal analysis is applied. If sampling was made in quadrats without taking the spatial position of individuals, the sum of the ranks of the species in the smallest quadrats could be used to form the SRS. 

I use the coefficient of determination ($R^2$) as a descriptive measure of goodness of fit [@Borda-de-Agua2002]. The source code in C++ to perform multifractal analysis is available at <https://github.com/lsaravia/mfsba>. 

### Generalized dimension relationship with spatial patterns and SADs

I simulated species spatial patterns with different SAD's to demonstrate how $D_q$ is related to them. First I used a uniform SAD, in this case all species have the same densities and the same number of individuals. To add a degree of stochasticity I take the number of individuals of each species from a Poisson distribution with the same mean. I distributed them in bands over a spatial grid so they form a regular spatial pattern. Each position of the grid is occupied by one individual and I choose the number of species to exactly divide the side of the grid so all species are strips with approximately the same width (Figure 1).
I used square grids with sides of 256 and 512 sites which contain 65536 and 262144 individuals respectively, and 8, 64 and 256 species. I calculated $D_q$ for the regular pattern, then I randomized the positions of species to compare $D_q$ obtained with these two extreme cases.
The second SAD I used is a Logseries [@Fisher1943] with the same number of species and the same sides as previously. I used the R package untb [@Hankin2007] to calculate the density for each species, basically they use Poisson distribution with the expected Logseries abundances as means. I then build the regular pattern with strips of species, but as species have different abundances the widths for each species are different (Figure 1). Then I estimated $D_q$ for the regular and randomized patterns. I simulated 10 spatial patterns for each case and calculated the mean and standard deviation of $D_q$. 

### Spatially explicit model

To simulate more realistic patterns of species-abundance-area relationships I used a stochastic spatially explicit model. I developed a stochastic cellular automata [@Molofsky2004] that can switch between a neutral or hierarchical competition representing a continuum between niche and neutral communities [@Gravel2006]. In a neutral model individuals do not interact, and have all the same mortality, colonization rates and dispersal distances. In spite of these gross simplifications neutral models are capable of predict several real community patterns [@Rosindell2011]. At the other end of the continuum are niche communities represented by a hierarchical competition model [@Tilman1994]. In this case species have differences that imply a competitive hierarchy in which some species are always better than others, then they produce competitive exclusion [@Chave2002]. Basically I added to the neutral model a probability of replacement $\rho$. When $\rho=1$ more competitive species always replace less competitive and the model behaves as pure hierarchical. When $\rho=0$, there is no replacement of species and the model is completely neutral. A more thorough description of the model is given in appendix A and the C++ source code of the model is available at <https://github/lasaravia/neutral> and figshare <http://dx.doi.org/10.6084/m9.figshare.969692>. 

Following a classical neutral scheme the model has a metacommunity: a regional collection of communities, from which migration occurs at a rate $m$. Species can also disperse locally and I assume an exponential dispersal kernel with average dispersal distance $d$. Other model parameters are the mortality rate $\mu$, the number of species in the metacommunity and also the size of the community, this last is represented as the *side* of the grid used in the simulations. I use a logseries SAD for the metacommunity, defined by the maximum number of individuals (*side x side*) and the number of species [@Fisher1943].  

The values of the parameters were in the range of estimated for BCI from the existing literature [@Anand2010; @Condit2002; @Etienne2007]. I performed 50 simulations for each combination of parameters given in Table 1. To compute the power I made comparisons of communities with different levels of $\rho$, representing more neutral or hierarchical communities, the other parameters were kept constant. I also made comparisons between repetitions with the same $\rho$ to calculate the type I error. 


### Statistical comparison of methods

I analyzed the performance of two kinds of methods to differentiate communities. The first consist of a set of points or curves: species abundance distributions (SAD), generalized dimensions $D_q^{SAD}$ and $D_q^{SRS}$. For these I used a non-parametric test related to the Kolmogorov-Smirnov test: the Anderson-Darling (AD) test [@Feigelson2012]. This test measure the differences between the empirical distribution functions (EDF) of two datasets as a weighted sum of square deviations between the EDFs. In extensive simulations the AD test has proven more sensitive than the Kolmogorov-Smirnov test [@Stephens1974]. I calculated p-values using randomization with 1000 repetitions. I use the package kSamples [@Scholz2012] in the R statistical statistical language [@RCoreTeam2014], the scripts for all the analysis are available at fighsare (DOI:XXX)	. 

The second kind are based on a single dimension or power exponent: the SAR exponent and the information dimension.  The SAR exponent is part of the $D_q^{SAD}$ spectra when $q=0$ [@Borda-de-Agua2002]. An equivalent single number measure from $D_q^{SRS}$ is the information dimension [@Ricotta2000;@Chappell2001], that is the $D_q^{SRS}$ when $q=1$. I calculate the power of these with a T-test using the standard deviation (SD) obtained from the box-counting regressions. These SD are obtained with autocorrelated data because small squares are nested within big squares  (See Multifractal Analysis). The consequence is that the SD may be underestimated, but the slopes estimates are still unbiased [@Kutner2005]. This should result in an increased type I error rate and also in a spurious increase in power. 

### Calculation of power and type I error

I simulated communities with different degree of neutral/hierarchical structure given by the parameter $\rho$ of the model. The power of a test is the ability to reject the null hypothesis ($H_0$) when it is false. The significance level to reject $H_0$ was set a priori at $\alpha$ = 0.05 in all cases, and the rejection rate of each test was calculated as the proportion of P values that were less than or equal to $\alpha$. To estimate power I used independent simulations of communities (50 repetitions) with the same parameters except $\rho$, and compute the proportion of rejection.

The type I error is the probability to reject $H_0$ when is true (false positive). In our simulations, $H_0$ is true if two simulated communities have the same $\rho$ (and also are equal in the other parameters). To estimate type I error I compare independent simulations of communities with the same set of parameters (50 repetitions) and compute the proportion of rejection.  


# Results

Generalized dimension ($D_q$) can be interpreted like SAR power law exponent, if it is bigger the change in the number of species is greater when change the scale of observation to a larger area. $D_q$ express the change of the quantity under study when scale changes but modulated by $q$. When $q$ is positive the terms of the sums (equations 2 & 5) with more abundant species have more weight, and become even more important when $q$ is greater. When $q$ is negative we have the opposite pattern: less abundant species have more weight in the sum, so $D_q$ reflects the change of rare species. When $q$ is greater in absolute value, $D_q$ is driven by more and more extreme values thus $D_q$ will have a bigger variance.  Here I present most figures with a range of $q$ from -24 to 24 but for statistical comparisons I use a smaller range from -10 to 10 to avoid large variances.

I calculated two versions of $D_q$: a) the original definition due to Borda-de-Água [-@Borda-de-Agua2002] where $D_q$ measure the change in SAD as we change scale ($D_q^{SAD}$). b) $D_q$ based in SRS, where we measure the change in the spatial distribution of species' ranks as scale changes ($D_q^{SRS}$). $D_q$ measures the rate of change with scale from a baseline that is defined by $D_0$. When we study SAD, $D_0^{SAD}$ is the SAR exponent and its value is around 0.5. A spatial distribution of species that duplicates its number with a duplication of the side of the area studied has a value of exactly 0.5. When we study SRS the $D_0^{SRS}$ is the fractal dimension of the species rank surface, in my simulations the individuals completely fill the available space thus $D_0^{SRS}$ is equal to 2.

For the uniform SAD we expected that $D_0^{SAD}$ around 0.5 and a symmetric pattern with small deviations around this value, as all species have the same abundance and occupy the same area. The symmetric pattern is not observed in the regular cases (Figure 2) because negative part ($q<0$) analyzes numbers close to 0 and the logarithm enhances the differences between small numbers [@Laurie2011] thus the difference $\Delta D_q = \mid D_q - D_0\mid$ is greater for $q<0$. 

Theoretically $D_q$ should be a decreasing or constant and this is not observed in $D_q^{SAD}$ for the randomized spatial patterns with less species. This is because when changing scales there is a point in which no new species are found, the scaling relationship breaks. Figure 3 shows an example of $D_q$ fitted linear relationships for 64 species and a side=256 sites. The scaling for a randomized pattern $D_q^{SAD}$ breaks at 2.4, equivalent to an area of 256 units. The scaling for regular pattern $D_q^{SAD}$ shows oscillations around the fitted line but no evidence of breaks. When the number of species is higher (256) the $D_q^{SAD}$ for regular spatial pattern is closely similar to the randomized one (Figure 2), this happens because new species appear in the whole range of scales used. 

Histograms of the $R^2$ (Figure 4) indicate the presence of poor fits or a scaling break. The $D_q^{SAD}$ for randomized patterns and uniform SAD have the lowest $R^2$ of all cases. Based in all simulations a rule of thumb can be derived: 90% of $D_q$ should have an $R^2$ of 0.6 or greater and 50% should have a $R^2$ of 0.9 or greater (Appendix table 1), if not you should check the plots of the fits (Figure 3). Several patterns fail to comply this rule: all the uniform randomized patterns and the logseries randomized with 8 species metacommunity (Appendix table 1 and Appendix figures 2-4). 

The $D_q^{SAD}$ for logseries have a more symmetric pattern than for uniform SAD (Figure 2), and they have better fits with higher $R^2$  (Figure 3,4). Comparing regular and randomized spatial patterns $D_q^{SAD}$ curves are superposed or inside the SD of the other. Thus it seems that $D_q^{SAD}$ can not distinguish between them (only considering the cases where the fits are good). Moreover the range of $D_q^{SAD}$ does not change very much with the number of species, $D_q^{SAD}$ seems to depend mostly on the SAD used to generate the spatial pattern.

For $D_q^{SRS}$ the theoretical decreasing pattern is fulfilled in all cases, no anomalies were observed (Figure 2). An asymmetric pattern, like the previous case, is observed but $D_q^{SRS}$ is around 2. The asymmetry is more pronounced for patterns with uniform SAD than for logseries SAD. This is because logseries SAD have one very abundant species, several less abundant and rare species that are scattered through the pattern (Figure 1). Thus the abundant species dominates the spatial pattern and in some cases produces a greater $\Delta D_q = \mid D_q - D_0\mid$ in the positive side of the plot (Figure 2, 8 Species). 

The uniform SAD produces $D_q^{SRS}$ with higher $\Delta D_q$ for regular pattern in the $q<0$ side. This is because in the regular pattern the species are aggregated, in the randomized pattern there is no aggregation so $D_q^{SRS}$ is closer to two. Thus $D_q^{SRS}$ for regular and randomized are more different in the negative side, and more similar in the positive side. For logseries SAD the differences in $D_q^{SRS}$ are similar at negative or positive sides of $q$.  In general $D_q^{SRS}$ curves for different spatial patterns and different SAD's are distinct, except in some cases for 8 species the curves are inside the SD of a different pattern.

The $R^2$ for $D_q^{SRS}$ are all greater than 0.9 and are all greater than $D_q^{SAD}$, and all comply with the rule of thumb (Appendix table 1). The linear trends are also better (Figure 3). An example of linear trends for different number of species and different SAD is shown in the appendix (Appendix figures 2-6). The same qualitative patterns of $D_q^{SAD}$ and $D_q^{SRS}$ are observed for simulations with side=512 (Appendix figure 1).



## Simulated Neutral communities

Examples of the patterns simulated by the Neutral/hierarchical model are shown in Figure 5. More hierarchical communities with are more dark. By definition hierarchical communities have more competitive species with lower index, and neutral communities have more abundant species with higher index numbers, determined by metacommunity abundance (see appendix model description). With more degree of competitive hierarchy one or few species dominate and several rare species are scattered over the landscape (Figure 6). This produces a mostly uniform pattern of dominant species with rare species distributed at random. In neutral communities the most abundant species are not so dominant (Figure 6), they left space for species with intermediate abundances, and that produces a pattern of several aggregated species. Aggregation is produced in this model only because dispersal is mainly near the parent. 

For both estimated $D_q$ the $R^2$ are very good, $D_q^{SRS}$ shows always $R^2$>0.9 and $D_q^{SAD}$ have in almost all cases $R^2$>0.6 and a 50% or more of the cases greater than 0.9 (Appendix table 2). Thus both satisfy the rule of thumb described previously.

There are two groups of $D_q^{SAD}$: one composed of neutral like communities for $\rho < 0.1$ and another composed of more hierarchical ones for $\rho >0.1$. The curves for hierarchical communities are more separated for negative $q$ than for positive $q$. In neutral communities this pattern is inverted, for positive $q$ the curves are more different. This reflects the patterns in SAD: hierarchical communities have one or few relatively abundant species and this is why $D_q^{SAD}$ reach 0 quickly, no new abundant species are found when scale changes. Neutral communities have more species with intermediate densities and this produces $D_q^{SAD} > 0$ in the positive side. 

In theory $D_q$ have a constant value when $q$ tends infinity (negative or positive). $D_q^{SAD}$ spectra quickly reaches a constant maximum for negative $q$ and a minimum for positive $q$, this pattern is more pronounced with hierarchical communities because they tend to have two types of species: dominant ones reflected in positive side and rare species in the negative. When communities are more neutral ($\rho < 0.1$) and there are more species with intermediate densities, $D_q^{SAD}$ tends to reach the asymptotic values more slowly in the negative side.  

For $D_q^{SRS}$ a similar groups of neutral or hierarchical communities are also present. We previously saw the $D_q^{SRS}$ is more related to the spatial pattern than $D_q^{SAD}$ thus we can interpret $D_q^{SRS}$ in terms of randomness and aggregation of species.  For hierarchical communities $D_q^{SRS}$ in the negative side is very close to 2, that is the dimension of a uniform surface, rare species exert a very low influence on dominant species that have a uniform distribution. For neutral communities there are more species with low to medium densities and they have greater aggregation thus $D_q^{SRS}$ is higher. 

When $q$ is positive lower values of $D_q^{SRS}$ means more intensity of spatial pattern. Communities with $\rho=1$ are the most hierarchical, with one dominant species and and a few very rare species (Figure 6). For these communities $D_q^{SRS}$ is closer to 2, that represents the uniform distribution of dominant species. When the metacommunity have more species the local community also have more species (Appendix table 3) and  $D_q^{SRS}$ starts to differentiate from 2 at lower $q$
$D_q^{SRS}$ for the intermediate hierarchical ($\rho=0.1$) starts higher than neutral at $q$ near 0, but  crosses neutral curves and ends in the lowest place.  These communities have more species that also are more abundant but still with few individuals, this forms very sharp peaks in the SRS and produces a $D_q^{SRS}$ farther from 2. The curvature of $D_q^{SRS}$ is more pronounced when there are more species. 
For $\rho$ less than 0.1 communities are more neutral and have more species with similar densities, these form softer valleys and peaks that result in a $D_q^{SRS}$ intermediate between the two hierarchical. Simulations with side=512 have similar patterns for $D_q$ (Appendix figure 7).


## Statistical Power and type I errors 

To calculate the power of the methods I compared communities with different $\rho$, thus the alternative hypothesis was true. Instead to estimate type I error we need to compare different runs of communities simulated with identical parameters.

For $D_q^{SAD}$ and $D_q^{SRS}$ I have the possibility to use different ranges of $q$. High values of $q$ in absolute terms should produce $D_q$ with high variances resulting in higher spread of values obtained in different simulation runs. Ranges of $q$ between -10 and 10 or narrower are generally used [@Saravia2012; @Wei2013;@Laurie2011;@Yakimov2008] but some times the applied range was wider [@Saravia2012a]. I started using a $q$ range of -24 to 24, I found for this range that type I error rates are, in all cases, higher than the nominal significance level $\alpha = 0.5$ (Appendix table 4). A statistical test is valid if the type I error is lower or equal to $\alpha$ [@Edgington1995], thus to assure the validity for these methods a narrower range should be used. I then use a $q$ range between -10 and 10.

Using only one dimension of the spectra ($D_0^{SAD}$ and $D_1^{SRS}$) result in a power generally below 0.5 (Table 2) and the type I error is around 0.4, much greater than $\alpha$. These high type I error values are expected due to the presence of spatial autocorrelations in the dependent variable [@Legendre2002]. Parameter estimates can be corrected in different ways [@Legendre2012], but these procedures should not increase the power of $D_0^{SAD}$ and $D_1^{SRS}$.  

For communities with lower species numbers (11 species in the metacommunity) the comparisons made with SAD have a constant low power across $\Delta \rho$, so no matter how different are the communities as the points used in the test are the number of species the power is low (Table 3). Generalized dimensions $D_q^{SAD}$ and $D_q^{SRS}$ in contrast have a high power but Type I error is also greater than $\alpha$. One way to alleviate this problem is to check for a coincidence of the two methods SAD & $D_q$, another way would be to rise the number of points used inside the $q$ range, because $D_q$ could be calculated for any real number. I used 21 points (Table 3) but that could be increased, the only restriction is the additional computational time required. In simulated communities with more species (86 & 341 species metacommunity) the type I error falls below $\alpha$ for all the methods and overall SAD is slightly more powerful (table 3). 

Differences between communities ($\Delta \rho$) influence power. With $\Delta \rho$ less than 0.09 the power is in most cases below 0.5. Except for SAD in some communities with high species metacommunity, these are comparisons with a higher number of points (circa 100) so this results in a greater power. 
For differences greater or equal than 0.09 the power is high (over 0.75) in most cases. The exception is $\Delta \rho=0.9$, this only happens when we compare two hierarchical communities. In these cases SAD but particularly $D_q^{SAD}$ have less power, below 0.25 in some cases. Note that with $\Delta \rho<0.1$ the communities compared are more neutral with a similar number of species and SADs, comparison with $\Delta \rho>=0.1$ are between neutral and hierarchical (except for 0.9) communities different number of species and SAD. 

   
# Discussion

While $D_q^{SAD}$ measures the change in species abundance distribution with scale, $D_q^{SRS}$ represent the change in the spatial distribution of ranks of species. Thus $D_q^{SRS} is related to the spatial pattern of species and to its abundance distribution. $D_q^{SAD}$ also reflect changes in spatial pattern but my results suggest that it can not distinguish between regular and randomized spatial patterns. In contrast $D_q^{SRS}$ curves differ clearly between these patterns. 

All $D_q$ can be interpreted in terms of $q$ that modulates the weight of abundant and rare species in the distribution. $D_q$ for positive $q$ reflect more abundant species or dominance patterns in SAD, while $D_q$ for negative $q$ represent rare species patterns. An alternative way to analyze $D_q$ would be to split species in ranges of abundances and calculate $D_0^{SAD}$ or $D_1^{SRS}$, this was done for biomass and forest height spatial analysis [prefix@CitationKey] but for species distributions has drawbacks. First the species spatial distribution is analyzed as a whole and it is quite possible that the complete set of species fits very well but one or more single species do not [@Sizling2004a]. Second rare species represent a few points in space thus the estimation of $D_q$ will have a high uncertainty.  

When I compare between competitive hierarchical communities the number of species is relatively low, SAD and $D_q^{SAD}$ have a low power but $D_q^{SRS}$ maintains a high power. This highlights the ability of $D_q^{SRS}$ to detect differences in spatial patterns of rare species. 

Both $D_q$ seem to depend on the number of species and SAD, when communities are very similar in both the power to detect differences is low in most cases. 


In neutral models SAR exponent depends on speciation rate --in our case migration from metacommunity-- and local community size [@Rosindell2007; @Cencini2012]. I do not expect to find a great power using SAR exponent ($D_0^{SAD}$) because I did not vary migration parameter and did not made comparisons between different sizes. But I found high type I error rates for $D_0^{SAD}$ and the information dimension $D_1^{SRS}$. This means that the statistical methods should be improved applying a correction for autocorrelation to lower type I errors, also a greater number of boxes should be used to increase power. In most cases exist a range of different $D_q$ and that means that the distribution is a multifractal [@Stanley1988] thus species spatial distributions will not be well described by only one generalized dimension. To compare communities $D_q^{SRS}$ and $D_q^{SAD}$ represent an improvement over comparisons made with only one dimension like SAR exponent or information dimension. 

The species abundance distribution SAD is the most studied biodiversity pattern but it's generally studied at one scale, here I used the whole simulation area, and at this scale the power of SAD is at the same level that generalized dimensions. Several studies regard SAD as not very informative because many different models can produce the same patterns, but in this simulations SAD can differentiate models quite well, except for low species numbers where its power is low. 

With less species and a more equitative abundance distributions $D_q^{SAD}$ should not be applied. With less species and more dominant SAD patterns D_qSAD do not detect changes in communities.
Thus $D_q^{SAD}$ could be used to analyze diversity-area relationships but it is valid in a more limited set of cases than $D_q^{SRS}$   



 
In summary $D_q^{SRS}$ have always better fits and can be applied in all the cases simulated here. It maintains a high power and when $D_q^{SAD}$ and SAD fail. But when the number of species is around 6 some care have to be taken and the number of $q$ used should be greater than 21.   
 

# Acknowledgments

I am grateful to the National University of General Sarmiento for financial support and to Graeme Ruxton for his English revision. 

# References  