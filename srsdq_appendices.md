# Appendixes

## Online Appendix A: Model Description 

This model is a stochastic cellular automata (CA) or also called interactive particle system [1]. In these kind of models space is discretized into a grid and only one individual can occupy a particular position. Each position represents an area fixed by the investigator to represent the real system. 

As this is a neutral model all individuals have the same parameters, besides they should belong to different species [2]. The only difference between species is that they have a different frequency in the metacommunity and also can have different abundances in the local community. The size of the community is given by *J = dimX* x *dimY*, where *dimX* and *dimY* are the dimension of the grid. Thus we have a maximum number of individuals in a given area. 

There are four processes included in the model: death, local dispersal, and migration, starting with an empty site the following events can happen:

(1) With probability *m* an individual of a species *i* can migrate from the metacommunity at a rate proportional to its frequency $X_i$ in the metacommunity.

(2) When the grid is not full individuals give birth to a new individual that disperse to the neighboorhood with a dispersal kernel. Here I use an exponential kernel but an inverse power kernel [3] is also implemented in the software.

(2) Individuals die a rate $\mu$

(5) When an individual dies it is replaced by a migrant from metacommunity with probability *m* and with probability *1-m* or by an individual of the neighborhood. The neighborhood is established using the dispersal kernel with average distance *d*. Once the grid is full it keeps full, because when an individual dies is immediately replace by another. This is called the zero-sum assumption. 

The dispersal kernels used are:

1. Exponential $p(x) = \lambda e^{-\lambda x}$ with $mean=1/\lambda$ where $x\ge 0$.

2. Power law $p(x) =  \frac{\alpha -1}{x_{min}} \left(\frac{x}{x_{min}} \right)^{-\alpha}$ with $mean =\frac{\alpha-1}{\alpha-2}x_{min}$ where $\alpha > 1$ and $x \ge x_{min}$. In all cases I used $x_{min} = 1$.
   

The C++ source code of the model is available at <http:github/lasaravia/neutral> including some additional features like non-saturated model and a hierarchical competition model.

Neutral: rule (5) determines that no species has an advantage 
over another and as all species have the same rates the model is 
termed neutral. 

Hierarchical: rule (5) is modified, an individual of a species 
k always replaces and individual of a species k+1 with probability 
ρ. Thus a hierarchical ordering of species is established. When 
this probability is zero the neutral model is recovered. 


## Appendix A References

1. Durrett R, Levin SA (1994) Stochastic spatial models: a user’s guide to ecological aplications. Philosophical transactions of the Royal Society of London Series B 343: 329–350.

2. Hubbell SP (2001) The unified neutral theory of biodiversity and biogeography. Princeton University Press. 375p.

3. Marco DE, Montemurro MA, Cannas SA (2011) Comparing short and long-distance dispersal: modelling and field case studies. Ecography 34: 671–682. doi:10.1111/j.1600-0587.2010.06477.x.
\newpage
