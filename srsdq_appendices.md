# Appendices

## Online Appendix A: Model Description 

This model represent a continuum between hierarchical and neutral model in the same spirit as [@a2]. The model is a stochastic cellular automata (CA) or also called interactive particle system [@a1]. In these kind of models space is discretized into a grid and only one individual can occupy a particular position. Each position represents an area fixed by the investigator to represent the real system. 

In this model all individuals have almost the same parameters, besides they should belong to different species [@a3], and each species is assigned with a number. There are only two possible differences between species: 

* They have a different frequency in the metacommunity and also can have different abundances in the local community.

* Hierarchical competition:  species with lower numbers have a probability to replace species with higher numbers as in [@a5]. 

The size of the community is given by *J = dimX* x *dimY*, where *dimX* and *dimY* are the dimension of the grid. Thus *J* is the maximum number of individuals in a given area. 

There are four processes included in the model: death, local dispersal, and migration, starting with an empty site the following events can happen:

(1) With probability *m* an individual of a species *i* can migrate from the metacommunity at a rate proportional to its frequency $X_i$ in the metacommunity.

(2) When the grid is not full, individuals give birth with rate 1 to a new individual that disperse to the neighborhood with a dispersal kernel. Here I use an exponential kernel but an inverse power kernel [@a4] is also implemented in the software.

(2) Individuals die a rate $\mu$

(5) When an individual dies it is replaced by a migrant from metacommunity with probability $m$ and with probability $1-m$ by an individual from the neighborhood. The neighborhood is established using the dispersal kernel with average distance $d$. Once the grid is full it keeps full, because when an individual dies is immediately replaced by another. This is called the zero-sum assumption in neutral models. 

(6) If the individual does not die it can be replaced by an individual from the metacommunity or neighborhood as in (5), but an individual of species with number $k$ can replace and individual of a species $k+1$ with probability $\rho$. Thus a hierarchical ordering of species is established. When this probability is zero the model behavior is neutral.

The possible dispersal kernels are:

1. Exponential $p(x) = \lambda e^{-\lambda x}$ with $mean=1/\lambda$ where $x\ge 0$.

2. Power law as defined in [@a4], $p(x) =  \frac{\alpha -1}{x_{min}} \left(\frac{x}{x_{min}} \right)^{-\alpha}$ with $mean =\frac{\alpha-1}{\alpha-2}x_{min}$ where $\alpha > 1$ and $x \ge x_{min}$. In all cases I used $x_{min} = 1$.
   
The C++ source code of the model is available at <https://github.com/lsaravia/neutral> including some additional features like non-saturated and biomass simulation of individuals.




## Appendix A References


@a1. Durrett R, Levin SA (1994) Stochastic spatial models: a user’s guide to ecological aplications. Philosophical transactions of the Royal Society of London Series B 343: 329–350.

@a2. Gravel D, Canham CD, Beaudet M, Messier C (2006) Reconciling niche and neutrality: the continuum hypothesis. Ecol Lett 9: 399–409. doi:10.1111/j.1461-0248.2006.00884.x.

@a3. Hubbell SP (2001) The unified neutral theory of biodiversity and biogeography. Princeton University Press. 375p.

@a4. Marco DE, Montemurro MA, Cannas SA (2011) Comparing short and long-distance dispersal: modelling and field case studies. Ecography 34: 671–682. doi:10.1111/j.1600-0587.2010.06477.x.

@a5. Tilman D (1994) Competition and biodiversity in spatially structured habitats. Ecology 75: 2–16.

\newpage

## Appendix B Tables

+------------+-----------+------------+-----------+------+-----------+-----------+
| $D_q$ Type |    SAD    |  Spatial   | Number of | Side | Freq > 60 | Freq > 90 |
|            |           |  Pattern   |  Species  |      |           |           |
+============+===========+============+===========+======+===========+===========+
| DqSAD      | Logseries | Randomized |         8 |  256 |      0.87 |      0.45 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |         8 |  512 |      0.75 |      0.46 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  256 |      0.99 |      0.51 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  512 |      0.91 |      0.50 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  256 |      1.00 |      0.51 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  512 |      1.00 |      0.50 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           | Regular    |         8 |  256 |      0.90 |      0.63 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |         8 |  512 |      0.82 |      0.54 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  256 |      0.98 |      0.64 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  512 |      0.97 |      0.57 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  256 |      1.00 |      0.85 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  512 |      1.00 |      0.72 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            | Uniform   | Randomized |         8 |  256 |      0.43 |      0.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |         8 |  512 |      0.34 |      0.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  256 |      0.69 |      0.18 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  512 |      0.57 |      0.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  256 |      1.00 |      0.33 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  512 |      1.00 |      0.23 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           | Regular    |         8 |  256 |      0.99 |      0.30 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |         8 |  512 |      0.53 |      0.38 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  256 |      1.00 |      0.79 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  512 |      1.00 |      0.76 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  256 |      1.00 |      0.82 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  512 |      1.00 |      0.77 |
+------------+-----------+------------+-----------+------+-----------+-----------+
| DqSRS      | Logseries | Randomized |         8 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |         8 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           | Regular    |         8 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |         8 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            | Uniform   | Randomized |         8 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |         8 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           | Regular    |         8 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |         8 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |        64 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  256 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+
|            |           |            |       256 |  512 |      1.00 |      1.00 |
+------------+-----------+------------+-----------+------+-----------+-----------+

Table: Proportion of $R^2$ that are greater than 0.6 and 0.9. The $R^2$ are from linear regressions used to estimate generalized dimensions $D_q$ in simulated spatial patterns. $D_q$ was estimated for different spatial patterns: *Regular* & *Randomized*, with two species abundance distributions: *Logseries* & *Uniform*, and using two approaches: *DqSRS* and *DqSAD* (See main text). *Side* is the side of grid used for simulate the spatial pattern.

\newpage


+------------+------+---------------+---------+--------+---------+---------+
| $D_q$ type | Side | Metacommunity |   Mean  | $\rho$ | Freq>60 | Freq>90 |
|            |      |    Species    | Species |        |         |         |
+============+======+===============+=========+========+=========+=========+
| DqSAD      |  256 |            11 |     6.1 |  0.000 |    0.98 |    0.81 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     6.1 |  0.001 |    1.00 |    0.81 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     6.3 |  0.010 |    0.97 |    0.72 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     6.8 |  0.100 |    0.97 |    0.64 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     5.0 |  1.000 |    0.90 |    0.49 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |            86 |    41.9 |  0.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    42.8 |  0.001 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    45.7 |  0.010 |    1.00 |    0.72 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    27.9 |  0.100 |    0.96 |    0.60 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    10.7 |  1.000 |    0.89 |    0.49 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |           341 |   109.9 |  0.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |   112.2 |  0.001 |    1.00 |    0.99 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |   121.1 |  0.010 |    1.00 |    0.70 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    51.0 |  0.100 |    0.98 |    0.65 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    13.9 |  1.000 |    0.91 |    0.51 |
+------------+------+---------------+---------+--------+---------+---------+
| DqSAD      |  512 |            11 |     6.8 |  0.000 |    1.00 |    0.81 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     7.0 |  0.001 |    0.99 |    0.81 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     7.3 |  0.010 |    0.95 |    0.80 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     7.9 |  0.100 |    0.87 |    0.53 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     3.7 |  1.000 |    0.83 |    0.54 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |            86 |    44.5 |  0.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    45.1 |  0.001 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    46.6 |  0.010 |    1.00 |    0.64 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    42.8 |  0.100 |    0.85 |    0.52 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    17.8 |  1.000 |    0.83 |    0.49 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |           341 |   146.5 |  0.000 |    1.00 |    0.99 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |   147.8 |  0.001 |    1.00 |    0.99 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |   159.9 |  0.010 |    1.00 |    0.60 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    99.5 |  0.100 |    0.87 |    0.52 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    36.7 |  1.000 |    0.84 |    0.50 |
+------------+------+---------------+---------+--------+---------+---------+
| DqSRS      |  256 |            11 |     6.1 |  0.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     6.1 |  0.001 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     6.3 |  0.010 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     6.8 |  0.100 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     5.0 |  1.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |            86 |    41.9 |  0.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    42.8 |  0.001 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    45.7 |  0.010 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    27.9 |  0.100 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    10.7 |  1.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |           341 |   109.9 |  0.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |   112.2 |  0.001 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |   121.1 |  0.010 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    51.0 |  0.100 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    13.9 |  1.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
| DqSRS      |  512 |            11 |     6.8 |  0.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     7.0 |  0.001 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     7.3 |  0.010 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     7.9 |  0.100 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |     3.7 |  1.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |            86 |    44.5 |  0.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    45.1 |  0.001 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    46.6 |  0.010 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    42.8 |  0.100 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    17.8 |  1.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |           341 |   146.5 |  0.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |   147.8 |  0.001 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |   159.9 |  0.010 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    99.5 |  0.100 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+
|            |      |               |    36.7 |  1.000 |    1.00 |    1.00 |
+------------+------+---------------+---------+--------+---------+---------+

Table: Proportion of $R^2$ that are greater than 0.6 and 0.9. The $R^2$ are from linear regressions used to estimate generalized dimensions $D_q$ in simulated neutral spatial patterns. $D_q$ was estimated using two approaches: *DqSRS* and *DqSAD* (See main text) for communities with different degree of neutrality given by the parameter *Replacement*.  I used a Logseries metacommunity with 11, 86 and 341 species and *Side* is the side of the simulation grid.

\newpage
   
+------+---------------+-------------+----------+
| Side | Metacommunity | Replacement | Mean No. |
|      |  No. Species  |    $\rho$   | Species  |
+======+===============+=============+==========+
|  256 |            11 |       0.000 |      6.1 |
+------+---------------+-------------+----------+
|  256 |            11 |       0.001 |      6.1 |
+------+---------------+-------------+----------+
|  256 |            11 |       0.010 |      6.3 |
+------+---------------+-------------+----------+
|  256 |            11 |       0.100 |      6.8 |
+------+---------------+-------------+----------+
|  256 |            11 |       1.000 |      5.0 |
+------+---------------+-------------+----------+
|  256 |            86 |       0.000 |     41.9 |
+------+---------------+-------------+----------+
|  256 |            86 |       0.001 |     42.8 |
+------+---------------+-------------+----------+
|  256 |            86 |       0.010 |     45.7 |
+------+---------------+-------------+----------+
|  256 |            86 |       0.100 |     27.9 |
+------+---------------+-------------+----------+
|  256 |            86 |       1.000 |     10.7 |
+------+---------------+-------------+----------+
|  256 |           341 |       0.000 |    109.9 |
+------+---------------+-------------+----------+
|  256 |           341 |       0.001 |    112.2 |
+------+---------------+-------------+----------+
|  256 |           341 |       0.010 |    121.1 |
+------+---------------+-------------+----------+
|  256 |           341 |       0.100 |     51.0 |
+------+---------------+-------------+----------+
|  256 |           341 |       1.000 |     13.9 |
+------+---------------+-------------+----------+
|  512 |            11 |       0.000 |      6.8 |
+------+---------------+-------------+----------+
|  512 |            11 |       0.001 |      7.0 |
+------+---------------+-------------+----------+
|  512 |            11 |       0.010 |      7.3 |
+------+---------------+-------------+----------+
|  512 |            11 |       0.100 |      7.9 |
+------+---------------+-------------+----------+
|  512 |            11 |       1.000 |      3.7 |
+------+---------------+-------------+----------+
|  512 |            86 |       0.000 |     44.5 |
+------+---------------+-------------+----------+
|  512 |            86 |       0.001 |     45.1 |
+------+---------------+-------------+----------+
|  512 |            86 |       0.010 |     46.6 |
+------+---------------+-------------+----------+
|  512 |            86 |       0.100 |     42.8 |
+------+---------------+-------------+----------+
|  512 |            86 |       1.000 |     17.8 |
+------+---------------+-------------+----------+
|  512 |           341 |       0.000 |    146.5 |
+------+---------------+-------------+----------+
|  512 |           341 |       0.001 |    147.8 |
+------+---------------+-------------+----------+
|  512 |           341 |       0.010 |    159.9 |
+------+---------------+-------------+----------+
|  512 |           341 |       0.100 |     99.5 |
+------+---------------+-------------+----------+
|  512 |           341 |       1.000 |     36.7 |
+------+---------------+-------------+----------+

Table: Mean number of species for simulated communities with different degree of neutrality given by the parameter *Replacement*.  I used a Logseries metacommunity with 11, 86 and 341 species and *Side* is the side of the simulation grid.

\newpage         


+------+---------------+----------+-------+-------+----------+--------+-----------+
| Side | Metacommunity | Mean No. |  Type | Power | n(Power) | Type I | n(Type I) |
|      |  No. Species  | Species  |       |       |          | Error  |           |
+======+===============+==========+=======+=======+==========+========+===========+
|  512 |            11 |     5.54 | SAD   | 0.239 |     1000 |      0 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  512 |            11 |     5.54 | DqSRS | 0.676 |     1000 |  0.058 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  512 |            11 |     5.54 | DqSAD | 0.964 |     1000 |   0.44 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  512 |            86 |       34 | SAD   | 0.687 |     1000 |  0.013 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  512 |            86 |       34 | DqSRS |   0.7 |     1000 |  0.089 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  512 |            86 |       34 | DqSAD | 0.761 |     1000 |  0.218 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  512 |           341 |      108 | SAD   | 0.756 |     1000 |  0.009 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  512 |           341 |      108 | DqSRS |   0.7 |     1000 |  0.089 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  512 |           341 |      108 | DqSAD | 0.885 |     1000 |  0.191 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  256 |            11 |     5.32 | SAD   | 0.012 |     1000 |      0 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  256 |            11 |     5.32 | DqSRS | 0.638 |     1000 |  0.116 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  256 |            11 |     5.32 | DqSAD | 0.865 |     1000 |  0.431 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  256 |            86 |    31.84 | SAD   | 0.662 |     1000 |  0.018 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  256 |            86 |    31.84 | DqSRS | 0.702 |     1000 |  0.116 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  256 |            86 |    31.84 | DqSAD | 0.914 |     1000 |  0.360 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  256 |           341 |     78.3 | SAD   | 0.793 |     1000 |  0.044 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  256 |           341 |     78.3 | DqSRS | 0.703 |     1000 |  0.191 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+
|  256 |           341 |     78.3 | DqSAD | 0.909 |     1000 |  0.289 |       225 |
+------+---------------+----------+-------+-------+----------+--------+-----------+

Table: Power and Type I error rate for SAD, DqSAD and DqSRS. The number of points used for SAD is 
the number of species found in the compared communities, and the number of points for multifractal spectra are the q used for the estimation ranging -24 to 24 (n=35). The complete set used is 
q={-24,-22,-20,-18,-16,-14,-12,-10,-8,-6,-4,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4,6,8,10,12,14,16,18,20,22,24
}.

\newpage         

## Appendix figures


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{Dq_UniLog_NumSp_512.png}
\caption{Generalized dimension spectra $D_q$ of simulated species spatial patterns. The points are means of 10 repetitions of simulated patterns using a spatial grid of side=512. A Logseries or uniform species abundance distribution were used, with 8,64 and 256 species. Two forms of generalized dimensions were estimated: DqSRS, from species rank surface $D_q^{SRS}$ and DqSAD, estimated from species abundance distribution $D_q^{SAD}$. Two different spatial patterns were used: Regular, the species are distributed in vertical bands, Randomized the spatial distribution of species was randomized.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{DqFit_Unif_8_256.png}
\caption{Linear fit from generalized dimension ($D_q$) estimation showing a range of $q$ from -4 to 4. The spatial grid has a side=256 occupied with 8 species with a uniform abundance distribution. Two different spatial patterns were used: a) Regular, a regular spatial pattern with species distributed in vertical bands of equal width. b) Randomized, the positions of species in the grid are randomized. Two kinds of generalized dimension were estimated: DqSRS corresponds the fit of $D_q^{SRS}$ (see text) and DqSAD is the fit from the estimation of $D_q^{SAD}$ (see text). $Z_q(\epsilon)$ corresponds to the partition function calculated for a box with side $\epsilon$ in the DqSRS case. In DqSAD case, $\epsilon$ represent the area of the box used.}
\end{figure}. 

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{DqFit_Unif_256_256.png}
\caption{Linear fit from generalized dimension ($D_q$) estimation showing a range of $q$ from -4 to 4. The spatial grid has a side=256 occupied with 256 species with a uniform abundance distribution. Two different spatial patterns were used: a) Regular, a regular spatial pattern with species distributed in vertical bands of equal width. b) Randomized, the positions of species in the grid are randomized. Two kinds of generalized dimension were estimated: DqSRS corresponds the fit of $D_q^{SRS}$ (see text) and DqSAD is the fit from the estimation of $D_q^{SAD}$ (see text). $Z_q(\epsilon)$ corresponds to the partition function calculated for a box with side $\epsilon$ in the DqSRS case. In DqSAD case, $\epsilon$ represent the area of the box used.}
\end{figure}. 

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{DqFit_Logser_8_256.png}
\caption{Linear fit from generalized dimension ($D_q$) estimation showing a range of $q$ from -4 to 4. The spatial grid has a side=256 occupied with 8 species with a uniform logseries distribution. Two different spatial patterns were used: a) Regular, a regular spatial pattern with species distributed in vertical bands of equal width. b) Randomized, the positions of species in the grid are randomized. Two kinds of generalized dimension were estimated: DqSRS corresponds the fit of $D_q^{SRS}$ (see text) and DqSAD is the fit from the estimation of $D_q^{SAD}$ (see text). $Z_q(\epsilon)$ corresponds to the partition function calculated for a box with side $\epsilon$ in the DqSRS case. In DqSAD case, $\epsilon$ represent the area of the box used.}
\end{figure}. 


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{DqFit_Logser_64_256.png}
\caption{Linear fit from generalized dimension ($D_q$) estimation showing a range of $q$ from -4 to 4. The spatial grid has a side=256 occupied with 64 species with a uniform logseries distribution. Two different spatial patterns were used: a) Regular, a regular spatial pattern with species distributed in vertical bands of equal width. b) Randomized, the positions of species in the grid are randomized. Two kinds of generalized dimension were estimated: DqSRS corresponds the fit of $D_q^{SRS}$ (see text) and DqSAD is the fit from the estimation of $D_q^{SAD}$ (see text). $Z_q(\epsilon)$ corresponds to the partition function calculated for a box with side $\epsilon$ in the DqSRS case. In DqSAD case, $\epsilon$ represent the area of the box used.}
\end{figure}. 


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{DqFit_Logser_256_256.png}
\caption{Linear fit from generalized dimension ($D_q$) estimation showing a range of $q$ from -4 to 4. The spatial grid has a side=256 occupied with 256 species with a uniform logseries distribution. Two different spatial patterns were used: a) Regular, a regular spatial pattern with species distributed in vertical bands of equal width. b) Randomized, the positions of species in the grid are randomized. Two kinds of generalized dimension were estimated: DqSRS corresponds the fit of $D_q^{SRS}$ (see text) and DqSAD is the fit from the estimation of $D_q^{SAD}$ (see text). $Z_q(\epsilon)$ corresponds to the partition function calculated for a box with side $\epsilon$ in the DqSRS case. In DqSAD case, $\epsilon$ represent the area of the box used.}
\end{figure}. 


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{Dq_NeuLog_NumSp_512.png}
\caption{Generalized dimension spectra $D_q$ of spatial patterns generated with a spatial neutral/hierarchical model. \textit{Replacement} is the parameter that determines the degree of neutrality. When this parameter is 0 the model is completely neutral and there is no competitive replacement of species. When \textit{Replacement} is 1 competitive superior species always replaces inferior ones and the model is completely hierarchical. The simulations use a metacommunity with a logseries abundance distribution with 11,86 and 341 species. The simulation grid side is 512, and the other parameters are: MortalityRate=0.2, DispersalDistance=0.4 (2.5 grid units), ColonizationRate=0.001.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{powADq10_DifR_512.png}
\caption{Power of the Anderson-Darling to test the hipothesis that $D_q$ curves and SAD of simulated neutral/hierarchical model are different. The parameter $Replacement$ determines the degree of neutrality. When this parameter is 0 the model is completely neutral and there is no competitive replacement of species. When \textit{Replacement} is 1 competitive superior species always replaces inferior ones and the model is completely hierarchical. $\Delta Replacement$ is the difference in the parameter, when is lower the communities are more similar. Simulations use a metacommunity with a logseries abundance distribution with 11,86 and 341 species. The simulation grid side is 512, and the other parameters are given in the main text.}
\end{figure}
