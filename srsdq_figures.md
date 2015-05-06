# Manuscript tables


---------------------------------------------------------------------------------------------
    Method      Name - Meaning   
--------------- ----------------------------------------------------------------------------- 
SAD             Species abundance distribution.

$D_q^{SAD}$     Generalized dimension spectra based on SAD: Characterize the scaling of 
                species abundances in space.
                     
$D_q^{SRS}$     Generalized dimension spectra based on SRS: Characterize the scaling of
                the spatial distribution of the ranks of species (SRS) derived from SAD.

$D_0^{SAD}$     The power exponent of the species-area relationship. This is part of 
                $D_q^{SAD}$ and characterize the scaling of richness.

$D_1^{SRS}$     Information dimension of SRS. This is part of $D_q^{SRS}$ and characterize
                the scaling of Shannon diversity index calculated on the spatial 
                distribution of species ranks.                 
---------------------------------------------------------------------------------------------

Table: Methods used to calculate the power to compare simulated communities with different degree of neutrality.

\newpage


------------------------------------------------------
  Side   No. Species   $\mu$   $d$    $m$     $\rho$  
------- ------------- ------- ----- -------- --------- 
   256            11     0.2    25    0.001        1  

   512            86                             0.1  

                 341                            0.01  

                                               0.001  

                                                   0  
------------------------------------------------------

Table: Parameters values used in the simulations of the neutral-hierarchical model 

\newpage

------------------------------------------------------------------
  Side   Metacommunity   Mean No.       Type      Power   Type I  
          No. Species    Species                          Error   
------- --------------- ---------- ------------- ------- ---------
   512              11       5.96   $D_1^{SRS}$   0.512    0.434  

                                    $D_0^{SAD}$   0.498    0.494  

                    86      36.54   $D_1^{SRS}$   0.521    0.430  

                                    $D_0^{SAD}$   0.445    0.426  

                   341     111.31   $D_1^{SRS}$   0.497    0.342  

                                    $D_0^{SAD}$   0.494    0.436  

   256              11       5.90   $D_1^{SRS}$   0.491    0.408  

                                    $D_0^{SAD}$   0.471    0.424  

                    86      32.27   $D_1^{SRS}$   0.501    0.447  

                                    $D_0^{SAD}$   0.474    0.388  

                   341      76.57   $D_1^{SRS}$   0.490    0.389  

                                    $D_0^{SAD}$   0.443    0.363  
------------------------------------------------------------------


Table: Power and Type I error rate for T-test comparison of a single dimension of the generalized spectra: the SAR exponent ($D_0^{SAD}$) and information dimension of the species rank surface ($D_1^{SRS}$). The test use the standard deviation obtained in the regressions to fit generalized dimensions. The number of comparisons to calculate the power is n=25000, and for type I error n=6125.

\newpage

------------------------------------------------------------------
  Side   Metacommunity   Mean No.       Type      Power   Type I  
          No. Species    Species                          Error   
------- --------------- ---------- ------------- ------- ---------
   512              11       5.96   SAD           0.115    0.025  

                                    $D_q^{SRS}$   0.720    0.102  

                                    $D_q^{SAD}$   0.568    0.212  

                    86      36.54   SAD           0.697    0.009  

                                    $D_q^{SRS}$   0.680    0.014  

                                    $D_q^{SAD}$   0.616    0.011  

                   341     111.31   SAD           0.830    0.039  

                                    $D_q^{SRS}$   0.688    0.000  

                                    $D_q^{SAD}$   0.609    0.017  

   256              11       5.90   SAD           0.175    0.000  

                                    $D_q^{SRS}$   0.654    0.068  

                                    $D_q^{SAD}$   0.704    0.204  

                    86      32.27   SAD           0.675    0.019  

                                    $D_q^{SRS}$   0.657    0.025  

                                    $D_q^{SAD}$   0.613    0.027  

                   341      76.57   SAD           0.799    0.035  

                                    $D_q^{SRS}$   0.670    0.030  

                                    $D_q^{SAD}$   0.610    0.048  
------------------------------------------------------------------

Table: Power and Type I error rate of Anderson-Darling statistic to test hypothesis of differences in:  species abundance distributions (SAD), generalized dimension based on SAD ($D_q^{SAD}$) and generalized dimension based on SRS ($D_q^{SRS}$). The power is calculated testing communities with different $\rho$ and type I error is calculated for communities with the same $\rho$. The p-values were estimated using 1000 randomizations. The number of points used for SAD comparisons is the number of species found in the communities. The number of points used for multifractal spectra correspond to the q in the range -10 to 10 (n=21), according to the following set q={-10,-8,-6,-4,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,4,6,8,10}. The number of comparisons for the power calculations were n=25000 except for SAD with side=500 & metacommunity species=11, where some comparison with less than 3 species were skipped (n=23800). For type I error the comparisons were n=6125, and the same exception applies (n=5846).


# Manuscript Figures

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{spat_UniLog64_256.png}
\caption{Spatial patterns generated with logseries and uniform species abundance distribution (all species have the same density) with 64 species and a grid with side=256. a) Regular: species are distributed in vertical bands. b) Randomized: the position of species are distributed at random in space.}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{Dq_UniLog_NumSp_256.png}
\caption{Generalized dimension spectra $D_q$ of simulated species spatial patterns. The points are means of 10 simulated patterns using a spatial grid of side=256. A logseries or uniform species abundance distribution were used, with 8,64 and 256 species. Two forms of generalized dimensions were estimated: DqSRS, from species rank surface $D_q^{SRS}$. DqSAD, estimated from species abundance distribution $D_q^{SAD}$. I use two spatial patterns: Regular, the species are distributed in vertical bands, Randomized the spatial distribution of species was randomized.}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{DqFit_Unif_64_256.png}
\caption{Linear fit of $Z_q(\epsilon)$ to estimate generalized dimensions ($D_q$) from simple species spatial patterns with 64 species and a uniform abundance distribution. The spatial grid has a side=256 sites and two different spatial patterns: Regular, a regular spatial pattern with species distributed in vertical bands of equal width. Randomized, the positions of species in the grid are randomized. Two kinds of generalized dimension were estimated: DqSRS corresponds the fit of $D_q^{SRS}$ (see text) and DqSAD is the fit from the estimation of $D_q^{SAD}$ (see text). $Z_q(\epsilon)$ corresponds to the partition function calculated for a box with side $\epsilon$, $q$ is the moment order.} 
\end{figure}. 


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{spat_NeuLog64_256.png}
\caption{Spatial patterns generated with a spatial neutral/hierarchical model. \textit{Rho} is the parameter that determines the degree of neutrality. When this parameter is 0 the model is completely neutral and there is no competitive replacement of species. When $\rho$ is 1 competitive superior species always replaces inferior ones and the model is completely hierarchical. \textit{Species} is the number of species actually present plot. The simulations use a metacommunity with a logseries abundance distribution with 86 species and a simulation grid side=256.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{RAD_NeuLog_spMeta_256.png}
\caption{Rank abundance diagrams of communities patterns generated with a spatial neutral/hierarchical model with different number of species in the metacommunity (labeled in each subfigure). \textit{$\rho$} is the parameter that determines the degree of neutrality. When this parameter is 0 the model is completely neutral and there is no competitive replacement of species. When $\rho$ is 1 competitive superior species always replaces inferior ones and the model is completely hierarchical. The simulations use a metacommunity with a logseries abundance distribution with 11, 86 and 341 species and a simulation grid side=256, the other parameters used were MortalityRate=0.2, DispersalDistance=0.4 (2.5 grid units), ColonizationRate=0.001. The ranks were calculated with averages of species densities over 50 simulations.}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{DqFit_NeuLog64_256.png}
\caption{Linear fit of $Z_q(\epsilon)$ to estimate generalized dimensions ($D_q$) from spatial patterns generated with a neutral/hierarchical model. \textit{Rho} is the parameter that determines the degree of neutrality. When this parameter is 0 the model is completely neutral and there is no competitive replacement of species. When $\rho$ is 1 competitive superior species always replaces inferior ones and the model is completely hierarchical. \textit{Species} is the number of species actually present plot. The simulations use a metacommunity with a logseries abundance distribution with 86 species and a simulation grid side=256 sites.}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{Dq_NeuLog_NumSp_256.png}
\caption{Generalized dimension spectra $D_q$ of spatial patterns generated with a spatial neutral/hierarchical model. \textit{Replacement} is the parameter $\rho$, that determines the degree of neutrality. When this parameter is 0 the model is completely neutral and there is no competitive replacement of species. When $\rho$ is 1 competitive superior species always replaces inferior ones and the model is completely hierarchical. The points are means and vertical lines are standard deviation of 50 simulated patterns. Simulations use a metacommunity with a logseries abundance distribution with 11,86 and 341 species. The simulation grid side is 256, and the other parameters are given in the main text.}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=6.5in]{powADq10_RplR_256.png}
\caption{Power of the Anderson-Darling test for the hypothesis of differences between simulated neutral/hierarchical communities. The test uses generalized dimensions curves calculated from SAD ($D_q^{SAD}$), generalized dimensions calculated from the species rank surfaces ($D_q^{SRS}$) and the species abundance distributions (SAD). The compared communities differ only in parameter $\rho$ (across panels) that determines the degree of neutrality/hierarchy. The number of comparisons to calculate the frequency is 2500 in all cases. Simulations use a metacommunity with a logseries abundance distribution with 11,86 and 341 species; a grid side of 256 sites, the other parameters are given in the main text.}
\end{figure}