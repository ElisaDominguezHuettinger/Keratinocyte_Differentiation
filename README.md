**README Keratinocyte differentiation motif**

**Manuscript**: History-dependent switch-like differentiation of keratinocytes in response to skin barrier damage
Authors: Elisa Domíguez-Hüttinger1,*, Eliezer Flores-Garza3, José Luis Caldú-Primo2, Harley Day3,
Abihail Roque Ramírez4, and Reiko J Tanaka3,+
*Affiliations:*
1Departamento de Biología Molecular y Biotecnología, Instituto de Investigaciones Biomédicas, Universidad Nacional Autónoma de México, Ciudad Universitaria, 04510, México, México. 
2Doctorado en Ciencias Biomédicas, Universidad Nacional Autónoma de México, Ciudad Universitaria, 04510, México, México.
3Department of Bioengineering, Imperial College London, South Kensington Campus, London SW7 2AZ
4Facultad de Ciencias, Universidad Nacional Autónoma de México, Ciudad Universitaria, 04510, México, México.
Correspondence: * elisa.dominguez@iibiomedicas.unam.mx; + r.tanaka@imperial.ac.uk


*The code is described in order of appearance of the main text*.

**Boolean analysis**

Figure S2

Results of the Boolean model of the regulatory network for keratinocyte differentiation. The model attractors under (A) basal and (B) high calcium conditions. (C) The size of the basins of attraction corresponding to the differentiated state under basal and high calcium conditions. (D) The model attractors under an asynchronous update regime coincide with the fixed-point attractors of the synchronous update regime.

•	Supplementary Equation 1 (Boolean model of keratinocyte differentiation) is provided in the .txt file: "Boolean_Model_for_Keratinocyte_Differentiation.txt"
•	The code to reproduce Supplementary Figure S2 was written in R studio 1.4.1717, using the BoolNet package: "Figure_S2_Boolean_Analysis.r" 

**ODE analysis**
Software: 
Matlab R2022a

*Functions called by the files that reproduce the figures:*

•	The keratinocyte differentiation motif, described as (Eq. 1) in the main text, is declared as the function Keratinocyte_Differentiation_ODE_Model_Integral_Np63.m.
•	The stable steady states of the model (Eq. 1) are computed as a function of model parameter values with: Keratinocyte_Differentiation_ODE_Model_SS_Int.m. 
•	This function calls the vpasolve function of Matlab to compute all the steady states, filters out the steady states that are either negative or complex, evaluates the stability of the remaining steady states by substituting the corresponding state variables in the Jacobian matrix, and filters out the unstable steady states to retrieve only the positive, real and stable steady states as a function of model parameters.
•	The unstable steady states are computed in Sepparatrix_Keratinocyte_Differentiation_ODE_Model_SS_Int.m, following the same procedure as for the stable steady states, but this time keeping only the unstable ones. The function is called “Sepparatrix…” because we use these unstable steady states to compute the separatrix that separates the basins of attraction. 
•	The cost function measuring the difference between the model simulations and the experimental data is given in CostFunction_Keratinocytes_ INT.m. It is a function of model parameters and experimental data.

*Files to reproduce the figures of the manuscript (main text):*

 Figure 2B: 

The parameter optimization of the keratinocyte differentiation motif is performed by running Figure2B_Optimization.m. This function calls a global optimization algorithm to find the optimal parameter values that minimize the difference between the model simulations and the experimental data.

The resulting optimal parameters are stored in the txt file: Parameter values.txt

Figure 2C

The validation of the keratinocyte differentiation motif is performed by running the function: 
Figure2C_Validation_ReversibilityAssays. This function numerically integrates the keratinocyte differentiation motif with the optimal model parameters and plots the resulting dynamic trajectories alongside the validation experimental dataset.
Figure 3: 
 
Figure3_BifurcationAnalysis.m. produces the phase-plane representation of the nullclines for different Calcium concentrations by computing a closed expression which we declare as an anonymous function. 
The same code also produces the corresponding bifurcation diagrams, by calling vpasolve function to compute the steady states for varying calcium concentrations and evaluate their linear stability through the eigenvalues of the corresponding Jacobian matrix.

Figure 4A:

Bifurcation diagrams for the 3 state variables (Stat3, Np63 and TDM) with three levels of NFkB (0, 0.1 and 0.25). The bifurcation diagrams are computed following the same computational procedure as in Figure 3.

Figure4A_Bifurcation_diagrams_NFkB.m 
    
Figure 4B

The threshold calcium concentration C+ that switches TDM on, as a function of NFkB (variation of NFkB: 0:0.01:1) is computed with the code: Figure4B.m, that calculates bifurcation diagrams as in figure 4A, stores and finally plots the threshold values C+.



*Files to reproduce the figures of the manuscript (supplementary material):*

Section S5: Analytical derivation of necessary conditions for bistability of the keratinocyte differentiation motif


The step-by-step symbolic analysis of the nullclines through curve sketching can be retrieved from Section_S5_Nullcline_Analysis.m and Section_S5_Nullcline_Analysis_Phase_Space.m

Figure S3:

The dynamic, quantitative and high-throughput experimental data from Toufighi et al. is plotted with Figure_S3_Experimental_Data.m

Figure S4:

Modelled dynamics of Np63 show a steady decrease in its expression which is qualitatively consistent with the observed reduction in Np63 expression measured in the two independent experiments of Toufighi et al. and Lena et al. Figure_S4_Validation_Np63.m

Figure S5:

The model simulations alongside a plot of the experimental data from Borowiec et al 2013 can be performed with Figure_S5_Model_Validation_decay.m


Figure S6

Code to retrieve the symbolic nullcline analysis that shows that adding an inhibition of Stat3 by Np63 as a Np63-mediated non-linear degradation term of Stat3 changes the Stat3 nullcline from a linear to a saturated curve that converges to an asymptote that is inversely proportional to the Np63-dependent degradation rate of Stat3:  Figure_S6A_Kurvendiscussion_Nullclines_With_Negative_Np63_on_Stat3.m

Code showing that the fit to the experimental data is also invariant to adding a Np63-dependent degradation rate of Stat3 to our original model: Figure_S6B_Optimization_With_Negative_Np63_on_Stat3.m

