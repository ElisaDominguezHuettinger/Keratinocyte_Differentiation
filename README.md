
(I) Boolean analysis - Written in R studio 1.4.1717

- Boolean network: **Boolean_Model_for_Keratinocyte_Differentiation_19May2022.txt**
- Reproduce figure 2: **Figure_2_Boolean_Analysis.r**

(II) ODE analysis - Main text - written in Matlab R2022a

(i) Functions called by the files that reproduce the figures:

- ODE model as a function of time, state and parameters: **Keratinocyte_Differentiation_ODE_Model_Integral_Np63.m**
- Retrieves **stable** steady states of the model as a function of parameters: 
**Keratinocyte_Differentiation_ODE_Model_SS_Int.m**" 
- Retrieves **unstable** steady states of the model as a function of parameters: 
**Sepparatrix_Keratinocyte_Differentiation_ODE_Model_SS_Int.m**
- Cost function (model vs. data), as function of parameters and experimental data: **CostFunction_Keratinocytes_Sep2023_INT.m**

(ii) Files to reproduce the figures:

- Fig 3B: **Figure3B_Optimization.m**
- Fig 3C: **Figure3C_Validation_ReversibilityAssays.m**

- Fig 4:  **Figure4_BifurcationAnalysis.m**

- Fig 5A and 5C: **Figure5AandC.m**
- Fig 5B **Figure5B.m**
- Fi  5D **Figure5D.m** 

(III) ODE analysis - Supplementary material - written in Matlab R2022a

- In-line figures in supplementary section 3:  **Fig_Sup_a1.m** and  **Curve_Sketching_Nullclines.m**
- Figure S3: ** "Fig_Sup_a2.m**; for an extended version: **Plot_Normalized_Toufighi_data.m**, reads data from:  "Names_of_genes.mat", "selected_genes_normalized.mat", "selected_genes_raw.mat", and "times_in_hours.mat".
- Figure S5: **Fig_Sup_a3.m** 

(IV) Raw data extracted from papers (including data for model calibration and validation)  - Excel

-**Supplementary_Table_1_Experimental_data_TerminalDifferentiationMarkers.xlsx**
Database of the dynamical responses of terminal differentiation marker expression in response to a step incrase in calcium.


(V) Optimal parameters of the ODE model

- **Parameter values.txt**
