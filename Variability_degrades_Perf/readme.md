# These scripts are intended to reproduce the results of Figure 2,3,4 and S1. 
# By running **VarDegradesPerf.m**, it will run its required functions dependencies and will generate the following in the MATLAB local current directory:

1. **testing accuracies for the 8 models** (Homogeneous model with all 3 network parameters fixed, to the Random model; with all params are varrying), see Fig2 in the paper, these .mat data will be named: "*test_p_ra**.mat", 
   1. 'test_p_ra'--> the Random model
   1. 'test_p_raH --> Homoegenous W and N,variable theta
   1. 'test_p_ra_Fixedtheta'-->Random variable N and W, fixed theta
   1. 'test_p_raH_FixedTheta' --> Homogeneous W,N and theta. **all params fixed**
   1. 'test_p_ra_varW_FixedN_FixedTheta' --> Variable W, fixed N &theta
   1. 'test_p_ra_varW_FixedN'--> Variable W and theta, fixed N. 
   1. 'test_p_ra_varN_FixedW_FixedTheta'--> Variable N, fixed W and theta
   1. 'test_p_ra_varN_FixedW'--> Variable N and theta, fixed W.

   These matrices calculate the models accuracies function of 4 variables: they are 5D matrices, {number of input odors **N** x random network trial (**fly**) x        **noise in the input** x **learning rate to update the KC-MBON output weights** x **determinancy in the decision making (C)**} 

   The bars in **Fig2 B** are at the peak performance conditions, (**N** =100, **noise in the input** =0.1, **peak learning rates for each model**, and **C=10**).      The performances trends were plotted though as function of each of these 4 variables while holding the other 3 constant in **Fig2 C1-C4.** 
                                          
1.  **dimensionality** and **Lifetime sparsity** for the Homogenous and Random models, *dim_HF, dim_S*,*std_H_LTSpar* and *std_S_LTSpar* giving the data in **Fig 3E and Fig4**, respectively.
1. **DBI between odors pairs, DBI between 'good' and 'bad' odor classes** for the Homogenous and Random models. 
1. **Angular distance between odors pairs, Angular Distance between 'good' and 'bad' odor classes** for the Homogenous and Random models.

## Finally, for replicating the results in **FigS1** which uses the original Hallem-Carlson data and not fictitious odors derived from it : in the **VarDegradesPerf.m** code one needs to:
1. change line 20 &21 to: **modss=1; mulOd=110;**
2. By commenting line **60** and uncommenting **63**. 

