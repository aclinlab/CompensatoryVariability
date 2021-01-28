These scripts are intended to reproduce the results of Figure 2,3,4 and S1. 
By running **VarDegradesPerf.m**, it will run its required functions dependencies and will generate the following in the MATLAB local current directory:

(a) **testing accuracies for the 8 models** (Homogeneous model with all 3 network parameters fixed, to the Random model; with all params are varrying), see Fig2 in the paper, these .mat data will be named: "*test_p_ra**.mat", 
                                          'test_p_ra'--> the Random model
                                          'test_p_raH --> Homoegenous W and N,variable theta
                                          'test_p_ra_Fixedtheta'-->Random variable N and W, fixed theta
                                          'test_p_raH_FixedTheta' --> Homogeneous W,N and theta. **all params fixed**
                                          'test_p_ra_varW_FixedN_FixedTheta' --> Variable W, fixed N &theta
                                          'test_p_ra_varW_FixedN'--> Variable W and theta, fixed N. 
                                          'test_p_ra_varN_FixedW_FixedTheta'--> Variable N, fixed W and theta
                                          'test_p_ra_varN_FixedW'--> Variable N and theta, fixed W.
       these matrices calculate the models accuracies function of 4 parameters, they are 5D, [number of input odors x random network trial (fly) x noise in the input x learning rate in the Output weights x determinancy in the decision making] 
       The bars in **Fig2 B** are the average at the peak performance conditions, (Number of input odors =100, noise in the input =0.1, peak learning rates for each model, and high determinancy in the decision making (C=10)). The performances trends were plotted though as function of each of these 4 variables while holding the other 3 constant in **Fig2 C1-C4.** 
                                          
(b) **dimensionality** and **Lifetime sparsity** for the Homogenous and Random models, *dim_HF, dim_S*,*std_H_LTSpar* and *std_S_LTSpar* giving the data in **Fig 3E and Fig4**, respectively.
(c) **DBI between odors pairs, DBI between 'good' and 'bad' odor classes** for the Homogenous and Random models. 
(d) **Angular distance between odors pairs, Angular Distance between 'good' and 'bad' odor classes** for the Homogenous and Random models.




