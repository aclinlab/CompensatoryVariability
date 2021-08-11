## These scripts can be used to reproduce the results of Figure 2,3,S2 and S1. 

Running **VarDegradesPerf.m** will run its required functions and generate the following in your current MATLAB directory:

1. **Testing accuracies for the 8 models** (Homogeneous model with all 3 network parameters fixed, to the Random model; with all params are varying), see Fig2 in the paper, these .mat data will be named: "*test_p_ra***.mat", 
   .mat file name | corresponding model
   ------------ | -------------
   'test_p_ra' | the Random model
   'test_p_raH | Homoegenous W and N,variable theta
   'test_p_ra_Fixedtheta'| Random variable N and W, fixed theta
   'test_p_raH_FixedTheta' | Homogeneous W,N and theta. **all params fixed**
   'test_p_ra_varW_FixedN_FixedTheta' | Variable W, fixed N &theta
   'test_p_ra_varW_FixedN'| Variable W and theta, fixed N. 
   'test_p_ra_varN_FixedW_FixedTheta'| Variable N, fixed W and theta
   'test_p_ra_varN_FixedW' | Variable N and theta, fixed W.

   These matrices calculate the models' accuracies as a function of 4 variables: they are 5D matrices with the following dimensions:
   1. number of input odors (N)
   2. random network trials (each fly with different input connectivity) 
   3. noise in the input
   4. learning rate for updating the KC-MBON output weights
   5. determinancy in the decision making (c) 

   The bars in **Fig 2B** are at the peak performance conditions, (N = 100, noise in the input = 1, peak learning rates for each model, and c = 10).      The performance trends in **Fig2 C1-C4** were plotted as a function of each of these 4 variables while holding the other 3 constant. 
                                          
1.  **Dimensionality** (Fig 3E) and **lifetime sparsity** (Fig 4) for the Homogenous and Random models are in *dim_HF, dim_S*,*std_H_LTSpar* and *std_S_LTSpar*, respectively.
1. **DBI between odor pairs and DBI between 'good' and 'bad' odor classes** for the Homogenous and Random models. 
1. **Angular distance between odor pairs and Angular Distance between 'good' and 'bad' odor classes** for the Homogenous and Random models.

Finally, to replicate the results in **Fig S1** which uses the original Hallem-Carlson data and not fictitious odors derived from it, in the **VarDegradesPerf.m** code one needs to:

1. change line 20 &21 to: **modss=1; mulOd=110;**
2. By commenting line **60** and uncommenting **63**. 

