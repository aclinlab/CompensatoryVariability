## These scripts can be used to reproduce the results of Figure 2, 3, S1 and S2. 

Running **VarDegradesPerf.m** will run its required functions and generate, in your current MATLAB directory, **Testing accuracies for the 8 models** (Homogeneous model with all 3 network parameters fixed, to the Random model; with all params are varying), see Fig2 in the paper, these .mat data will be named: "*test_p_ra***.mat", 
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

**Fig 2E** is produced by `Perf_AtNonSparseCodingRegime_Random_HomogModels.m` and `Perf_AtSparseCodingRegime_APLNotZero_Random_HomogModels.m`, which also produced the metrics in **Fig 3** and **Fig S2** (angular distance, dimensionality, lifetime sparseness, valence specificity)

Finally, to replicate the results in **Fig S1** which uses the original Hallem-Carlson data and not fictitious odors derived from it, use `varDegradesPerf_S1_HallemOslenInp.m`
