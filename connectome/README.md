### How to reproduce the analysis from Fig. 6 and S5 of Abdelrahman et al., 2021

1. Download the raw connectome data: Run `getAllKCskelsInputsROIsv1.1.py` using Python 3 to download 1927 KCs (the skeletons, their synaptic inputs and their synaptic outputs). (You need to install neuprint-python first: https://github.com/connectome-neuprint/neuprint-python.) NB: this relies on ROI annotations which are from https://storage.cloud.google.com/hemibrain/v1.1/hemibrain-v1.1-primary-roi-segmentation.tar.gz. Creates the following files (you might want to modify the Python script to put these files in their own folder):
    * `KCids_v1.1.csv`
    * `KCskel0.csv` ... `KCskel1926.csv`
    * `KCinputs0.csv` ... `KCinputs1926.csv`
    * `KCoutputs0.csv` ... `KCsoutputs1926.csv`

2. Run `createAllKCs.m`. This reads in all the skeleton and input/output .csv files and creates 1927 `KCskel` objects which will be saved individually as `kc0.mat` - `kc1926.mat`.

3. Run `measureKCfeatures.m`. NB: this excludes KCs with truncated skeletons (stored in `KCsToExclude.csv`).

4. Run `sortSubtypesAndDisplayCorrelation.m`. This generates graphs for Fig. 6G-I and Fig. S5.

5. Run `plotManyCorrelations.m`. This generates grids for Fig. 6K and calculates Holm-Bonferroni-corrected p-values for the correlations.

NB: the file `neurSkel.m` is slightly modified from the code from https://github.com/aclinlab/amin-et-al-2020 to store the ROIs of skeleton nodes. The class `KCskel` extends the class `neurSkel`.