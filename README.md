**If you find our web application (https://shiny.tricities.wsu.edu/bacteriocin-prediction/) useful, please consider citing the following papers:**

• Akhter, S. and Miller, J.H., BPAGS: A web application for bacteriocin prediction via feature evaluation using alternating decision tree, genetic algorithm, and linear support vector classifier. Frontiers in Bioinformatics, 3, p.1284705.

• Akhter, S., & Miller, J. H., BaPreS: a software tool for predicting bacteriocins using an optimal set of features. BMC bioinformatics, 24(1), 313.

**Scripts and Datasets:**

**feature_extract.R** – Used to extract features from protein sequences.

**pearsons_corr.R** - To perform Pearson correlation analysis.

**shap.R** - To measure each feature contribution.

**positiveBac.faa and negativeBac.faa** - Bacteriocin (positive) and non-bacteriocin (negative) FASTA sequences.

**cvfs.py** -To implement cross-validation feature evaluation, please download the cvfs.csv and cvfs.py files and save them to your desired location. Then, open the command prompt. Next, copy the path where you saved the cvfs.py and cvfs.csv files and paste it into the command prompt. Finally, run the following command:cvfs.py -i cvfs.csv -o C2P5P0.4 -c 2 -e 5 -p 0.4

**hfe.ipynb** - To implement hypergraph-based feature evaluation, please use the jupyter notebook and this file.

**pearsoncorrelation.csv** - CVS file for Pearson's correlation coefficient.


