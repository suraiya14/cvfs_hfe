**If you find our web application (https://shiny.tricities.wsu.edu/bacteriocin-prediction/) useful, please consider citing the following papers:**

• Akhter, S. and Miller, J.H., BPAGS: A web application for bacteriocin prediction via feature evaluation using alternating decision tree, genetic algorithm, and linear support vector classifier. Frontiers in Bioinformatics, 3, p.1284705.

• Akhter, S., & Miller, J. H., BaPreS: a software tool for predicting bacteriocins using an optimal set of features. BMC bioinformatics, 24(1), 313.

Scripts and Datasets:

feature_extract.R – Used to extract features from protein sequences.

pearson's_corr.R - To perform Pearson correlation analysis.

shap.R - To measure each feature contribution.

positiveBac.faa and negativeBac.faa - Bacteriocin (positive) and non-bacteriocin (negative) FASTA sequences.

Below the link is a Python script where I discretized data into categorical features:

https://colab.research.google.com/drive/1CGNsH4fa3TgtCyd2-5U4v8ZgSQGD7-2j#scrollTo=vH0VWd3Zs_qd&line=4&uniqifier=1.

cvfe - To implement cross-validation feature evaluation, please see the below link: https://github.com/mingren0130/CVFS_code.

hfe - To implement hypergraph-based feature evaluation, please see the below link: https://github.com/hypper-team/hypper.

Pearsoncorrelation.csv - CVS file for Pearson's correlation coefficient.

positiveBac.faa and negativeBac.faa - Bacteriocin (positive) and non-bacteriocin (negative) FASTA sequences.
