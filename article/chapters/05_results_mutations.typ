== 3.4 Response Prediction from Cell Fractions
Predicting patient response from raw deconvolution sub-cluster fractions using Random Forest models yielded moderate predictive accuracy in large cohorts (e.g. `Combined-Melanoma` ROC AUC = *0.707*, @fig-deconv-pred). In smaller cohorts with elevated noise, embedding *Univariate ANOVA Feature Selection* directly within the cross-validation loop substantially improved performance and mitigated variance, boosting the ROC AUC in breast cancer (`Anders`, $N=31$) from 0.518 to *0.735* and the Matthews Correlation Coefficient (MCC) from -0.142 to *0.392* (@fig-supp-s11).

#figure(
  image("../figures/deconvolution/deconv_prediction_advanced_grid_kmeans_subcluster_res_0.5.png", width: 90%),
  caption: [Performance metrics (Mean ± SD across 5 random seeds) of Random Forest response prediction using fine-grained deconvolution fractions.]
) <fig-deconv-pred>

== 3.5 Somatic Mutations & _TGM6_ Alterations
Training Random Forest classifiers on binary somatic mutation indicators achieved high predictive performance across tumor types (e.g., cross-validated ROC AUC = *0.738* in Rosenberg bladder, *0.726* in Riaz melanoma). In comparison, artificial neural networks (Adaline and MLP-4) exhibited severe overfitting when trained on high-dimensional genomic feature sets (validation AUCs plummeting to 0.30–0.40), but achieved competitive stability when restricted to top univariate-filtered genes (@fig-supp-s12).

Screening specific candidate mutations across melanoma cohorts identified a marked clinical association between *somatic mutations in _TGM6_* (transglutaminase 6) and checkpoint response:
- *Response Rate*: 13 out of 14 _TGM6_-mutated patients responded (*92.9%*), compared with 66 responders out of 205 wildtype patients (*32.2%*).
- *Statistical Association*: Two-sided Fisher's Exact Test yielded $p = 8.719 times 10^(-6)$ with an Odds Ratio of *27.38* (95% CI: 3.52 – 212.8).
- *TMB Linkage*: _TGM6_-mutant tumors displayed significantly higher overall TMB: mean 68.63 mut/Mb versus 15.33 mut/Mb in wildtype tumors (Mann-Whitney U $p = 2.638 times 10^(-5)$).
- *Independent Prognostic Contribution*: Adding _TGM6_ mutation status to a TMB-only model improved cross-validated ROC AUC from 0.5867 to *0.6045*.
- *mRNA Decoupling*: Importantly, _TGM6_ transcript abundance was unassociated with response (Mann-Whitney U $p = 0.2060$), confirming that the biomarker value resides specifically in the genomic alteration rather than steady-state gene expression.
