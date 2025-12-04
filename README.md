## Optimizing Parametric Survival Model Selection: An Analysis Using Constructed Individual Patient Data Extracted from Oncology Trials From 2014 to 2023
<br>

### Summary

PS: Codes in this repository are used to reproduce this analysis.<br>
Original data and colour figures are also provided in this repository.<br>
<br>
Any questions or technical consulting:<br>
Email: 3220040596@stu.cpu.edu.cn <br>
<br>

### Description 

#### *1_Control.R (in Code folder):*
This script is used to read the original Individual Patient Data (IPD), call the main survival extrapolation function, and perform survival model fitting and extrapolation separately for the intervention and control groups. The resulting outputs are saved as external files.<br>

#### *2_Extrapolation function.R (in Code folder):*
This is the core script of the analysis. It is based on the methodology proposed by Kearns et al. (Medical Decision Making, 2019; 39(7):867–878), 
and implements both standard parametric and flexible parametric survival models for extrapolation. 
The script defines the main survival extrapolation function (Surv_analysis), which fits and extrapolates long-term survival curves using a range of models, 
including Exponential, Weibull, Gamma, Log-normal, Gompertz, Log-logistic, Generalized Gamma, Fractional Polynomial (FP), Restricted Cubic Spline (RCS), Royston-Parmar Spline (RP), 
Generalized Additive Model (GAM), and Parametric Mixture Models (PMM). The goodness-of-fit metrics, including Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC), 
are automatically calculated for model comparison.<br>

#### *3_Generate results.R (in Code folder):*
This script calculates various extrapolation performance metrics based on the extrapolated survival curves, including Restricted Mean Survival Time (RMST), 
Mean RMST (MRMST), Restricted Mean Survival Difference (RMSD), Mean Absolute Error (MAE, also referred to as ARMST), and Mean Absolute Percentage Error (MAPE). 
In this study, only RMSD and MAE were used for evaluating extrapolation error in the subsequent modeling; the other metrics were calculated but not included in further analysis.
Note: The results correspond to the models with the lowest AIC/BIC values. However, visual inspection is recommended to assess whether the extrapolated curves are clinically plausible. 
If not, the selected model should be excluded from the candidate model list, and the code rerun accordingly.<br>

#### *Example IPD.xls (in Code folder):*
This file contains example IPD used in the analysis, derived from the ATTRACTION-3 trial, including both pre-update and post-update overall survival (OS) data. Six datasets are presented:<br>
* Treatment group, pre-update IPD (with Minimum Number-at-Risk considered)<br>
* Treatment group, pre-update IPD (without Minimum Number-at-Risk)<br>
* Treatment group, post-update IPD (with Minimum Number-at-Risk)<br>
* Control group, pre-update IPD (with Minimum Number-at-Risk)<br>
* Control group, pre-update IPD (without Minimum Number-at-Risk)<br>
* Control group, post-update IPD (with Minimum Number-at-Risk)<br>

#### *Data.rar:*
This zip file contained all the reconstructed IPD data used in this study.

