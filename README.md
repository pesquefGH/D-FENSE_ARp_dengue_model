# D-FENSE_ARp_dengue_model
Repository of the ARp model for the D-FENSE team 
1) D-FENSE (Dynamics for Epidemic Surveillance and Evaluation)

Paulo Antonio Andrade Esquef (National Laboratory for Scientific Reearch)

Add other tem members 


2) Repository Structure


Aggregated_Data: raw data (CSV spreadsheets) of dengue cases per state used model estimation. 

DFense_ARp: dengue cases prediction based on high-order AR(p) model

  |_ validation1:  material related to validation 1 challenge

      |_ matlab: Matlab scritps needed to run (run_batch_v1_predictor_ARp.m) the simulation and generate the CSV and PDF files, related to dengue case predictions for each state. CSV files are stored in planilhas and related plots (in PDF) are stored in plots  

      |_ planilhas: stores CSV files, one for each state, with predictions of dengue cases

      |_ plots: stores PDF files, one for each state, with 4 subplots related to predictions of dengue cases: median prediction, 50%, 80%, 90%, and 95% preduction intervals.


  |_ validation2:  material related to validation 2 challenge

      |_ matlab: Matlab scritps needed to run (run_batch_v2_predictor_ARp.m) the simulation and generate the CSV and PDF files, related to dengue case predictions for each state. CSV files are stored in planilhas and related plots (in PDF) are stored in plots  

      |_ planilhas: stores CSV files, one for each state, with predictions of dengue cases

      |_ plots: stores PDF files, one for each state, with 4 subplots related to predictions of dengue cases: median prediction, 50%, 80%, 90%, and 95% preduction intervals.

  |_ validation3:  material related to validation 2 challenge

      |_ matlab: Matlab scritps needed to run (run_batch_v3_predictor_ARp.m) the simulation and generate the CSV and PDF files, related to dengue case predictions for each state. CSV files are stored in planilhas and related plots (in PDF) are stored in plots  

      |_ planilhas: stores CSV files, one for each state, with predictions of dengue cases

      |_ plots: stores PDF files, one for each state, with 4 subplots related to predictions of dengue cases: median prediction, 50%, 80%, 90%, and 95% preduction intervals.



3) Libraries and Dependencies

matlab functions: readtable.m, buffer.m (Signal Processing Toolbox), armcov.m (Signal Processing Toolbox), filter.m (Signal Processing Toolbox), filter2.m (Signal Processing Toolbox), ssa_modPE.m (Singular Spectral Analysis - Smoothing Filter, included in the folder 'matlab').

 

4) Data and Variables

Only the time-series of raw number of dengue cases per state along epidemic weeks have been used. Data are
available from:
 
https://github.com/americocunhajr/D-FENSE/tree/main/DengueSprint2025_DataAggregated


5) Model Training

DFense_ARp: for each state (UF), the log2 mapping of time-series of raw dengue cases, in the defined range for each validation, have been used to estimate an AR(p), p=92 (experimentaly chosen), via function armcov.m. Initial conditions for the AR(p) model at epidemic week (EW) 25 of 2022/23/24 have been obtained by a simple scheme of inverse filtering of the time-series, followed by direct filtering of the modeling error. The modeling error sequence has been organized in a matrix with 52 columns, each matrix row containing a modeling error sequence related to one year. Assuming a zero-mean Gaussian White distribution for the modeling error ensemble, the standard deviation of a typical model excitation has been estimated. Then, a Monte Carlo simulation with 10000 runs has been carried out, to generate the dengue cases predictions: the AR(p) and initial conditions were fixed, only the model excitation have been drawn from a Gaussian distribution. Each of these model excitations have 79 samples, covering a forecast from EW 26 of a given year to EW 52 of the subsequent year. Then, the attained results have been mapped back to the original amplitude domain (via the inverse of the log2 function). From the set of these 10000 case predictions, the median, lower- and upper-bounds of the 50%, 80%, 0%, 90%, and 95% prediction intervals are calculated. Finally, the resulting curves are smoothed out via an SSA (Singular Spectral Analysis) filter and cropped out to be in the range from EW 41 of a given year to EW 40 of the subsequent year.       
    
     
6) References

None.

7) Data Usage Restriction

None. 

DFense_ARp: from the trained/estimated model (see section 5 above): we run a Monte Carlo simulation with 10000 runs to generate the dengue cases predictions: the AR(p) and initial conditions were fixed, only the model excitation have been drawn from a zero-mean Gaussian distribution, whose standard deviation has been estimated from the modeling error. Each of these artificially generated model excitations have 79 samples, covering a forecast range from EW 26 of a given year to EW 52 of the subsequent year. Then, the attained results have been mapped back to the original amplitude domain (via the inverse of the log2 function, 2^(predictions)). From the set of these 10000 case predictions, the median, lower- and upper-bounds of the 50%, 80%, 0%, 90%, and 95% prediction intervals have been calculated. Finally, the resulting curves are smoothed out via an SSA (Singular Spectral Analysis) filter and cropped out to be in the range from EW 41 of a given year to EW 40 of the subsequent year. 

8) Predictive Uncertainty

From the set of 10000 case predictions (for each state and each validation), we used matlab function prctile.m (percentiles of a sample) to obtain the median, as well as the lower- and upper bounds of 50%, 80%, 90%, and 95% prediction intervals. The median of the case predictions is the 50% percentile. The lower-bounds for the 50%, 80%, 90%, and 95% prediction intervals are, respectively, the 25%, 10%, 5%, and 2.5% percentiles. The upper-bounds for the 50%, 80%, 90%, and 95% prediction intervals are, respectively, the 75%, 90%, 95%, and 97.5% percentiles.

Summary:
median prediction: 50% percentile
50% prediction interval: from 25% percentile to 75% percentile  
80% prediction interval: from 10% percentile to 90% percentile
90% prediction interval: from 5% percentile to 95% percentile
95% prediction interval: from 2.5% percentile to 97.5% percentile
