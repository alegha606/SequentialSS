
/* ************************************************************************** */
/*        0. IMPORTING IN DATA AND CREATING AN EMPTY POSTFILE DATASET         */
/* ************************************************************************** */


clear all
set more off
set maxvar 5000

/* 0.0 Set working directory */

cd "C:\Users\User XYZ\Stata Results Tables_Unpenalised Log Reg + Bootstrap Uniform Shrinkage\"

/* 0.1 Setting the seed before beginning to allow reproducibility of results */

set seed 1235	

/* 0.2 Importing in the original MIMIC dataset */

use "mimic_cohort_analysis_ethnicity", clear

* Note: the MIMIC dataset used in this study is freely available and relates to the following paper: 
* Johnson AEW, Pollard TJ, Shen L, Lehman L-wH, Feng M, Ghassemi M, et al. MIMIC-III, a freely accessible critical care database. Scientific Data. 2016;3(1):160035.

/* 0.3 Keep only the outcome (development of acute kidney injury) and the 6 predictors of interest */

keep AKI bicarbonate_mean creatinine_mean hemoglobin_mean bun_mean potassium_mean sysbp_mean

/* 0.4 Create a uniform random variable with the same number of observations as 
the number of patients in the MIMIC dataset (this variable will be used to 
randomly select a new batch of 'recruited' patients from the MIMIC dataset) */

generate random = runiform()

/* 0.41 Sort by this new random variable and save the dataset */

sort random

drop random
 
save "mimic_reducedvars", replace

/* 0.5 Create a dataset to post the final summary output results from the 
sequential process into (see end of code).

- n will be stored (100, 200, ... 3000), as well as each performance measure
of interest (e.g. bootstrap corrected calibration slope, optimism in AUC, 
mean 95% uncertainty interval width etc) */ 

* 0.51 Naming the postfile dataset

tempname FinalDataset

/* 0.52 Setting up the structure of an initially empty dataset called 'FinalDataset' 
by defining which vars it contains */

postfile `FinalDataset'  n ///
UIwidth_mean UIwidth_pctile2dot5 UIwidth_pctile97dot5 ///
Delta_mean Delta_pctile2dot5 Delta_pctile97dot5 ///
MAPE_mean MAPE_pctile2dot5 MAPE_pctile97dot5 ///
probmisclass_mean probmisclass_pctile2dot5 probmisclass_pctile97dot5 ///
avg_EVPI ///
using "FinalDataset", replace 

/* ************************************************************************** */
/*         1. START OF FIRST LOOP: DIFFERENT SAMPLE SIZES 'N'                 */
/* ************************************************************************** */

/* 1.0 Now adding a loop for different sample sizes. 

The idea is that we want to mimic sequential recruitment of batches of patients, 
starting from an initial dataset of n=100, then increasing in increments of 
100 each time, and stopping at a final SS of n=3000. 

Note: for convenience (to speed up computational time), it's easier to do this 
process in reverse (starting with a randomly selected dataset of n=3000, and 
then removing 100 patients at a time, until we end up with just n=100 in the last 
iteration of the loop). 

Also note: due to model convergence issues at n = 100, the learning curves were 
started from Ninitial = 200

*/

forvalues j = 3000(-100)200 {

/* 1.00 Reloading the dataset called 'mimic_reducedvars' that contains 
the current pool of patients that could potentially be randomly selected from 
the MIMIC data for this particular iteration of the sample size loop */
	
use "mimic_reducedvars", replace


/* 1.001 As 201 EVPI values will be created per sample size iteration 
(see further details below) in this unpenalised logistic regression 
+ bootstrap uniform shrinkage modelling strategy, these 201 EVPI values are 
now stored iteratively in a new dataset that is created in this section of the 
code. Then, an average EVPI (from these 201 EVPI values) will later be calculated 
for a given sample size.

(Full details of the mechanism used in this bootstrap uniform shrinkage modelling 
strategy: for a given iteration of the learning curve generation process, 
a development sample of fixed size is chosen and used to develop the initial model.
Then, Harrell's bootstrap approach is applied to produce a bootstrap-corrected calibration 
slope, which is used as a shrinkage factor to produce a bootstrap penalised 'original' model.

To generate instability risk estimates for a given individual, B=200 
bootstrap samples are drawn randomly (with replacement) from the original 
development data (to produce B=200 'new' development datasets), 
then Harrell's bootstrap is applied to produce a penalised model for each of 
these B=200 'new' development datasets, and finally these penalised models are 
applied back to the original development data to give 200 instability risk 
estimates for each individual.

Thus, for the 'original' bootstrap penalised model, as well as the B=200 
penalised models based on the B=200 'new' bootstrap development datasets, a 
separate EVPI value will be calculated each time. We now store this iteratively 
with the following code...) */

/* 1.001 Naming the postfile dataset */

tempname BSShrinkageFactor_SS`j'

/* 1.002 Setting up the structure of an initially empty dataset called 
'Temp_BSShrinkageFactor_`i'' by defining which vars it contains */

postfile `BSShrinkageFactor_SS`j'' EVPI_temp /// 
using "BSShrinkageFactor_SS`j'", replace 

/* 1.01 Creating a local value for n, so we can 'post' it to the final output 
dataset later on */

local n = `j'

/* ************************************************************************** */
/* 1.1 RANDOMLY SELECTING A BATCH OF NEW PATIENTS FOR (OR TO REMOVE FROM)     */
/* INITIAL MODEL DEVELOPMENT DATASET                                          */ 
/* ************************************************************************** */

generate insample = _n <= `j'

keep if insample == 1

/* 1.11 Dropping the insample variable now as it is no longer needed within this loop */

drop insample

/* 1.12 Saving the updated dataset containing 'initial dev data', ready for 
use in analysis in next step */

save "mimic_reducedvars", replace

/* 1.13 Also saving another copy of this dataset for a particular 
sample size iteration, this 'MasterDevData' copy of the dataset will be updated 
throughout the rest of the code below, but 'mimic_reducedvars' will not.

Having a copy of the dataset here allows the processes below to run without 
issue */

save "MasterDevData", replace




/* ************************************************************************** */
/* 2. START OF SECOND (NEW) LOOP FOR INSTABILITY: 'B1' BOOTSTRAP REPS TO:     */
/*                                                                            */
/* i) PRODUCE FIRST A BOOTSTRAP SHRUNKEN PENALISED MODEL BASED ON THE ORIGINAL*/
/* DEV DATA (WHERE 'F'= 1 OF THE INSTABILITY LOOP), AND CRUCIALLY THIS THEN   */
/* GIVES US AN ORIGINAL RISK ESTIMATE PER INDIVIDUAL                          */
/*                                                                            */
/* ii) THEN, FOR 'F'= 2 TO 201 OF THE INSTABILITY LOOP, WE PRODUCE B1=200     */
/* RANDOM BOOTSTRAP SAMPLES FROM THE DEV DATA USED ABOVE, TO THEN BE ABLE TO  */
/* CREATE 200 SETS OF BOOTSTRAP SHRUNKEN/PENALISED MODELS, WHICH ALLOWS US TO */
/* WORK OUT B=200 (INSTABILITY) RISK ESTIMATES PER INDIVIDUAL                 */
/*                                                                            */
/* SEE COMMENTS IN CODE IN SECTION 1.001 FOR FURTHER DETAILS                  */
/*                                                                            */
/* ************************************************************************** */

qui forvalues f = 1/201 { // instability loop 

use "MasterDevData", replace

/* 2.01 Adding a visual counter of the sequential SS iteration (n), and 
stability loop iteration */

di as err "Currently on  n = " `j'
di as err "And on instability loop = " `f'

/* 2.1 Randomly selecting a bootstrap sample from the dev data 

Note: this random bootstrap sample is only selected for `f'>1, so that for the 
first iteration of the instability loop only (i.e. f=1), all of the dev data is 
used (i.e. no bootstrap sample is taken of the dev data when f=1). 

When performing instability calculations, pr_original (the probability estimate 
based on the original model) is then generated in the loop f=1, and compared to 
the 200 sets of probabilities generated from the 200 bootstrap shrunken models 
(when f = 2 to 201), based on 200 random bootstrap samples of the dev data).

*/

if `f'>1 {
			bsample
		}



		



/* ************************************************************************** */
/* 2.2 FITTING THE CHOSEN MODELLING STRATEGY ONTO CURRENT DEVELOPMENT DATASET,*/
/* AND STORING OVERALL MODEL PERFORMANCE MEASURES, AS WELL AS INDIVIDUAL RISK */
/* PREDICTIONS */
/* ************************************************************************** */

/* 2.20 Apply the initial model which is the base-case unpenalised model to the 
development dataset: 

We are fitting a logistic regression model to determine the probability of 
developing AKI (outcome) based on 6 (continuous) predictors 
(linear predictor-outcome relationships are assumed, and all 6 variables 
are forced into the model). */

* Changed max iterations from 300 to 500 to improve model convergence rate

qui logistic AKI bicarbonate_mean creatinine_mean hemoglobin_mean ///
bun_mean potassium_mean sysbp_mean, iterate(500)

/* 2.21 Calculate the predicted probabilities (let's name them 'pr_noshrinkage') 
and the value of the linear predictor (i.e. 'pr_noshrinkage' but on logit-p scale 
rather than p-scale, let's name: 'lp_noshrinkage') for each individual using 
the predict command. */

predict pr_noshrinkage, pr

gen lp_noshrinkage = logit(pr_noshrinkage)

/* 2.22 Also assessing the models apparent calibration via the calibration slope
 (will need this 'apparent' calibration later to be able to adjust the model 
 by the optimism in the c-slope from the third and final loop for an internal 
 validation bootstrap) */

qui logistic AKI lp_noshrinkage, coef
local cslope_orig = _b[lp_noshrinkage]

/* 2.3 Save predictions and original data in a new dataset */

save "BSshrinkage_originaldevdatawithmodelresults", replace












/* ************************************************************************** */
/* 3. START OF THIRD (FINAL AND INNERMOST) LOOP: 'B' BOOTSTRAP REPS TO:       */
/*                                                                            */
/* i) PERFORM INTERNAL VALIDATION OF THE (RANDOM BOOTSTRAP) DATA SAMPLE       */
/* SELECTED IN THE CURRENT ITERATION OF THE INSTABILITY LOOP `f', SO THAT CAN */
/* THEN GENERATE OPTIMISM ADJUSTED C-SLOPE, TO BE ABLE TO PRODUCE A BOOTSTRAP */
/* PENALISED MODEL FOR THE PARTICULAR INSTABILITY LOOP 'f' (RANDOM BOOTSTRAP) */
/* DATA SAMPLE MODEL)                                                         */
/* I.E. WE WANT TO APPLY BOOTSTRAP SHRINKAGE TO THE MODEL IN CODE SECTION 2.2 */
/*                                                                            */
/* NOTE: THIS FINAL BOOTSTRAP LOOP IS REPEATED (IN TOTAL IT IS RUN F=201 TIMES,*/
/* EACH TIME WITH B=200 REPS), SO THAT 'DOUBLE' BOOTSTRAPS ARE APPLIED IN THIS */
/* MODELLING STRATEGY CODE                                                    */ 
/*                                                                            */
/* ************************************************************************** */

/*
	Note, for the internal validation results (overall model performance measures),
	we store these in a matrix (B=200 bootstrap reps, with just one overall model 
	performance measure: calibration slope (don't need c-statistic, as only need c-slope 
	to apply bootstrap shrinkage later on), this is calculated with bootstrap 
	model first applied to bootstrap sample (i.e. 'apparent' performance) and 
	then bootstrap model applied to original model development dataset from 
	step 2 of code (i.e. 'test' performance)
	
	3x measures relating to EVPI for each of the B=200 bootstrap reps are also 
	stored.
	
	So this matrix contains 5x overall performance measures in total, and 200 rows)	
	
*/

matrix results = J(200,5,.)

qui forvalues k = 1/200 {

/* 3.01 Loading the current batch of development data (this could be the 
original development dataset, i.e. if f=1, or a bootstrap sample of the 
original development dataset, i.e. if f>1) and just keeping the 
outcome and 6 predictors (i.e. removing the predicted probabilities from the 
original model) */

use "BSshrinkage_originaldevdatawithmodelresults", clear
keep AKI bicarbonate_mean creatinine_mean hemoglobin_mean bun_mean potassium_mean sysbp_mean

* 3.02 Taking a bootstrap sample (random sample with replacement) of same size
bsample

* 3.03 Iteration indicator to help the user know the progress
nois _dots `k' 0

/* ************************************************************************** */

/* 3.1 Appling the initial model which is the base-case model to bootstrap dataset: 

Fitting a logistic regression model to determine probability of AKI based on 
6 continuous predictors (linear relationship for all vars with outcome, all 6 vars 
forced in, unpenalised regression). */

* changed max iterations from 300 to 500 to improve model convergence rate

logistic AKI bicarbonate_mean creatinine_mean hemoglobin_mean ///
bun_mean potassium_mean sysbp_mean, iterate(500)
		
/* 3.11 Calculate the predicted apparent bootstrap probabilities (let's name them 'pr') 
and the value of the linear predictor (i.e. 'pr' but on logit-p scale 
rather than p-scale, let's name: 'lp') for each individual using 
the predict command */

predict pr, pr

gen lp = logit(pr)

/* 3.12 Also assessing the models apparent calibration via the 
calibration slope.*/

logistic AKI lp, coef
matrix results[`k',1] = _b[lp]

/* ************************************************************************** */

/* 3.2 Now testing performance of the bootstrap model in original sample
	NB: to do this we must first fit the model to the bootstrap data again, so
	that Stata has the coefficients stored in memory */

* changed max iterations from 300 to 500 to improve model convergence rate
logistic AKI bicarbonate_mean creatinine_mean hemoglobin_mean ///
bun_mean potassium_mean sysbp_mean, iterate(500)

/* 3.21 Load the current batch of development data again */
use "BSshrinkage_originaldevdatawithmodelresults", clear

/* 3.22 Predict probabilities & linear predictor from the bootstrap model 
 applied to the original dataset */

predict pr, pr
gen lp = logit(pr)

/* 3.23 Also assessing the models test calibration by via the calibration slope.*/

logistic AKI lp, coef
matrix results[`k',2] = _b[lp]

/* 3.24 Now performing intermediary calculations to enable bootstrap-based computation 
of the Expected Value of Perfect Information (EVPI) */

/* 3.241 First estimating 'NB_All' for this current kth bootstrap rep 
(NB_All is the net benefit of administering a treatment (for example, closer 
monitoring of renal function) to all patients regardless of their underlying risk) 

For full mechanism and corresponding equations see: 
Sadatsafavi M, Yoon Lee T, Gustafson P. Uncertainty and the Value of Information in Risk Prediction Modeling. Medical Decision Making. 2022;42(5):661-71.
*/

gen NB_All_`k' = ((pr-(1-pr)*0.1)/(0.9))

qui su NB_All_`k', d

matrix results[`k',3] = r(mean)

/* 3.242 Then estimating 'NB_Model' (NB_Model is the net benefit of using the 
proposed model to decide which patients to treat)
 */

gen prorig_greater_z_indicator = 0
replace prorig_greater_z_indicator = 1 if pr_noshrinkage > 0.1 


gen NB_Model_`k' = (prorig_greater_z_indicator*((pr-(1-pr)*0.1)/(0.9)))

qui su NB_Model_`k', d

matrix results[`k',4] = r(mean)

/* 3.243 Finally estimating 'NB_Max' (NB_Max is the net benefit of using the 
`correct' model to decide which patients to treat) */

gen prbootstrap_greater_z_indicator = 0
replace prbootstrap_greater_z_indicator = 1 if pr > 0.1 

gen NB_Max_`k' = (prbootstrap_greater_z_indicator*((pr-(1-pr)*0.1)/(0.9)))

qui su NB_Max_`k', d

matrix results[`k',5] = r(mean)

}



/* ************************************************************************** */
/* 4. USING BOOTSTRAP RESULTS TO COMPUTE OVERALL POPULATION-LEVEL STABILITY   */
/* MEASURES:                                                                  */
/*                                                                            */
/* 1) OPTIMISM ADJUSTED C-SLOPE                                               */
/*                                                                            */
/* 2) EVPI                                                                    */
/*                                                                            */
/* ************************************************************************** */

/* 4.0 Rename columns of the matrix used to store bootstrap overall model performance results */

mat colnames results = cslope_boot cslope_test NB_All NB_Model NB_Max

/* 4.01 Clear and load these matrix results into a Stata dataset */

clear
svmat results, n(col)

/* 4.02 We now have a dataset of 200 estimates (from the B=200 bootstrap replications)
of apparent and test performance. So now we can calculate the optimism in the calibration slope */

gen optimism_slope = cslope_boot - cslope_test


/* 4.03 Summarise the optimism in the calibration slope
NB: this is simply the mean of the apparent bootstrap performance minus the test performance
(averaged across the B=200 bootstrap replications) */

/* Summarising and storing the mean optimism in the calibration slope (this will 
be needed to adjust the current development model by bootstrap uniform shrinkage 
in the next steps) */

qui su optimism_slope
local cslope_opt = r(mean)

/* 4.04 Now calculate the optimism adjusted performance of the original model by
	subtracting the optimism from the original models apparent performance in
	the original dataset, and store the optimism adjusted c-slope (as this is 
	one of my overall performance measure outputs I want to use as stopping rule 
	criteria; the other one is optimism in AUC, which I already have in above 
	step as 'cslope_opt') */
	
local cslope_bscorrected = `cslope_orig' - `cslope_opt'

di as err "Shrinkage factor = " `cslope_bscorrected'

/* 4.1 EVPI MEASURES */

/* 4.11 First calculating the ENBall, ENBModel, and ENBMax by taking the 
mean of the 200 bootstrap values for these measures.

As in Section 3.24, for full mechanism and corresponding equations see: 
Sadatsafavi M, Yoon Lee T, Gustafson P. Uncertainty and the Value of Information in Risk Prediction Modeling. Medical Decision Making. 2022;42(5):661-71.
 */

/* Calculating ENB_All 
(the expected value of the net benefit of treating all) */

qui su NB_All 
local ENB_All = r(mean)

/* Calculating ENB_Model 
(the expected value of the net benefit of using the proposed model to decide which patients to treat) */

qui su NB_Model 
local ENB_Model = r(mean)

/* Calculating ENB_Max 
(the expected value of the net benefit of using the 'correct' model to decide which patients to treat) */

qui su NB_Max 
local ENB_Max = r(mean)

/* 4.12 Finally, calculating EVPI.

This is the expected gain in net benefit of using the 'correct' model to decide which 
patients to treat instead of the proposed model (or a model which treats no one 
if this has a greater net benefit than the proposed model)

For further information see: Sadatsafavi M, Yoon Lee T, Gustafson P. Uncertainty and the Value of Information in Risk Prediction Modeling. Medical Decision Making. 2022;42(5):661-71.
*/

local EVPI_temp = `ENB_Max'-max(0,`ENB_Model',`ENB_All')

/* 4.2 Now 'posting' the 201 EVPI values into the dataset called 'Temp_BSShrinkageFactor_`i'' */

post `BSShrinkageFactor_SS`j'' (`EVPI_temp')

/* 4.3 Finally, we can now applying the bootstrap uniform shrinkage penalisation 
to the current unpenalised logistic regression model from code Section 2.2.

To do this, we use the bootstrap corrected c-slope as the shrinkage term 
to the original model (by original model, we mean in this context the base-case 
log-reg model fitted to the (bootstrap) sample of the dev data of the `f' 
iteration of the instability loop from code Section 2.2) */

/* 4.31 First need to fit the base-case model again to repeat step 2.2 (which applies 
the model to the (bootstrap) sample of the dev data) so that it is in Stata's memory */

use "BSshrinkage_originaldevdatawithmodelresults", clear

logistic AKI bicarbonate_mean creatinine_mean hemoglobin_mean ///
bun_mean potassium_mean sysbp_mean, iterate(500)

* 4.32 Now we store the constant/intercept term to update post shrinkage

local orig_cons = _b[_cons]

/* 4.33 Then we apply bootstrap shrinkage to our final model by multiplying the 
beta coefficient terms by the bootstrap corrected c-slope.

 However the linear predictor of a logistic model includes the constant/intercept term
 Which must be re-estimated post shrinkage of the beta terms.
 
 Therefore we must first remove this element from the linear predictor variable and
 then multiply the betas (lp) by the shrinkage factor */
 
gen lp_nocons = lp_noshrinkage - `orig_cons'

gen BSshrunken_lp_nocons = lp_nocons*`cslope_bscorrected'	

/* 4.34 Now to re-estimate a new intercept value allowing for the shrunken betas
 we fit a logistic model with the shrunken linear predictor set as an offset to 'fix' these
 shrunken beta coefficients.
 
 We then add this additional intercept term to the linear predictor to obtain
 our final shrunken model linear predictor */
 
logistic AKI, offset(BSshrunken_lp_nocons) coef
global new_cons = _b[_cons]


/* 4.35 Now we have the shrinkage term cslope_bscorrected for the fth instability 
loop bootstrap rep,
and the constant term (adjusted after shrinkage applied to beta coeffs) new_cons
for the fth instability bootstrap rep...

So we have everything we need to be able to calculate the individual shrunken 
coeffs terms (betas),
first just need to re-run the model again so it is in Statas memory to be able to 
extract original unshrunken betas, then we apply shrinkage to them 

*/

logistic AKI bicarbonate_mean creatinine_mean hemoglobin_mean ///
bun_mean potassium_mean sysbp_mean, iterate(500)

/* 4.36 Now finally can calculate and store the bootstrap shrunken beta terms from 
the bootstrap uniform shrinkage model, as we need them to be able to apply the 
bootstrap shrinkage model generated from the bootstrap sample  onto the original development data.

This step is needed as we don't have the model already stored in Stata memory,
so we ave to manually construct it. */

global bicarbonate_mean = _b[bicarbonate_mean]*`cslope_bscorrected'
global creatinine_mean = _b[creatinine_mean]*`cslope_bscorrected'
global hemoglobin_mean = _b[hemoglobin_mean]*`cslope_bscorrected'
global bun_mean = _b[bun_mean]*`cslope_bscorrected'
global potassium_mean = _b[potassium_mean]*`cslope_bscorrected'
global sysbp_mean = _b[sysbp_mean]*`cslope_bscorrected'

********************************************************************************
/* ************************************************************************** */

/* 4.4 Final step: applying the bootstrap uniform shrinkage model, to the actual 
original dev data to be able to finally close of this iteration `f' of the 
instability loop 

Note: pr_instability1 represents pr_original (as this is the original dev data),
whereas pr_instability2 - pr_instability201 are the 200 bootstrap estimates of the 
pr based on a bootstrap sample of the dev data (so will need to amend instability 
codes below accordingly)
*/

* 4.41 Load the current batch of development data again 

use "MasterDevData", replace

/* 4.42 Predict probabilities (manually now) & linear predictor from the bootstrap 
shrinkage model in the original dev dataset and save, alongside original model predictions */
* bootstrap model 1 predictions are pr1, model 2 is pr2, and so forth

* Now can use this code to apply the above HS shrinkage model to any dataset 
* and work out the predicted logp and p values ...

gen logp_instability`f' = 0

replace logp_instability`f' = $new_cons ///
+ bicarbonate_mean*$bicarbonate_mean + creatinine_mean*$creatinine_mean ///
+ hemoglobin_mean*$hemoglobin_mean + bun_mean*$bun_mean + ///
potassium_mean*$potassium_mean + sysbp_mean*$sysbp_mean 

gen pr_instability`f' = exp(logp_instability`f')/(1+exp(logp_instability`f'))

/* 4.43 Saving the dataset at each iteration of the instability loop (it's fine 
that it will get overwritten across iterations, as the final iteration before 
the loop closes for a particular sample size, will contain the 201th iteration 
of the instability loop, which is the final one that contains everything we need 
to produce summaries for a particular SS).

Note: we need to save, because further down in Section 5, for the instability 
metric summaries, we need to be able to open up (separately for each instability 
metric) the dataset that contains the pr_original and 200 bootstrap estimates of 
pr (from instability loop)

*/

save "MasterDevData", replace

}

/* 4.5 Renaming pr and logp variables, for ease, to make codes below on instability 
and pr risk regions etc run easily with old naming convention */

use "MasterDevData", replace

rename pr_instability1 pr_original

rename logp_instability1 logpr_original


forvalues h=2/201 {
	
local g = `h'-1

rename 	pr_instability`h' pr`g'

rename 	logp_instability`h' logp`g'

}

/* 4.6 Re-saving the final master dataset for the current sample size 
iteration */

save "MasterDevData", replace


/* 4.7 Creating an average EVPI across the f=201 iterations of the instability 
loop too, so that have a unique EVPI value to store for summarisation for a 
particular value of the SS loop */

/* 4.71 First, terminating the postfile dataset process (to create a dataset) 
for the 201 bootstrap EVPI values */

postclose `BSShrinkageFactor_SS`j'' 

/* 4.72 Then using the dataset created of 201 bootstrap EVPI values, and 
storing the mean EVPI as a local variable, to store for a particular SS in 
the final dataset */

use "BSShrinkageFactor_SS`j'.dta", replace

su EVPI_temp
local avg_EVPI = r(mean)

/* ************************************************************************** */
/* 5. FINALLY, USING BOOTSTRAP RESULTS TO COMPUTE INDIVIDUAL-LEVEL STABILITY  */
/* OF PREDICTIONS AND CLASSIFICATIONS (ACROSS THE B=200 BOOTSTRAP SAMPLES)    */
/*                                                                            */
/* ************************************************************************** */

/* 5.1 95% UNCERTAINTY INTERVAL (UI) WIDTH AND DELTA STATISTIC */

/* 5.11 For each individual, first record the 2.5th and 97.5th percentiles of 
predictions across the B=200 bootstrap models. */

use "MasterDevData", replace
order *, sequential
egen lower = rowpctile(pr1-pr200), p(2.5)
egen upper = rowpctile(pr1-pr200), p(97.5)

/* 5.12 Now calculating the following individual-level prediction instability 
measures:

1) Mean 95% Uncertainty Interval (UI) Width, 
2) Delta = max of (original point risk estimate (from original model developed in original dataset) - lower 95% UI value) 
and (upper 95% UI value - original point risk estimate (from original model developed in original dataset)) 

First calculate these statistics at the individual level */

gen UIwidth = upper - lower

gen pointriskminuslower = pr_original - lower
gen upperminuspointrisk = upper - pr_original

gen delta = max(pointriskminuslower, upperminuspointrisk)

/* 5.13 Now summarise mean UI width and mean delta (this is what I'll use to 
contribute towards instability-based stopping rule criteria), and median alongside
2.5th and 97.5th percentiles (as I'll actually use these to show on LCs?) */

qui su UIwidth, d 
local UIwidth_mean = r(mean)

_pctile UIwidth, p(2.5)
local UIwidth_pctile2dot5 = r(r1)

_pctile UIwidth, p(97.5)
local UIwidth_pctile97dot5 = r(r1)


qui su delta, d 
local Delta_mean = r(mean)

_pctile delta, p(2.5)
local Delta_pctile2dot5 = r(r1)

_pctile delta, p(97.5)
local Delta_pctile97dot5 = r(r1)


/* ************************************************************************** */

/* 5.2 MEAN ABSOLUTE PREDICTION ERROR (MAPE) */

/* 5.21 First working out absolute prediction 'error' based on p-scale */

use "MasterDevData", replace
qui forvalues k=1/200 {
	gen error_bs`k' = pr`k' - pr_original
	gen abs_error_bs`k' = abs(error_bs`k')
}
/* 5.22 Now work out the mean of this absolute prediction error (i.e. the MAPE).
Note: need to order the variables sequentially by name, to avoid rowmean issues */

order *, sequential
qui egen MAPE = rowmean(abs_error_bs1-abs_error_bs200)

/* 5.23 Estimating the mean of the MAPE, and 2.5th and 97.5th percentiles */

qui su MAPE, d 
local MAPE_mean = r(mean)

_pctile MAPE, p(2.5)
local MAPE_pctile2dot5 = r(r1)

_pctile MAPE, p(97.5)
local MAPE_pctile97dot5 = r(r1)

/* ************************************************************************** */

/* 5.3 CLASSIFICATION INSTABILITY */

/* 5.31 Examining the change in classifications based on a chosen risk threshold */

use "MasterDevData", replace
order *, sequential
qui forvalues k=1/200 {
	gen pr`k'_change = 0
	* assuming a risk threshold of >= 0.1 is relevant for defining 'high' risk
	replace pr`k'_change = 1 if pr`k' >= 0.100 & pr_orig < 0.100
	replace pr`k'_change = 1 if pr`k' < 0.100 & pr_orig >= 0.100
} 
qui egen pr_change = rowmean(pr1_change - pr200_change)

/* 5.32 Estimating the mean probability of misclassification (based on a 10% risk threshold),
alongside the 2.5th and 97.5th percentiles */

qui su pr_change, d 
local probmisclass_mean = r(mean)

_pctile pr_change, p(2.5)
local probmisclass_pctile2dot5 = r(r1)

_pctile pr_change, p(97.5)
local probmisclass_pctile97dot5 = r(r1)

/* ************************************************************************** */

/* 6.0 Finally, 'posting' all the overall model and individual-level performance 
measure values (mean, 2.5th percentile, 97.5th percentile) obtained from 
bootstrapping for this particular sample size iteration into the final output dataset */

post `FinalDataset' (`n') (`UIwidth_mean') (`UIwidth_pctile2dot5') (`UIwidth_pctile97dot5') (`Delta_mean') (`Delta_pctile2dot5') (`Delta_pctile97dot5') (`MAPE_mean') (`MAPE_pctile2dot5') (`MAPE_pctile97dot5') (`probmisclass_mean') (`probmisclass_pctile2dot5') (`probmisclass_pctile97dot5') (`avg_EVPI')

}

/* Once the sample size loop is complete (i.e. we have iteratively increased the 
sample size from 100 to 3000 in increments of 100 and stored the relevant 
performance measures each time, we can save the final dataset */

postclose `FinalDataset'


/* ************************************************************************** */
/* ************************************************************************** */
/* ************************************************************************** */
/* ONCE THE ABOVE CODE HAS RUN AND DATA HAS BEEN GENERATED, RUN THE CODES     */
/* BELOW TO GENERATE LEARNING CURVES AND CALCULATE WHEN THE SPECIFIED STOPPING*/
/* RULE CRITERIA IS FIRST MET                                                 */
/* ************************************************************************** */ 
/* ************************************************************************** */
/* ************************************************************************** */




/* ************************************************************************** */
/* 7.0: GENERATING LEARNING CURVES FOR EACH PERFORMANCE MEASURE OF INTEREST   */
/* ************************************************************************** */

/* 7.01 First open the final learning curve dataset */

use "FinalDataset", clear

/* 7.1 Learning Curve for: Mean UI width */

gen StoppingRuleCriteria = 0.1	

twoway line UIwidth_mean UIwidth_pctile2dot5 UIwidth_pctile97dot5 StoppingRuleCriteria n, ///
ytitle(Uncertainty Interval Width) xtitle(n) legend(order(1 "Mean" 2 "2.5th Percentile" 3 "97.5th Percentile" 4 "Stopping Rule Criteria")) ///
lcolor (midblue gray gray black) lpat(solid solid solid dash) lwidth(*4 *1 *1 *1) legend(size(*0.75)) xlabel(0(500)3000) ylabel(0.0(.1)0.5)

/* 7.2 Learning Curve for: Mean Delta */

drop StoppingRuleCriteria

gen StoppingRuleCriteria = 0.05	

twoway line Delta_mean Delta_pctile2dot5 Delta_pctile97dot5 StoppingRuleCriteria n, ///
ytitle(Delta) xtitle(n) legend(order(1 "Mean" 2 "2.5th Percentile" 3 "97.5th Percentile" 4 "Stopping Rule Criteria")) ///
lcolor (midblue gray gray black) lpat(solid solid solid dash) lwidth(*4 *1 *1 *1) legend(size(*0.75)) xlabel(0(500)3000) ylabel(0.0(.05)0.35)

/* 7.3 Learning Curve for: Mean MAPE */

drop StoppingRuleCriteria

gen StoppingRuleCriteria = 0.02	

twoway line MAPE_mean MAPE_pctile2dot5 MAPE_pctile97dot5 StoppingRuleCriteria n, ///
ytitle(MAPE (on p-scale)) xtitle(n) legend(order(1 "Mean" 2 "2.5th Percentile" 3 "97.5th Percentile" 4 "Stopping Rule Criteria")) ///
lcolor (midblue gray gray black) lpat(solid solid solid dash) lwidth(*4 *1 *1 *1) legend(size(*0.75)) xlabel(0(500)3000) ylabel(0.0(.02)0.1)				

/* 7.4 Learning Curve for: Mean Probability of Misclassification (10% threshold) */

drop StoppingRuleCriteria

gen StoppingRuleCriteria = 0.1	

twoway line probmisclass_mean probmisclass_pctile2dot5 probmisclass_pctile97dot5 StoppingRuleCriteria n, ///
ytitle(Probability of Misclassification (10% threshold)) xtitle(n) legend(order(1 "Mean" 2 "2.5th Percentile" 3 "97.5th Percentile" 4 "Stopping Rule Criteria")) ///
lcolor (midblue gray gray black) lpat(solid solid solid dash) lwidth(*4 *1 *1 *1) legend(size(*0.75)) xlabel(0(500)3000) ylabel(0.0(0.1)0.6)

/* 7.5 Learning Curve for: Mean EVPI */

drop StoppingRuleCriteria

gen StoppingRuleCriteria = 0.001

twoway line avg_EVPI StoppingRuleCriteria n, ///
ytitle(Mean of Average EVPI) xtitle(n) legend(order(1 "Mean" 2 "Stopping Rule Criteria")) ///
lcolor (midblue black) lpat(solid dash) lwidth(*4 *1) legend(size(*0.75)) xlabel(0(500)3000) ylabel(0.0(.001)0.006)


/* ************************************************************************** */
/* 8.0: NOW CALCULATING THE SAMPLE SIZE AT WHICH THE STOPPING RULE CRITERIA   */
/* ARE FIRST MET (FOR CHOICE OF PERFORMANCE MEASURE AND ASSOCIATED            */
/* STOPPING RULE), AND SUSTAINED OVER TWO CONSECUTIVE SAMPLE SIZE INCREMENTS  */
/* (AS WELL AS 3 AND 5 CONSECUTIVE SAMPLE SIZE INCREMENTS, RESPECTIVELY, AS   */
/* A SENSITIVITY ANALYSIS)                                                    */ 
/* ************************************************************************** */

/* ************************************************************************** */
/* 8.1 Minimum Required Sample Size Based on Stopping Rule Criteria of:       */
/* Mean UI width <=0.10                                                       */
/* ************************************************************************** */

/* 8.11 First open the final learning curve dataset, and just keep the variables 
for 'n' (sample size increment) and the performance measure of interest */

use "FinalDataset", clear

keep n UIwidth_mean

/* 8.12 Re-sort data so that n goes from smallest to largest */

sort n

/* 8.13 Now assess the performance measure value for the next sample size increment 
(i.e. n = current SS + 100),
so can check if stopping rule criteria sustained over two consecutive SS increments
(or further increments as necessary, here I'm checking for 2, 3 and 5, 
consecutive increments, respectively), as this is what I require for each 
stopping rule to be considered as 'met' */

gen NextSSMeasureValue = UIwidth_mean[_n+1]
gen ThirdConsecutiveSSMeasureValue = UIwidth_mean[_n+2]
gen FourthConsecutiveSSMeasureValue = UIwidth_mean[_n+3]
gen FifthConsecutiveSSMeasureValue = UIwidth_mean[_n+4]

/* 8.14 Create a variable to indicate whether the particular stopping rule criteria
has been met over two consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetTwice = 0 
replace StoppingRuleMetTwice = 1 if ((UIwidth_mean <= 0.1) & (NextSSMeasureValue <= 0.1))

summ n if StoppingRuleMetTwice == 1

di as err "Minimum sequential SS based on stopping rule of mean UI width<=0.1 and two consecutive occurrences is: " r(min)

/* 8.15 Create a variable to indicate whether the particular stopping rule criteria
has been met over three consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetThreeTimes = 0 
replace StoppingRuleMetThreeTimes = 1 if ((UIwidth_mean <= 0.1) & (NextSSMeasureValue <= 0.1) & (ThirdConsecutiveSSMeasureValue <= 0.1))

summ n if StoppingRuleMetThreeTimes == 1

di as err "Minimum sequential SS based on stopping rule of mean UI width<=0.1 and three consecutive occurrences is: " r(min)

/* 8.16 Create a variable to indicate whether the particular stopping rule criteria
has been met over five consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetFiveTimes = 0 
replace StoppingRuleMetFiveTimes = 1 if ((UIwidth_mean <= 0.1) & (NextSSMeasureValue <= 0.1) & (ThirdConsecutiveSSMeasureValue <= 0.1) & (FourthConsecutiveSSMeasureValue <= 0.1) & (FifthConsecutiveSSMeasureValue <= 0.1))

summ n if StoppingRuleMetFiveTimes == 1

di as err "Minimum sequential SS based on stopping rule of mean UI width<=0.1 and five consecutive occurrences is: " r(min)

/* ************************************************************************** */
/* 8.2 Minimum Required Sample Size Based on Stopping Rule Criteria of:       */
/* mean Delta <=0.05                                                          */
/* ************************************************************************** */

/* 8.21 First open the final learning curve dataset, and just keep the variables 
for 'n' (sample size increment) and the performance measure of interest */

use "FinalDataset", clear

keep n Delta_mean

/* 8.22 Re-sort data so that n goes from smallest to largest */

sort n

/* 8.23 Now assess the performance measure value for the next sample size increment 
(i.e. n = current SS + 100),
so can check if stopping rule criteria sustained over two consecutive SS increments
(or further increments as necessary, here I'm checking for 2, 3 and 5, 
consecutive increments, respectively), as this is what I require for each 
stopping rule to be considered as 'met' */

gen NextSSMeasureValue = Delta_mean[_n+1]
gen ThirdConsecutiveSSMeasureValue = Delta_mean[_n+2]
gen FourthConsecutiveSSMeasureValue = Delta_mean[_n+3]
gen FifthConsecutiveSSMeasureValue = Delta_mean[_n+4]

/* 8.24 Create a variable to indicate whether the particular stopping rule criteria
has been met over two consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetTwice = 0 
replace StoppingRuleMetTwice = 1 if ((Delta_mean <= 0.05) & (NextSSMeasureValue <= 0.05))

summ n if StoppingRuleMetTwice == 1

di as err "Minimum sequential SS based on stopping rule of mean Delta<=0.05 and two consecutive occurrences is: " r(min)

/* 8.25 Create a variable to indicate whether the particular stopping rule criteria
has been met over three consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetThreeTimes = 0 
replace StoppingRuleMetThreeTimes = 1 if ((Delta_mean <= 0.05) & (NextSSMeasureValue <= 0.05) & (ThirdConsecutiveSSMeasureValue <= 0.05))

summ n if StoppingRuleMetThreeTimes == 1

di as err "Minimum sequential SS based on stopping rule of mean Delta<=0.05 and three consecutive occurrences is: " r(min)

/* 8.26 Create a variable to indicate whether the particular stopping rule criteria
has been met over five consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetFiveTimes = 0 
replace StoppingRuleMetFiveTimes = 1 if ((Delta_mean <= 0.05) & (NextSSMeasureValue <= 0.05) & (ThirdConsecutiveSSMeasureValue <= 0.05) & (FourthConsecutiveSSMeasureValue <= 0.05) & (FifthConsecutiveSSMeasureValue <= 0.05))

summ n if StoppingRuleMetFiveTimes == 1

di as err "Minimum sequential SS based on stopping rule of mean Delta<=0.05 and five consecutive occurrences is: " r(min)

/* ************************************************************************** */
/* 8.3 Minimum Required Sample Size Based on Stopping Rule Criteria of:       */
/* mean MAPE <=0.02                                                           */
/* ************************************************************************** */

/* 8.31 First open the final learning curve dataset, and just keep the variables 
for 'n' (sample size increment) and the performance measure of interest */

use "FinalDataset", clear

keep n MAPE_mean

/* 8.32 Re-sort data so that n goes from smallest to largest */

sort n

/* 8.33 Now assess the performance measure value for the next sample size increment 
(i.e. n = current SS + 100),
so can check if stopping rule criteria sustained over two consecutive SS increments
(or further increments as necessary, here I'm checking for 2, 3 and 5, 
consecutive increments, respectively), as this is what I require for each 
stopping rule to be considered as 'met' */

gen NextSSMeasureValue = MAPE_mean[_n+1]
gen ThirdConsecutiveSSMeasureValue = MAPE_mean[_n+2]
gen FourthConsecutiveSSMeasureValue = MAPE_mean[_n+3]
gen FifthConsecutiveSSMeasureValue = MAPE_mean[_n+4]

/* 8.34 Create a variable to indicate whether the particular stopping rule criteria
has been met over two consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetTwice = 0 
replace StoppingRuleMetTwice = 1 if ((MAPE_mean <= 0.02) & (NextSSMeasureValue <= 0.02))

summ n if StoppingRuleMetTwice == 1

di as err "Minimum sequential SS based on stopping rule of mean MAPE<=0.02 and two consecutive occurrences is: " r(min)

/* 8.35 Create a variable to indicate whether the particular stopping rule criteria
has been met over three consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetThreeTimes = 0 
replace StoppingRuleMetThreeTimes = 1 if ((MAPE_mean <= 0.02) & (NextSSMeasureValue <= 0.02) & (ThirdConsecutiveSSMeasureValue <= 0.02))

summ n if StoppingRuleMetThreeTimes == 1

di as err "Minimum sequential SS based on stopping rule of mean MAPE<=0.02 and three consecutive occurrences is: " r(min)

/* 8.36 Create a variable to indicate whether the particular stopping rule criteria
has been met over five consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetFiveTimes = 0 
replace StoppingRuleMetFiveTimes = 1 if ((MAPE_mean <= 0.02) & (NextSSMeasureValue <= 0.02) & (ThirdConsecutiveSSMeasureValue <= 0.02) & (FourthConsecutiveSSMeasureValue <= 0.02) & (FifthConsecutiveSSMeasureValue <= 0.02))

summ n if StoppingRuleMetFiveTimes == 1

di as err "Minimum sequential SS based on stopping rule of mean MAPE<=0.02 and five consecutive occurrences is: " r(min)


/* ************************************************************************** */
/* 8.4 Minimum Required Sample Size Based on Stopping Rule Criteria of:       */
/* Mean Probability of Misclassification <=0.1                                */
/* ************************************************************************** */

/* 8.41 First open the final learning curve dataset, and just keep the variables 
for 'n' (sample size increment) and the performance measure of interest */

use "FinalDataset", clear

keep n probmisclass_mean

/* 8.42 Re-sort data so that n goes from smallest to largest */

sort n

/* 8.43 Now assess the performance measure value for the next sample size increment 
(i.e. n = current SS + 100),
so can check if stopping rule criteria sustained over two consecutive SS increments
(or further increments as necessary, here I'm checking for 2, 3 and 5, 
consecutive increments, respectively), as this is what I require for each 
stopping rule to be considered as 'met' */

gen NextSSMeasureValue = probmisclass_mean[_n+1]
gen ThirdConsecutiveSSMeasureValue = probmisclass_mean[_n+2]
gen FourthConsecutiveSSMeasureValue = probmisclass_mean[_n+3]
gen FifthConsecutiveSSMeasureValue = probmisclass_mean[_n+4]

/* 8.44 Create a variable to indicate whether the particular stopping rule criteria
has been met over two consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetTwice = 0 
replace StoppingRuleMetTwice = 1 if ((probmisclass_mean <= 0.1) & (NextSSMeasureValue <= 0.1))

summ n if StoppingRuleMetTwice == 1

di as err "Minimum sequential SS based on stopping rule of Mean Probability of Misclassification<=0.1 and two consecutive occurrences is: " r(min)

/* 8.45 Create a variable to indicate whether the particular stopping rule criteria
has been met over three consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetThreeTimes = 0 
replace StoppingRuleMetThreeTimes = 1 if ((probmisclass_mean <= 0.1) & (NextSSMeasureValue <= 0.1) & (ThirdConsecutiveSSMeasureValue <= 0.1))

summ n if StoppingRuleMetThreeTimes == 1

di as err "Minimum sequential SS based on stopping rule of Mean Probability of Misclassification<=0.1 and three consecutive occurrences is: " r(min)

/* 8.46 Create a variable to indicate whether the particular stopping rule criteria
has been met over five consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetFiveTimes = 0 
replace StoppingRuleMetFiveTimes = 1 if ((probmisclass_mean <= 0.1) & (NextSSMeasureValue <= 0.1) & (ThirdConsecutiveSSMeasureValue <= 0.1) & (FourthConsecutiveSSMeasureValue <= 0.1) & (FifthConsecutiveSSMeasureValue <= 0.1))

summ n if StoppingRuleMetFiveTimes == 1

di as err "Minimum sequential SS based on stopping rule of Mean Probability of Misclassification<=0.1 and five consecutive occurrences is: " r(min)



/* ************************************************************************** */
/* 8.5 Minimum Required Sample Size Based on Stopping Rule Criteria of:       */
/* Mean EVPI <=0.001                                                          */
/* ************************************************************************** */

/* 8.51 First open the final learning curve dataset, and just keep the variables 
for 'n' (sample size increment) and the performance measure of interest */

use "FinalDataset", clear

keep n avg_EVPI

/* 8.52 Re-sort data so that n goes from smallest to largest */

sort n

/* 8.53 Now assess the performance measure value for the next sample size increment 
(i.e. n = current SS + 100),
so can check if stopping rule criteria sustained over two consecutive SS increments
(or further increments as necessary, here I'm checking for 2, 3 and 5, 
consecutive increments, respectively), as this is what I require for each 
stopping rule to be considered as 'met' */

gen NextSSMeasureValue = avg_EVPI[_n+1]
gen ThirdConsecutiveSSMeasureValue = avg_EVPI[_n+2]
gen FourthConsecutiveSSMeasureValue = avg_EVPI[_n+3]
gen FifthConsecutiveSSMeasureValue = avg_EVPI[_n+4]

/* 8.54 Create a variable to indicate whether the particular stopping rule criteria
has been met over two consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetTwice = 0 
replace StoppingRuleMetTwice = 1 if ((avg_EVPI <= 0.001) & (NextSSMeasureValue <= 0.001))

summ n if StoppingRuleMetTwice == 1

di as err "Minimum sequential SS based on stopping rule of Mean EVPI<=0.001 and two consecutive occurrences is: " r(min)

/* 8.55 Create a variable to indicate whether the particular stopping rule criteria
has been met over three consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetThreeTimes = 0 
replace StoppingRuleMetThreeTimes = 1 if ((avg_EVPI <= 0.001) & (NextSSMeasureValue <= 0.001) & (ThirdConsecutiveSSMeasureValue <= 0.001))

summ n if StoppingRuleMetThreeTimes == 1

di as err "Minimum sequential SS based on stopping rule of Mean EVPI<=0.001 and three consecutive occurrences is: " r(min)

/* 8.56 Create a variable to indicate whether the particular stopping rule criteria
has been met over five consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetFiveTimes = 0 
replace StoppingRuleMetFiveTimes = 1 if ((avg_EVPI <= 0.001) & (NextSSMeasureValue <= 0.001) & (ThirdConsecutiveSSMeasureValue <= 0.001) & (FourthConsecutiveSSMeasureValue <= 0.001) & (FifthConsecutiveSSMeasureValue <= 0.001))

summ n if StoppingRuleMetFiveTimes == 1

di as err "Minimum sequential SS based on stopping rule of Mean EVPI<=0.001 and five consecutive occurrences is: " r(min)
























