
/* ************************************************************************** */
/*        0. IMPORTING IN DATA AND CREATING AN EMPTY POSTFILE DATASET         */
/* ************************************************************************** */


clear all
set more off
set maxvar 5000

/* 0.0 Set working directory */

cd "C:\Users\User XYZ\Stata Results Tables_Logistic regression with LASSO penalty\"

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
EVPI ///
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
iteration of the loop). */

forvalues j = 3000(-100)100 {

/* 1.00 Reloading the dataset called 'mimic_reducedvars' that contains 
the current pool of patients that could potentially be randomly selected from 
the MIMIC data for this particular iteration of the sample size loop */
	
use "mimic_reducedvars", replace	

/* 1.01 Creating a local value for n, so we can 'post' it to the final output 
dataset later on */

local n = `j'

* 1.02 Adding a visual counter of the sequential SS iteration (n)

di as err "Currently on n = " `j'

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

/* ************************************************************************** */
/* 1.2 FITTING THE CHOSEN MODELLING STRATEGY ONTO CURRENT DEVELOPMENT DATASET,*/
/* AND STORING OVERALL MODEL PERFORMANCE MEASURES, AS WELL AS INDIVIDUAL RISK */
/* PREDICTIONS */
/* ************************************************************************** */

/* 1.2 Apply the initial model which is the logistic regression with LASSO 
penalty model to the development dataset: 
We are fitting a logistic regression model to determine the probability of 
developing AKI (outcome) based on 6 (continuous) predictors 
(linear predictor-outcome relationships are assumed, and all 6 variables 
are forced into the model). */

* Changed max iterations from 300 to 500 to improve model convergence rate

qui capture lasso logit AKI bicarbonate_mean creatinine_mean hemoglobin_mean ///
bun_mean potassium_mean sysbp_mean, selection(cv)

/* 1.21 Calculate the predicted probabilities (let's name them 'pr_original') 
and the value of the linear predictor (i.e. 'pr_original' but on logit-p scale 
rather than p-scale, let's name: 'logpr_original') for each individual using 
the predict command. */

predict pr_original, pr

gen logpr_original = logit(pr_original)

/* 1.2 Save predictions and original data in a new dataset */

save "LASSO_originaldevdatawithmodelresults", replace

* ************************************************************************** */
/* 2. START OF SECOND (FINAL AND INNERMOST) LOOP: 'B' BOOTSTRAP REPS TO:      */
/*                                                                            */
/* i) GENERATE B=200 RISK ESTIMATES PER INDIVIDUAL (FOR INSTABILITY MEASURES) */
/* ************************************************************************** */

/*
	Note: we will store 3x measures relating to EVPI (these only relate to 
	the bootstrap model applied to the original development data); these will 
	be stored in a matrix (B=200 bootstrap reps, and the 3x EVPI intermediary 
	variables) 
	
	Additionally, we will store the individual risk estimates, separately in a 
	dataset (not matrix), of the n=`i' observations for each of the B=200 
	bootstrap replications. 
	
	Iteratively storing these individual risk estimates into the same dataset 
	created in step 1.4 of code, that contains the original model risk estimates: 
	'HS_originaldevdatawithmodelresults.dta'.
	
	So we at the end of this bootstrap, we are storing 
	(i) a matrix of EVPI intermediary measures, and
	(ii) a dataset of individual risk estimates across the B=200 bootstrap reps.

*/

matrix results = J(200,3,.)

qui forvalues k = 1/200 {

/* 2.01 Loading the current batch of development data and just keeping the 
outcome and 6 predictors (i.e. removing the predicted probabilities from the 
original model) */

use "LASSO_originaldevdatawithmodelresults", clear
keep AKI bicarbonate_mean creatinine_mean hemoglobin_mean bun_mean potassium_mean sysbp_mean

* 2.02 Taking a bootstrap sample (random sample with replacement) of same size
bsample

* 2.03 Iteration indicator to help the user know the progress
nois _dots `k' 0

/* ************************************************************************** */

/* 3.1 Apply the initial model which is the logistic regression with LASSO 
penalty model to the development dataset: 
We are fitting a logistic regression model to determine the probability of 
developing AKI (outcome) based on 6 (continuous) predictors 
(linear predictor-outcome relationships are assumed, and all 6 variables 
are forced into the model). */

* Changed max iterations from 300 to 500 to improve model convergence rate

capture lasso logit AKI bicarbonate_mean creatinine_mean hemoglobin_mean ///
bun_mean potassium_mean sysbp_mean, selection(cv)

/* ************************************************************************** */

/* 3.2 Now testing performance of the bootstrap model in original sample */

/* 3.21 Load the current batch of development data again */
use "LASSO_originaldevdatawithmodelresults", clear

/* 3.22 Predict probabilities & linear predictor from the bootstrap model 
 applied to the original dataset and save, alongside original model predictions.
 Bootstrap model 1 predictions are saved as 'pr1', model 2 as 'pr2', and so forth
 up to 'pr200' */

predict pr`k', pr
gen logp`k' = logit(pr`k')
save "LASSO_originaldevdatawithmodelresults", replace 

/* 3.23 Now performing intermediary calculations to enable bootstrap-based computation 
of the Expected Value of Perfect Information (EVPI) */

/* 3.231 First estimating 'NB_All' for this current kth bootstrap rep 
(NB_All is the net benefit of administering a treatment (for example, closer 
monitoring of renal function) to all patients regardless of their underlying risk) 

For full mechanism and corresponding equations see: 
Sadatsafavi M, Yoon Lee T, Gustafson P. Uncertainty and the Value of Information in Risk Prediction Modeling. Medical Decision Making. 2022;42(5):661-71.
*/

gen NB_All_`k' = ((pr`k'-(1-pr`k')*0.1)/(0.9))

qui su NB_All_`k', d

matrix results[`k',1] = r(mean)

/* 3.232 Then estimating 'NB_Model' (NB_Model is the net benefit of using the 
proposed model to decide which patients to treat)
 */

gen prorig_greater_z_indicator = 0
replace prorig_greater_z_indicator = 1 if pr_original > 0.1 


gen NB_Model_`k' = (prorig_greater_z_indicator*((pr`k'-(1-pr`k')*0.1)/(0.9)))

qui su NB_Model_`k', d

matrix results[`k',2] = r(mean)

/* 3.233 Finally estimating 'NB_Max' (NB_Max is the net benefit of using the 
`correct' model to decide which patients to treat) */

gen prbootstrap_greater_z_indicator = 0
replace prbootstrap_greater_z_indicator = 1 if pr`k' > 0.1 

gen NB_Max_`k' = (prbootstrap_greater_z_indicator*((pr`k'-(1-pr`k')*0.1)/(0.9)))

qui su NB_Max_`k', d

matrix results[`k',3] = r(mean)

}



/* ************************************************************************** */
/* 4. USING BOOTSTRAP RESULTS TO COMPUTE EVPI                                 */
/* (AN OVERALL POPULATION-LEVEL STABILITY MEASURE)                            */
/* ************************************************************************** */

/* 4.0 Rename columns of the matrix used to store bootstrap overall model performance results */

mat colnames results = NB_All NB_Model NB_Max

/* 4.01 Clear and load these matrix results into a Stata dataset */

clear
svmat results, n(col)

/* 4.1 EVPI MEASURES */

/* 4.11 First calculating the ENBall, ENBModel, and ENBMax by taking the 
mean of the 200 bootstrap values for these measures.

As in Section 3.25, for full mechanism and corresponding equations see: 
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

local EVPI = `ENB_Max'-max(0,`ENB_Model',`ENB_All')

/* ************************************************************************** */
/* 5. FINALLY, USING BOOTSTRAP RESULTS TO COMPUTE INDIVIDUAL-LEVEL STABILITY  */
/* OF PREDICTIONS AND CLASSIFICATIONS (ACROSS THE B=200 BOOTSTRAP SAMPLES)    */
/*                                                                            */
/* ************************************************************************** */

/* 5.1 95% UNCERTAINTY INTERVAL (UI) WIDTH AND DELTA STATISTIC */

/* 5.11 For each individual, first record the 2.5th and 97.5th percentiles of 
predictions across the B=200 bootstrap models. */

use "LASSO_originaldevdatawithmodelresults", replace
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

use "LASSO_originaldevdatawithmodelresults", replace
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

/* 5.31 First examining which of the B=200 bootstrap model predictions in the original 
dataset result in a classification that is different to the classification from 
the original prediction (of the original model applied to the original data)/

I.e. examining the change in classifications based on a chosen risk threshold */

use "LASSO_originaldevdatawithmodelresults", replace
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

post `FinalDataset' (`n') (`UIwidth_mean') (`UIwidth_pctile2dot5') (`UIwidth_pctile97dot5') (`Delta_mean') (`Delta_pctile2dot5') (`Delta_pctile97dot5') (`MAPE_mean') (`MAPE_pctile2dot5') (`MAPE_pctile97dot5') (`probmisclass_mean') (`probmisclass_pctile2dot5') (`probmisclass_pctile97dot5') (`EVPI')

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
lcolor (midblue gray gray black) lpat(solid solid solid dash) lwidth(*4 *1 *1 *1) legend(size(*0.75)) xlabel(0(500)3000) ylabel(0.0(.1)0.6)

/* 7.2 Learning Curve for: Mean Delta */

drop StoppingRuleCriteria

gen StoppingRuleCriteria = 0.05	

twoway line Delta_mean Delta_pctile2dot5 Delta_pctile97dot5 StoppingRuleCriteria n, ///
ytitle(Delta) xtitle(n) legend(order(1 "Mean" 2 "2.5th Percentile" 3 "97.5th Percentile" 4 "Stopping Rule Criteria")) ///
lcolor (midblue gray gray black) lpat(solid solid solid dash) lwidth(*4 *1 *1 *1) legend(size(*0.75)) xlabel(0(500)3000) ylabel(0.0(.05)0.5)

/* 7.3 Learning Curve for: Mean MAPE */

drop StoppingRuleCriteria

gen StoppingRuleCriteria = 0.02	

twoway line MAPE_mean MAPE_pctile2dot5 MAPE_pctile97dot5 StoppingRuleCriteria n, ///
ytitle(MAPE (on p-scale)) xtitle(n) legend(order(1 "Mean" 2 "2.5th Percentile" 3 "97.5th Percentile" 4 "Stopping Rule Criteria")) ///
lcolor (midblue gray gray black) lpat(solid solid solid dash) lwidth(*4 *1 *1 *1) legend(size(*0.75)) xlabel(0(500)3000) ylabel(0.0(.02)0.12)				

/* 7.4 Learning Curve for: Mean Probability of Misclassification (10% threshold) */

drop StoppingRuleCriteria

gen StoppingRuleCriteria = 0.1	

twoway line probmisclass_mean probmisclass_pctile2dot5 probmisclass_pctile97dot5 StoppingRuleCriteria n, ///
ytitle(Probability of Misclassification (10% threshold)) xtitle(n) legend(order(1 "Mean" 2 "2.5th Percentile" 3 "97.5th Percentile" 4 "Stopping Rule Criteria")) ///
lcolor (midblue gray gray black) lpat(solid solid solid dash) lwidth(*4 *1 *1 *1) legend(size(*0.75)) xlabel(0(500)3000) ylabel(0.0(0.1)0.7)

/* 7.5 Learning Curve for: Mean EVPI */

drop StoppingRuleCriteria

gen StoppingRuleCriteria = 0.001

twoway line EVPI StoppingRuleCriteria n, ///
ytitle(Mean EVPI) xtitle(n) legend(order(1 "Mean" 2 "Stopping Rule Criteria")) ///
lcolor (midblue black) lpat(solid dash) lwidth(*4 *1) legend(size(*0.75)) xlabel(0(500)3000) ylabel(0.0(.001)0.008)


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

keep n EVPI

/* 8.52 Re-sort data so that n goes from smallest to largest */

sort n

/* 8.53 Now assess the performance measure value for the next sample size increment 
(i.e. n = current SS + 100),
so can check if stopping rule criteria sustained over two consecutive SS increments
(or further increments as necessary, here I'm checking for 2, 3 and 5, 
consecutive increments, respectively), as this is what I require for each 
stopping rule to be considered as 'met' */

gen NextSSMeasureValue = EVPI[_n+1]
gen ThirdConsecutiveSSMeasureValue = EVPI[_n+2]
gen FourthConsecutiveSSMeasureValue = EVPI[_n+3]
gen FifthConsecutiveSSMeasureValue = EVPI[_n+4]

/* 8.54 Create a variable to indicate whether the particular stopping rule criteria
has been met over two consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetTwice = 0 
replace StoppingRuleMetTwice = 1 if ((EVPI <= 0.001) & (NextSSMeasureValue <= 0.001))

summ n if StoppingRuleMetTwice == 1

di as err "Minimum sequential SS based on stopping rule of Mean EVPI<=0.001 and two consecutive occurrences is: " r(min)

/* 8.55 Create a variable to indicate whether the particular stopping rule criteria
has been met over three consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetThreeTimes = 0 
replace StoppingRuleMetThreeTimes = 1 if ((EVPI <= 0.001) & (NextSSMeasureValue <= 0.001) & (ThirdConsecutiveSSMeasureValue <= 0.001))

summ n if StoppingRuleMetThreeTimes == 1

di as err "Minimum sequential SS based on stopping rule of Mean EVPI<=0.001 and three consecutive occurrences is: " r(min)

/* 8.56 Create a variable to indicate whether the particular stopping rule criteria
has been met over five consecutive sample size increments or not, and if so, 
when is the minimum sample size that this occurs */

gen StoppingRuleMetFiveTimes = 0 
replace StoppingRuleMetFiveTimes = 1 if ((EVPI <= 0.001) & (NextSSMeasureValue <= 0.001) & (ThirdConsecutiveSSMeasureValue <= 0.001) & (FourthConsecutiveSSMeasureValue <= 0.001) & (FifthConsecutiveSSMeasureValue <= 0.001))

summ n if StoppingRuleMetFiveTimes == 1

di as err "Minimum sequential SS based on stopping rule of Mean EVPI<=0.001 and five consecutive occurrences is: " r(min)



















