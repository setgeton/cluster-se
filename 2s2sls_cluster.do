vers 14.0

clear
set more off

set seed 13452

save results.dta, emptyok replace

drop _all

display "Start: $S_TIME  $S_DATE"

loc numiter = 10000
local nobs = 200
loc exp = 10 // create x obs per observation--> introduce cluster dependency
loc het = 1 // insert heteroskedasticity (1/0)
loc clu = 1 // insert cluster variance

quietly {

forvalues mciter=1/`numiter' {
noisily : if mod(`mciter',10)==0 di "`mciter' replications"

****GENERATE SYNTHETIC DATA***
set mem 1g
set obs `nobs' //create 600 individuals
gen id=_n //Create ID numbers
 //set random number seed for reproducibility


* unobservables, indep over individ
loc rhoA=0.3
loc rhoB=-0.2
matrix covmat=(1 , `rhoA', `rhoB' \ `rhoA',1,`rhoA'*`rhoB' \ `rhoB', `rhoA'*`rhoB',1)
drawnorm alpha1 alpha2 alpha3, cov(covmat)




* generate explanatory vars, individual specific
gen w=0+1 * uniform() 
gen z1=0+1 * uniform()
gen z2=0+1 * uniform() 
gen z3=0+1 * uniform()  

* expand in dimension t
expand `exp'
bys id: gen t=_n


* explanatory, individual and time specific component added
scalar define sigma_expl_noise =1
replace w = w + sigma_expl_noise * invnormal(uniform())
replace z1 = z1 + sigma_expl_noise * invnormal(uniform())
replace z2 = z2 + sigma_expl_noise * invnormal(uniform())
replace z3 = z3 + sigma_expl_noise * invnormal(uniform())


scalar define beta0 = 0
scalar define beta1 = 0.9
scalar define beta2 = -0.6
scalar define betawx1 = 0.2
scalar define betawx2 = 0.4
scalar define betawy = 0.2

* unobservables, indep over observ and time 
loc rho1=`rhoA'
loc rho2=`rhoB'
matrix covmat=(1 , `rho1', `rho2' \ `rho1',1,`rho1'*`rho2' \ `rho2', `rho1'*`rho2',1)
drawnorm u1 u2 u3, cov(covmat)



loc hetg1=0+`het'*1
loc hetg2=0+`het'*-1
loc hetg3=0+`het'*0.5
loc hetgc1=0+`het'*1
loc hetgc2=0+`het'*0.5
loc hetgc3=0+`het'*1
loc clu1=`clu'*1
loc clu2=`clu'*1
loc clu3=`clu'*1
*****generate x1, x2 and y
gen x1= -0.6*z1 - 1.60*z2 + 1.60*z3 +w*betawx1 + sqrt(exp(`hetg3'*z3+`hetg2'*z2+`hetg1'*z1))*u1 + `clu1'*sqrt(exp(`hetgc3'*z3+`hetgc2'*z2+`hetgc1'*z1))*alpha1
gen x2= 0.5*z1 - 0.50*z2 - 1.80*z3 +w*betawx2 + sqrt(exp(`hetg3'*z3+`hetg2'*z2+`hetg1'*z1))*u2 + `clu2'*sqrt(exp(`hetgc3'*z3+`hetgc2'*z2+`hetgc1'*z1))*alpha2
gen y = x1*beta1 + x2*beta2 + w*betawy + beta0 + sqrt(exp(`hetg3'*z3+`hetg2'*z2+`hetg1'*z1))*u3 + `clu3'*sqrt(exp(`hetgc3'*z3+`hetgc2'*z2+`hetgc1'*z1))*alpha3



***SAVE SAMPLE1
preserve  
keep y w z1 z2 z3 id
keep if id<`nobs'*0.4
sleep 200
save sample1.dta, replace //2stage data set
restore


***SAVE SAMPLE2
preserve
keep x1 x2 w z1 z2 z3 id
keep if id>=`nobs'*0.6
sleep 200
save sample2.dta, replace //1stage data set
restore


set matsize 1200



********ESTIMATION*************

use sample2.dta, clear
gen const = 1
*robust gmm
qui gmm (x1 - {xb1: z1 z2 z3 w const}) ///
(x2 - {xb2:z1 z2 z3 w const}), ///
instruments(1 2: z1 z2 z3 w) ///
winit(unadjusted,independent) onestep ///
deriv(1/xb1 = -1) ///
deriv(2/xb2 = -1)
mat Vx2het = e(V) /*Robust variance estimate of pix2*/
matlist Vx2het

**NEW: cluster-robust gmm***********
qui gmm (x1 - {xb1: z1 z2 z3 w const}) ///
(x2 - {xb2:z1 z2 z3 w const}), ///
instruments(1 2: z1 z2 z3 w) ///
winit(unadjusted,independent) onestep ///
deriv(1/xb1 = -1) ///
deriv(2/xb2 = -1) ///
vce(cluster id)
mat Vx2clu = e(V) /*CLUSTER-Robust variance estimate of pix2*/
matlist Vx2clu
*non-robust seemingly unrelated regressions
qui sureg (x1 x2 = z1 z2 z3 w )
mat Vx2hom = e(V) /*Non-robust variance estimate of pix2*/
matlist Vx2hom


use sample1.dta, clear
/*Generating predicted X*/
qui predict x1h, equation(x1)
qui predict x2h, equation(x2)
scalar kx = 2 /*Number of predicted variables, here x1 and x2*/
scalar ke = 2 /*Number of exogenous variables, here w and constant*/
qui reg y z1 z2 z3 w
mat Vy1hom = e(V)*e(df_r)/e(N) /*Non-robust variance estimate of piy1,*/
di e(df_r)
di e(df_m) 
di e(N_clust)
di e(N)
matlist e(V)
/*without degrees of freedom correction*/
reg y z1 z2 z3 w, rob
mat Vy1het = e(V)*e(df_r)/e(N) /*Robust variance estimate of piy1,*/
di e(df_r)
di e(df_m) 
di e(N_clust)
di e(N)
matlist e(V)
/*without degrees of freedom correction*/
reg y z1 z2 z3 w, cluster(id)
di e(df_r)
di e(df_m) 
di e(N_clust)
di e(N)
mat Vy1clu = e(V)*(e(df_r)/(e(df_r)+1)) // linreg correction: e(V)*((e(df_r)-e(df_m))/e(df_r))*(e(N)-e(df_m))/e(N) /*Cluster-Robust variance estimate of piy1,*/
matlist e(V)

/*without degrees of freedom correction*/

/*TS2SLS estimator*/
qui reg y x1h x2h w
mat b2s = e(b)
mat b2sx = b2s[1,1..kx]'  /*Selecting beta for predicted X only*/
/*Constructing C hat*/
qui reg z1 x1h x2h w
mat ch=e(b)
mat ch = ch'
qui reg z2 x1h x2h w
mat ch = ch,e(b)'
qui reg z3 x1h x2h w
mat ch = ch,e(b)'
mat ch = ch,(J(kx,ke,0)\I(ke)) /*Adjusting ch for the exogenous variables*/

/*Calculating non-robust standard errors*/
mat var1hom = ch*Vy1hom*ch' + (b2sx' # ch)*Vx2hom*(b2sx # ch')
mat seb2shom = vecdiag(cholesky(diag(vecdiag(var1hom))))'
/*Calculating robust standard errors*/
mat var1het = ch*Vy1het*ch' + (b2sx' # ch)*Vx2het*(b2sx # ch')
mat seb2shet = vecdiag(cholesky(diag(vecdiag(var1het))))'
/*NEW: Calculating cluster-robust standard errors*/
mat var1clu = ch*Vy1clu*ch' + (b2sx' # ch)*Vx2clu*(b2sx # ch')
mat seb2sclu = vecdiag(cholesky(diag(vecdiag(var1clu))))'

/*Displaying the results*/
mat results = b2s',seb2shom,seb2shet,seb2sclu
mat colnames results = b_ts2sls "hom-se" "rob-se" "cluster-se"
mat rownames results = x1 x2 w _cons
matlist results

gen b1=b2s[1,1]
gen b2=b2s[1,2]
gen bw=b2s[1,3]
gen sehom1=seb2shom[1,1]
gen sehom2=seb2shom[2,1]
gen sehom3=seb2shom[3,1]
gen sehet1=seb2shet[1,1]
gen sehet2=seb2shet[2,1]
gen sehet3=seb2shet[3,1]
gen seclu1=seb2sclu[1,1]
gen seclu2=seb2sclu[2,1]
gen seclu3=seb2sclu[3,1]
gen iter=`mciter'
keep if _n==1
keep b1-seclu3 iter
sleep 200
   append using results.dta
sleep 200
save results.dta, replace
drop _all
}

 } //quietly
sleep 200

**** ORGANIZE RESULTS/ COMPUTE STD DEV AND MEANS
use results.dta, clear

gen waldhom1=abs(b1-beta1)/sehom1>1.96
gen waldhom2=abs(b2-beta2)/sehom2>1.96
gen waldhom3=abs(bw-betawy)/sehom3>1.96

gen waldhet1=abs(b1-beta1)/sehet1>1.96
gen waldhet2=abs(b2-beta2)/sehet2>1.96
gen waldhet3=abs(bw-betawy)/sehet3>1.96

gen waldclu1=abs(b1-beta1)/seclu1>1.96
gen waldclu2=abs(b2-beta2)/seclu2>1.96
gen waldclu3=abs(bw-betawy)/seclu3>1.96


mat t1=J(3,8,.) //Defining empty matrix
mat colnames t1 = "bhat_ts2sls" "stddev" "hom-se" "rob-se" "cluster-se" "Wald-hom" "Wald-rob" "Wald-clu"
mat rownames t1 = x1 x2 w

	local a=1
	foreach v of varlist b1-bw {
		qui: sum `v'
		mat t1[`a',1]=r(mean)
		mat t1[`a',2]=r(sd)
		loc ++a
		}
	local a=1
	foreach v of varlist sehom* {
		sum `v', meanonly
		mat t1[`a',3]=r(mean)
		loc ++a
		}
	local a=1
	foreach v of varlist sehet* {
		sum `v', meanonly
		mat t1[`a',4]=r(mean)
		loc ++a
		}
	local a=1
	foreach v of varlist seclu* {
		sum `v', meanonly
		mat t1[`a',5]=r(mean)
		loc ++a
		}
		local a=1
	foreach v of varlist waldhom* {
		sum `v', meanonly
		mat t1[`a',6]=r(mean)
		loc ++a
		}
		local a=1
	foreach v of varlist waldhet* {
		sum `v', meanonly
		mat t1[`a',7]=r(mean)
		loc ++a
		}
		local a=1
	foreach v of varlist waldclu* {
		sum `v', meanonly
		mat t1[`a',8]=r(mean)
		loc ++a
		}	
	


****DISPLAY RESULTS**
mat li t1
