clear

vers 14.0

****GENERATE SYNTHETIC DATA***
set mem 1g
local nobs = 3000
set obs `nobs' //create 1000 individuals
gen id=_n //Create ID numbers
set seed 12152 //set random number seed for reproducibility

cd "H:\owncloud\health\Dofiles\2s2sls"

* unobservables, indep over individ
scalar define sigma_a = 1
bys id: gen alpha1 = sigma_a*invnormal(uniform())
bys id: gen alpha2 = sigma_a*invnormal(uniform())
bys id: gen alpha3 = sigma_a*invnormal(uniform())

* generate explanatory vars, individual specific
gen w=0+3 * uniform() 
gen z1=0+1 * uniform()
gen z2=0+1 * uniform() 
gen z3=0+1 * uniform()  
gen x1_0 =1 * uniform()
gen x2_0 =1 * uniform()

* expand in dimension t
expand 3
bys id: gen t=_n


* explanatory, individual and time specific component added
scalar define sigma_expl_noise =.3
replace w = w + sigma_expl_noise * invnormal(uniform())
replace z1 = z1 + sigma_expl_noise * invnormal(uniform())
replace z2 = z2 + sigma_expl_noise * invnormal(uniform())
replace z3 = z3 + sigma_expl_noise * invnormal(uniform())
replace x1_0 = x1_0 + sigma_expl_noise * invnormal(uniform()) //part of x1 not explained through instruments
replace x2_0 = x2_0 + sigma_expl_noise * invnormal(uniform()) //part of x2 not explained through instruments


scalar define beta0 = -2
scalar define beta1 = 2
scalar define beta2 = 2
scalar define betaw = 3.0

* unobservables, indep over observ and time, 
gen u1 = 1 * invnorm(uniform()) //stand norm distr, sd=1
gen u2 = 1 * invnorm(uniform()) //stand norm distr, sd=1
gen u3 = 1 * invnorm(uniform()) //stand norm distr, sd=1


*****generate x1, x2 and y
gen x1= 0.92*z1 - 0.82*z2 + 1.41*z3+ x1_0 + u1 + alpha1
gen x2= 0.53*z1 - 3.44*z2 - 0.95*z3 + x2_0 + u2 + alpha2
gen y = x1*beta1 + x2*beta2 + w*betaw + beta0 + u3 + alpha3



***SAVE SAMPLE1
preserve  
keep y w z1 z2 z3 id
save sample1.dta, replace //2stage data set
restore


***SAVE SAMPLE2
preserve
keep x1 x2 w z1 z2 z3 id
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
qui reg y z1 z2 z3 w, rob
mat Vy1het = e(V)*e(df_r)/e(N) /*Robust variance estimate of piy1,*/
di e(df_r)
di e(df_m) 
di e(N_clust)
di e(N)
matlist e(V)
/*without degrees of freedom correction*/
qui reg y z1 z2 z3 w, cluster(id)
mat Vy1clu = e(V)*((e(df_r)-e(df_m))/e(df_r))*(e(N)-e(df_m))/e(N) /*Cluster-Robust variance estimate of piy1,*/
di e(df_r)
di e(df_m) 
di e(N_clust)
di e(N)
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
