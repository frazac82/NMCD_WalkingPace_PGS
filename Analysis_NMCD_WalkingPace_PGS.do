
************************************
****Baseline characteristics********
************************************
cls
cd "...\Updated_analysis2\"
use "db4_cad", clear

/*baseline tables*/
baselinetable   											/*
*/	age0(cts tab("p50 (p5-p95)")) 							/*
*/  msbp(cts tab("p50 (p5-p95)")  medianformat(%5.0f))		/*
*/  ldl(cts  tab("p50 (p5-p95)"))							/*
*/	tws(cts  tab("p50 (p5-p95)")) 							/*
*/	age0(cts tab("p50 (p25-p75)")) 							/*
*/  msbp(cts tab("p50 (p25-p75)") medianformat(%5.0f))		/*
*/  ldl(cts  tab("p50 (p25-p75)"))							/*
*/	tws(cts  tab("p50 (p25-p75)")) 							/*
*/	, by(sex, total) meanformat(%5.2f) sdformat(%5.1f) exportexcel("...\Results_R1\R1_Table_1", replace)

/*baseline tables by wp*/
gen wp_sex = "wp" + string(wp) + "_sex" + string(sex)
sencode wp_sex, replace gsort(wp_sex)
baselinetable   											/*
*/	age0(cts tab("p50 (p5-p95)")) 							/*
*/  msbp(cts tab("p50 (p5-p95)")  medianformat(%5.0f))		/*
*/  ldl(cts  tab("p50 (p5-p95)"))							/*
*/	tws(cts  tab("p50 (p5-p95)")) 							/*
*/	age0(cts tab("p50 (p25-p75)")) 							/*
*/  msbp(cts tab("p50 (p25-p75)") medianformat(%5.0f))		/*
*/  ldl(cts  tab("p50 (p25-p75)"))							/*
*/	tws(cts  tab("p50 (p25-p75)")) 							/*
*/	, by(wp_sex, total) meanformat(%5.2f) sdformat(%5.2f) exportexcel("...\Results_R1\R1_TableS1", replace)



*################################################################################################################################################################################################
*################################################################################################################################################################################################
*ABSOLUTE PREDICTIONS - FLEXIBLE SURVIVAL REGRESSIONS - SURVIVAL CONTINUOUS MODIFIERS - 1 OUTCOME, 1 EXPOSURES, TWO SEXES --- EFFECT ACROSS STANDARD DEVIATIONS*

clear all
cls
display "$S_TIME  $S_DATE"
cd "...\Updated_analysis2\"
use "db4_cad", clear

tab wp, m
tab wp, m nolabel
tab sex, m
tab sex, m nolabel

foreach var of varlist wp smok {
	tab `var', gen(`var')
}
drop wp wp1 smok smok1

renames sa_pgs10 sa_pgs11 sa_pgs12 sa_pgs13 sa_pgs18 sa_pgs19 sa_pgs57 sa_pgs58 sa_pgs59 \ s10 s11 s12 s13 s18 s19 s57 s58 s59  /*simpler names for loop below*/
sum s10-s59

tempfile survs_co
cap postclose stats
postfile stats sex str10 prs npart events s1 lb1 ub1 s2 lb2 ub2 s3 lb3 ub3 d2 lbd2 ubd2 d3 lbd3 ubd3 xval using `survs_co'

stset cad_time, f(cad==1) scale(365.24) id(id)
gen timevar = 10 in 1

forvalues sex = 0/1 {
	foreach prs in s10 s11 s12 s13 s18 s19 s57 s58 s59 {
		
		tabstat cad if sex == `sex', statistics(count sum) save
		qui tabstatmat mx
		local npart  = mx[1,1]
		local events = mx[2,1]
	
		cap drop int_*
		gen int_`prs'wp2 = `prs'*wp2
		gen int_`prs'wp3 = `prs'*wp3

		qui stpm2 `prs' wp2 wp3 int_`prs'wp2 int_`prs'wp3 age0 tws smok2 smok3 msbp ldl drx frx if sex == `sex', df(4) scale(hazard)
		
		forvalues k = -3(0.1)3.1 {													/*for a standard normal curve, 99.73% of the area*/
			cap drop surv1* surv2* surv3* d2* d3*

			standsurv if sex == `sex', atvars(surv1 surv2 surv3) 					/*
				*/	at1(wp2 0 wp3 0 `prs' `k' int_`prs'wp2 0   int_`prs'wp3 0)  	/*
				*/	at2(wp2 1 wp3 0 `prs' `k' int_`prs'wp2 `k' int_`prs'wp3 0)  	/*
				*/	at3(wp2 0 wp3 1 `prs' `k' int_`prs'wp2 0   int_`prs'wp3 `k') 	/*
				*/	timevar(timevar) ci contrast(difference) contrastvars(d2 d3)
			
				local s1   = surv1[1]
				local s2   = surv2[1]
				local s3   = surv3[1]
				local lb1  = surv1_lci[1]
				local lb2  = surv2_lci[1]
				local lb3  = surv3_lci[1]
				local ub1  = surv1_uci[1]
				local ub2  = surv2_uci[1]
				local ub3  = surv3_uci[1]
				local d2   = d2[1]
				local lbd2 = d2_lci[1]
				local ubd2 = d2_uci[1]
				local d3   = d3[1]
				local lbd3 = d3_lci[1]
				local ubd3 = d3_uci[1]
				
				local xval = `k'
				
				post stats (`sex') ("`prs'") (`npart') (`events') (`s1') (`lb1') (`ub1') (`s2') (`lb2') (`ub2') (`s3') (`lb3') (`ub3') (`d2') (`lbd2') (`ubd2') (`d3') (`lbd3') (`ubd3') (`xval')
		}
	}
}
postclose stats
use `survs_co', clear
compress
save "Results\absurvival", replace
display "$S_TIME  $S_DATE"


*################################################################################################################################################################################################
*################################################################################################################################################################################################
*ABSOLUTE PREDICTIONS - FLEXIBLE SURVIVAL REGRESSIONS - SURVIVAL CONTINUOUS MODIFIERS - 1 OUTCOME, 1 EXPOSURES, TWO SEXES --- EFFECT TOP 20% GENETIC RISK VS BOTTOM 80%*

clear all
cls
display "$S_TIME  $S_DATE"
cd "...\Updated_analysis2\"
use "db4_cad", clear

tab wp, m
tab wp, m nolabel
tab sex, m
tab sex, m nolabel

foreach var of varlist wp smok {
	tab `var', gen(`var')
}
drop wp wp1 smok smok1

preserve
use "...\Updated_analysis2\db4_cad_originalPRS", replace
keep id a_*
foreach var of varlist a_pgs10 a_pgs11 a_pgs12 a_pgs13 a_pgs18 a_pgs19 a_pgs57 a_pgs58 a_pgs59 {
    cumul `var', gen(c`var')
	gen cd`var' = c`var'>=0.8
}
keep id cda_*
renames cda_pgs10 cda_pgs11 cda_pgs12 cda_pgs13 cda_pgs18 cda_pgs19 cda_pgs57 cda_pgs58 cda_pgs59 \ c10 c11 c12 c13 c18 c19 c57 c58 c59
mdesc
tempfile cumdis
save `cumdis', replace
restore

merge 1:1 id using `cumdis', update
drop _merge date0 cadm_e cadi_e op_e cad_date unr_eur sa_pgs10-sa_pgs59

tempfile survs_cat
cap postclose stats
postfile stats sex str10 prs npart events s1 lb1 ub1 s2 lb2 ub2 s3 lb3 ub3 d2 lbd2 ubd2 d3 lbd3 ubd3 xval using `survs_cat'

stset cad_time, f(cad==1) scale(365.24) id(id)
gen timevar = 10 in 1

forvalues sex = 0/1 {
	foreach prs in c10 c11 c12 c13 c18 c19 c57 c58 c59 {
		
		tabstat cad if sex == `sex', statistics(count sum) save
		qui tabstatmat mx
		local npart  = mx[1,1]
		local events = mx[2,1]
	
		cap drop int_*
		gen int_`prs'wp2 = `prs'*wp2
		gen int_`prs'wp3 = `prs'*wp3

		qui stpm2 `prs' wp2 wp3 int_`prs'wp2 int_`prs'wp3 age0 tws smok2 smok3 msbp ldl drx frx if sex == `sex', df(4) scale(hazard)
		
		forvalues k = 0(1)1 {													
			cap drop surv1* surv2* surv3* d2* d3*

			standsurv if sex == `sex', atvars(surv1 surv2 surv3) 					/*
				*/	at1(wp2 0 wp3 0 `prs' `k' int_`prs'wp2 0   int_`prs'wp3 0)  	/*
				*/	at2(wp2 1 wp3 0 `prs' `k' int_`prs'wp2 `k' int_`prs'wp3 0)  	/*
				*/	at3(wp2 0 wp3 1 `prs' `k' int_`prs'wp2 0   int_`prs'wp3 `k') 	/*
				*/	timevar(timevar) ci contrast(difference) contrastvars(d2 d3)
			
				local s1   = surv1[1]
				local s2   = surv2[1]
				local s3   = surv3[1]
				local lb1  = surv1_lci[1]
				local lb2  = surv2_lci[1]
				local lb3  = surv3_lci[1]
				local ub1  = surv1_uci[1]
				local ub2  = surv2_uci[1]
				local ub3  = surv3_uci[1]
				local d2   = d2[1]
				local lbd2 = d2_lci[1]
				local ubd2 = d2_uci[1]
				local d3   = d3[1]
				local lbd3 = d3_lci[1]
				local ubd3 = d3_uci[1]
				
				local xval = `k'
				
				post stats (`sex') ("`prs'") (`npart') (`events') (`s1') (`lb1') (`ub1') (`s2') (`lb2') (`ub2') (`s3') (`lb3') (`ub3') (`d2') (`lbd2') (`ubd2') (`d3') (`lbd3') (`ubd3') (`xval')
		}
	}
}
postclose stats
use `survs_cat', clear
compress
rename xval top20
rename sex male
save "Results\absurvival_topbottom", replace
display "$S_TIME  $S_DATE"


*##################################################################################################################################################################################################
*##################################################################################################################################################################################################
*C-INDEX*

clear all
cls
display "$S_TIME  $S_DATE"
cd "...\Updated_analysis2\"
use "db4_cad", clear

renames sa_pgs10 sa_pgs11 sa_pgs12 sa_pgs13 sa_pgs18 sa_pgs19 sa_pgs57 sa_pgs58 sa_pgs59 \ s10 s11 s12 s13 s18 s19 s57 s58 s59  /*simpler names for loop below*/
gen studydes = 1
gen studyid  = 1

global adj1 = ""
forvalues sex = 0/1 {
	capture predaddc age0 tws msbp ldl drx frx i.smok i.wp if sex == `sex', studyid(studyid) subjid(id) studydes(studydes) timevar(cad_time) failure(cad==1) maxadj(1) maxonly /*
	*/saving("Results\cindex_single_sex`sex'") replace
}
forvalues sex = 0/1 {
	foreach prs in s10 s11 s12 s13 s18 s19 s57 s58 s59 {
		capture predaddc `prs' if sex == `sex', studyid(studyid) subjid(id) studydes(studydes) timevar(cad_time) failure(cad==1) maxadj(1) maxonly /*
		*/saving("Results\cindex_`prs'_sex`sex'") replace
	}
}
forvalues sex = 0/1 {
	capture predaddc age0+tws+msbp+ldl+drx+frx+i.smok if sex == `sex', studyid(studyid) subjid(id) studydes(studydes) timevar(cad_time) failure(cad==1) maxadj(1) maxonly /*
	*/saving("Results\cindex_RF_sex`sex'") replace
}

global adj2 = "age0 tws msbp ldl drx frx i.smok"
forvalues sex = 0/1 {
	predaddc i.wp if sex == `sex', studyid(studyid) subjid(id) studydes(studydes) timevar(cad_time) failure(cad==1) maxadj(2) maxonly /*
	*/saving("Results\cindex_RFWP_sex`sex'") replace
}
forvalues sex = 0/1 {
	foreach prs in s10 s11 s12 s13 s18 s19 s57 s58 s59 {
		predaddc `prs' if sex == `sex', studyid(studyid) subjid(id) studydes(studydes) timevar(cad_time) failure(cad==1) maxadj(2) maxonly /*
		*/saving("Results\cindex_RFPRS_sex`sex'_`prs'") replace
	}
}

cd "...\Updated_analysis2\Results"

///
use "cindex_single_sex0.dta", clear
gen sex = "Women"
append using "cindex_single_sex1.dta"
replace sex = "Men" if sex == ""
gen m = "Single nogenetic factors"
keep adjvars n-m
replace adjvars = "Null" if adjvars == ""
drop p_diff1-p_diff8
compress
save "cindex_singlerf.dta", replace

///
clear
foreach prs in s10 s11 s12 s13 s18 s19 s57 s58 s59 {
	append using "cindex_`prs'_sex0"
}
gen sex = "Women"
foreach prs in s10 s11 s12 s13 s18 s19 s57 s58 s59 {
	append using "cindex_`prs'_sex1"
}
replace sex = "Men" if sex == ""
gen m = "Single genetic factors"
keep adjvars n-m
compress
replace adjvars = "Null" if adjvars == ""
duplicates drop
gen prs = substr(adjvars, 1, 3) if adjvars != "Null"
replace adjvars = prs if prs != ""
drop p_diff1 prs
compress
save "cindex_singlegen.dta", replace

///
use "cindex_RF_sex0.dta", clear
gen sex = "Women"
append using "cindex_RF_sex1.dta"
replace sex = "Men" if sex == ""
gen m = "All conventional risk factors"
keep adjvars n-m
drop if adjvars == ""
drop p_diff1
compress
save "cindex_convrf.dta", replace

///
use "cindex_RFWP_sex0.dta", clear
gen sex = "Women"
append using "cindex_RFWP_sex1.dta"
replace sex = "Men" if sex == ""
gen m = "All conventional risk factors + WP"
keep adjvars n-m
drop p_diff1-diff1_ci
compress
save "cindex_convrfwp.dta", replace

///
clear
foreach prs in s10 s11 s12 s13 s18 s19 s57 s58 s59 {
	append using "cindex_RFPRS_sex0_`prs'.dta"
}
gen sex = "Women"
foreach prs in s10 s11 s12 s13 s18 s19 s57 s58 s59 {
	append using "cindex_RFPRS_sex1_`prs'.dta"
}
replace sex = "Men" if sex == ""
gen m = "All conventional risk factors + PRS"
keep adjvars n-m
drop p_diff1-diff1_ci
duplicates drop
compress
save "cindex_convrfprs.dta", replace

///
use "cindex_singlerf.dta", clear
append using "cindex_singlegen.dta"
append using "cindex_convrf.dta"
append using "cindex_convrfwp.dta"
append using "cindex_convrfprs.dta"

gen reference = "Null" if (m == "Single nogenetic factors" | m == "Single genetic factors" | m == "All conventional risk factors"), before(adjvars)
replace reference = "All conventional risk factors" if reference == ""
replace adjvars = "Null" if diff1 == .
replace adjvars = subinstr(adjvars, "age0 tws msbp ldl drx frx i.smok ", "",.)
order sex m, first
drop n nfail 
gen lb95 = cindex - 1.96*secindex
gen ub95 = cindex + 1.96*secindex
gen d_lb95 = diff1 - 1.96*sediff1
gen d_ub95 = diff1 + 1.96*sediff1
replace adjvars = "age0 tws msbp ldl drx frx i.smok" if adjvars == "age0+tws+msbp+ldl+drx+frx+i.smok"
replace adjvars = "Diabetes" 	 if adjvars == "drx"
replace adjvars = "LDL" 		 if adjvars == "ldl"
replace adjvars = "FH MI" 		 if adjvars == "frx"
replace adjvars = "Deprivation"  if adjvars == "tws"
replace adjvars = "Smoking" 	 if adjvars == "i.smok"
replace adjvars = "Walking pace" if adjvars == "i.wp"
replace adjvars = "SBP" 		 if adjvars == "msbp"
replace adjvars = "Age" 		 if adjvars == "age0"
compress
save "results_cindex.dta", replace

///C-index adding WP + PGS-13 to traditional risk factors
clear all
cls
cd "...\Updated_analysis2\"
use "db4_cad", clear

renames sa_pgs10 sa_pgs11 sa_pgs12 sa_pgs13 sa_pgs18 sa_pgs19 sa_pgs57 sa_pgs58 sa_pgs59 \ s10 s11 s12 s13 s18 s19 s57 s58 s59  /*simpler names for below*/
gen studydes = 1
gen studyid  = 1

global adj1 = "age0 tws msbp ldl drx frx i.smok"
forvalues sex = 0/1 {
	predaddc i.wp+s13 if sex == `sex', studyid(studyid) subjid(id) studydes(studydes) timevar(cad_time) failure(cad==1) maxadj(1) /*
	*/saving("Results\cindex_WP_PGS13_sex`sex'") replace
}

use "Results\cindex_WP_PGS13_sex0", clear
gen sex = "Women"
append using "Results\cindex_WP_PGS13_sex1"
replace sex = "Men" if sex == ""
keep adj0vars adjvars cindex-p_diff sex
rename adj0vars reference
gen lb95 = cindex - 1.96*secindex
gen ub95 = cindex + 1.96*secindex
gen d_lb95 = diff1 - 1.96*sediff1
gen d_ub95 = diff1 + 1.96*sediff1
drop if diff1 == .
gen m = "All conventional risk factors + PRS & PGS-13"
replace reference = "All conventional risk factors"
replace adjvars = "Walking pace + s13"
order sex m reference adjvars cindex-sediff1 lb95 ub95 d_lb95 d_ub95
drop p_diff1
save "Results\results_cindex_WP_PGS13.dta", replace

*###############################################################################################################################################################
*###############################################################################################################################################################
**RESULTS - GRAPHICS

/////DISTRIBUTION SCORES - BY OUTCOME/////
cls
cd "...\Updated_analysis2\"
use "db4_cad", clear
keep id sex wp cad sa_pgs*
reshape long sa_, i(id) j(score) string
drop wp
sort sex score
sdecode sex, replace
replace sex = "Women" if sex == "Female"
replace sex = "Men" if sex == "Male"
split score, p("s")
gen comb = sex + ", " + upper(score1) + "S" + "-" + score2
drop score1 score2
sencode comb, replace

distplot sa_, over(cad) lcolor(green red) by(comb, cols(9) note("")) subtitle(, size(medsmall) color(black) fcolor(none) lcolor(black)) xsize(5) ysize(1.5) ytitle("Cumulative probability") yla(#6, angle(h) format(%7.1f)) xla(-6(2)6, grid gmax) xtitle("Standardised score") by(, legend(off)) name("dist_prs_cad", replace)
graph save "dist_prs_cad" "...\dist_prs_cad.gph", replace
graph close _all

/////DISTRIBUTION SCORES - BY EXPOSURE/////
cls
cd "...\Updated_analysis2\"
use "db4_cad", clear
keep id sex wp cad sa_pgs*
reshape long sa_, i(id) j(score) string
sort sex score
sdecode sex, replace
replace sex = "Women" if sex == "Female"
replace sex = "Men" if sex == "Male"
split score, p("s")
gen comb = sex + ", " + upper(score1) + "S" + "-" + score2
drop score1 score2
sencode comb, replace

distplot sa_, over(wp) lcolor(blue orange green) by(comb, cols(9) note("")) subtitle(, size(medsmall) color(black) fcolor(none) lcolor(black)) xsize(5) ysize(1.5) ytitle("Cumulative probability") yla(#6, angle(h) format(%7.1f)) xla(-6(2)6, grid gmax) xtitle("Standardised score") by(, legend(off)) name("dist_prs_wp", replace)
graph save "dist_prs_wp" "...\dist_prs_wp.gph", replace
graph close _all

/////GRAPH ABSOLUTE EFFECT////
clear all
cls
use "...\absurvival", clear
mdesc

foreach var of varlist s1-ubd3 {
	replace `var' = 100*`var'
}
split prs, p("s")
replace prs = "PGS-" + upper(prs2)
drop prs1 prs2

foreach var of varlist s1-ub3 {
	gen f_`var' = 100-`var'
}
drop s1-ub3
renames f_s1 f_lb1  f_ub1   \ f1 ub1 lb1				/*inversion of CI from survival to failure - split in rows for clarity*/ 
renames f_s2 f_lb2  f_ub2   \ f2 ub2 lb2
renames f_s3 f_lb3  f_ub3   \ f3 ub3 lb3
order f1 lb1 ub1 f2 lb2 ub2 f3 lb3 ub3 d2 lbd2 ubd2 d3 lbd3 ubd3, after(xval)
tostring sex, replace
replace sex = "Women" if sex == "0"
replace sex = "Men" if sex == "1"
sencode sex, replace

twoway (line f1 xval, 		sort lcolor(blue)      lwidth(medium)) 					///
	   (line f2 xval, 		sort lcolor(orange)    lwidth(medium)) 					///
	   (line f3 xval, 		sort lcolor(green)     lwidth(medium)) 					///
	   (rarea ub1 lb1 xval, sort fcolor(blue%30)   lwidth(none))					///
	   (rarea ub2 lb2 xval, sort fcolor(orange%30) lwidth(none))					///
	   (rarea ub3 lb3 xval, sort fcolor(green%30)  lwidth(none))					///
	   if (prs == "PGS-13" | prs == "PGS-18"),										///
	   ytitle("10-year risk (%)", size(small)) ymtick(##2)							///
	   ylabel(#10, labsize(small) angle(horizontal) format(%7.0f) grid gmin gmax) 	///
	   xtitle("Standardised score", size(small))									///
	   xlabel(#10, labsize(small) format(%7.0f) grid gmin gmax) xmtick(##2) 		///
	   by(, legend(off)) by(, graphregion(fcolor(white))) by(sex prs, cols(4) note("")) subtitle(, size(small) color(black) fcolor(none) lcolor(black)) xsize(5) ysize(1.5) ///
	   name("abs_main", replace) nodraw
	   
twoway (line f1 xval, 		sort lcolor(blue)      lwidth(medium)) 					///
	   (line f2 xval, 		sort lcolor(orange)    lwidth(medium)) 					///
	   (line f3 xval, 		sort lcolor(green)     lwidth(medium)) 					///
	   (rarea ub1 lb1 xval, sort fcolor(blue%30)   lwidth(none))					///
	   (rarea ub2 lb2 xval, sort fcolor(orange%30) lwidth(none))					///
	   (rarea ub3 lb3 xval, sort fcolor(green%30)  lwidth(none))					///
	   if (prs != "PGS-13" & prs != "PGS-18"), scale(0.8)							///
	   ytitle("10-year risk (%)", size(small)) ymtick(##2)							///
	   ylabel(#10, labsize(small) angle(horizontal) format(%7.0f) grid gmin gmax) 	///
	   xtitle("Standardised score", size(small))									///
	   xlabel(#10, labsize(small) format(%7.0f) grid gmin gmax) xmtick(##2) 		///
	   by(, legend(off)) by(, graphregion(fcolor(white))) by(sex prs, rows(2) note("")) subtitle(, size(small) color(black) fcolor(none) lcolor(black)) xsize(4) ysize(2) ///
	   name("abs_others", replace)
graph save "abs_others" "...\abs_others.gph", replace	
graph close _all

/////GRAPH ABSOLUTE EFFECT DIFFERENCE////
twoway (line d2 xval, 		  sort lcolor(orange)    lwidth(medium)) 				///
	   (line d3 xval, 		  sort lcolor(green)     lwidth(medium)) 				///
	   (rarea lbd2 ubd2 xval, sort fcolor(orange%30) lwidth(none))  				///
	   (rarea lbd3 ubd3 xval, sort fcolor(green%30)  lwidth(none))  				///
	   if (prs == "PGS-13" | prs == "PGS-18"),			 							///
	   ytitle("10-year risk difference (%)", size(small)) ymtick(##2)				///
	   ylabel(#10, labsize(small) angle(horizontal) format(%7.0f) grid gmin gmax)  	///
	   xtitle("Standardised score", size(small))									///
	   xlabel(#10, labsize(small) format(%7.0f) grid gmin gmax) xmtick(##2) 		///
	   by(, legend(off)) by(, graphregion(fcolor(white))) by(sex prs, cols(4) note("")) subtitle(, size(small) color(black) fcolor(none) lcolor(black)) xsize(5) ysize(1.5) ///
	   name("absdif_main", replace) nodraw
	   
twoway (line d2 xval, 		  sort lcolor(orange)    lwidth(medium)) 				///
	   (line d3 xval, 		  sort lcolor(green)     lwidth(medium)) 				///
	   (rarea lbd2 ubd2 xval, sort fcolor(orange%30) lwidth(none))  				///
	   (rarea lbd3 ubd3 xval, sort fcolor(green%30)  lwidth(none))  				///
	   if (prs != "PGS-13" & prs != "PGS-18"), scale(0.8)							///
	   ytitle("10-year risk difference (%)", size(small)) ymtick(##2)				///
	   ylabel(#10, labsize(small) angle(horizontal) format(%7.0f) grid gmin gmax)  	///
	   xtitle("Standardised score", size(small))									///
	   xlabel(#10, labsize(small) format(%7.0f) grid gmin gmax) xmtick(##2) 		///
	   by(, legend(off)) by(, graphregion(fcolor(white))) by(sex prs, rows(2) note("")) subtitle(, size(small) color(black) fcolor(none) lcolor(black)) xsize(4) ysize(2)  ///
	   name("absdif_others", replace)
graph save "absdif_others" "...\absdif_others.gph", replace	
graph close _all 

/////GRAPH ABSOLUTE EFFECT AND DIFFERENCE - COMBINED MAIN////
graph combine abs_main absdif_main, xcommon cols(1) imargin(zero) graphregion(fcolor(white) ifcolor(white)) plotregion(margin(zero) fcolor(white) ifcolor(white)) xsize(6) ysize(4) name("main", replace) nocopies
graph save "main" "...\main.gph", replace
graph close _all

/////GRAPH CROSS-COMPARISON ABSOLUTE////
clear all
cls
use "...\absurvival", clear
mdesc
keep if prs == "s13" | prs == "s18"
replace xval = round(xval, 0.1)
keep if xval == -2 | xval == 0 | xval == 2
drop npart events d2-ubd3
reshape long s lb ub, j(wpace) i(sex prs xval)
foreach var of varlist s-ub {
	replace `var' = 100*(1-`var')
}
renames s lb ub \ f_est f_ub f_lb
order f_lb, before(f_ub)
tostring xval, replace
replace xval = "High" if xval == "2"
replace xval = "Mean" if xval == "0"
replace xval = "Low"  if xval == "-2"
tostring wpace, replace
replace wpace = "Slow"    if wpace == "1"
replace wpace = "Average" if wpace == "2"
replace wpace = "Fast"    if wpace == "3"
tostring sex, replace
replace sex = "Men"   if sex == "1"
replace sex = "Women" if sex == "0"
split prs, p("s")
replace prs = "PGS-" + upper(prs2)
drop prs1 prs2
replace prs = prs + " [metaGRS_CAD]" if prs == "PGS-18"
replace prs = prs + " [GPS_CAD]"     if prs == "PGS-13"
gen comb = sex + ", " + prs
renames xval wpace \ xval1 xval2
gen xval = xval1 + "_" + xval2
gen seq = _n
bys sex prs (seq): gen xvn = _n-1
sort seq
replace xvn = xvn + 1 if xval1 == "Mean"
replace xvn = xvn + 2 if xval1 == "High"
sencode comb, replace
labmask xvn, values(xval)

twoway (scatter f_est xvn if xval2 == "Slow",        sort mcolor(blue)   mlcolor(none) msize(small)) 		///
	   (scatter f_est xvn if xval2 == "Average",     sort mcolor(orange) mlcolor(none) msize(small)) 		///
	   (scatter f_est xvn if xval2 == "Fast",        sort mcolor(green)  mlcolor(none) msize(small)) 		///
	   (rspike  f_lb f_ub xvn if xval2 == "Slow",    sort lcolor(blue)   lwidth(thin)) 						///
	   (rspike  f_lb f_ub xvn if xval2 == "Average", sort lcolor(orange) lwidth(thin)) 						///
	   (rspike  f_lb f_ub xvn if xval2 == "Fast",    sort lcolor(green)  lwidth(thin)) 						///
	   , ytitle("10-year risk (%)") ylabel(0(2)12, labsize(small) angle(horizontal)) xtitle("") ymtick(#2)	///
	   xlabel(0 1 2 4 5 6 8 9 10, labsize(vsmall) angle(90) valuelabel) 									///
	   by(, legend(off)) xsize(6) ysize(6) by(, graphregion(fcolor(white))) 								///
	   by(comb, note("")) subtitle(, size(small) fcolor(none) lcolor(black) lwidth(vthin))					///
	   name("crosscomp", replace)
graph save "crosscomp" "...\crosscomp.gph", replace
graph close _all

/////GRAPH CROSS-COMPARISON DIFFERENCE////
clear all
cls
use "...\absurvival", clear
mdesc
keep if prs == "s13" | prs == "s18"
replace xval = round(xval, 0.1)
keep if xval == -2 | xval == 0 | xval == 2
drop npart-ub3
tostring xval, replace
replace xval = "High" if xval == "2"
replace xval = "Mean" if xval == "0"
replace xval = "Low"  if xval == "-2"
replace prs = upper(prs)
split prs, p("S")
replace prs = "PGS-" + prs2
drop prs1 prs2
replace prs = prs + " [metaGRS_CAD]" if prs == "PGS-18"
replace prs = prs + " [GPS_CAD]"     if prs == "PGS-13"
tostring sex, replace
replace sex = "Men"   if sex == "1"
replace sex = "Women" if sex == "0"
gen comb = sex + ", " + prs
foreach var of varlist sex prs xval comb {
	sencode `var', replace
}
reshape long d lbd ubd, j(diff) i(sex prs xval)
tostring diff, replace
replace diff = "Slow vs Average" if diff == "2"
replace diff = "Slow vs Fast"    if diff == "3"
foreach var of varlist d-ubd {
	replace `var' = 100*`var'
}
gen seq = _n
bys sex prs (seq): gen xvn = _n-1
replace xvn = xvn + 1 if xval == 2
replace xvn = xvn + 2 if xval == 3
sdecode xval, replace
labmask xvn, values(xval)

twoway (scatter d xvn if diff == "Slow vs Average",       sort mcolor(orange) mlcolor(none) msize(small)) 	///
	   (scatter d xvn if diff == "Slow vs Fast",          sort mcolor(green)  mlcolor(none) msize(small)) 	///
	   (rspike  lbd ubd xvn if diff == "Slow vs Average", sort lcolor(orange) lwidth(thin)) 				///
	   (rspike  lbd ubd xvn if diff == "Slow vs Fast",    sort lcolor(green)  lwidth(thin)) 				///
	   , ytitle("10-year risk difference (%)") ylabel(0(1)5, labsize(vsmall) angle(horizontal) gmin gmax) 	///
	   xtitle("") xlabel(0 1 3 4 6 7, labsize(vsmall) angle(90) valuelabel) 								///
	   by(, legend(off)) xsize(6) ysize(6) by(, graphregion(fcolor(white))) 								///
	   by(comb, note("") ixaxes iylabel ixlabel iytitle ixtitle) 											///
	   subtitle(, size(small) fcolor(none) lcolor(black) lwidth(vthin))										///
	   name("crosscomp_diff", replace)
graph save "crosscomp_diff" "...\crosscomp_diff.gph", replace
graph close _all

/////GRAPH CROSS-COMPARISON ABSOLUTE: TOP-BOTTOM 20/80%///
*aggregate data for all PGS*
clear all
cls
use "...\absurvival_topbottom", clear
mdesc
foreach var of varlist s1-ubd3 {
	replace `var' = 100*`var'
}
split prs, p("c")
replace prs = "PGS-" + upper(prs2)
drop prs1 prs2 npart events
foreach var of varlist s1-ub3 {
	gen f_`var' = 100-`var'
}
drop s1-ub3
renames f_s1 f_lb1  f_ub1   \ f1 ub1 lb1				/*inversion of CI from survival to failure - split in rows for clarity*/ 
renames f_s2 f_lb2  f_ub2   \ f2 ub2 lb2
renames f_s3 f_lb3  f_ub3   \ f3 ub3 lb3
order f1 lb1 ub1 f2 lb2 ub2 f3 lb3 ub3 d2 lbd2 ubd2 d3 lbd3 ubd3, after(top20)
tostring male, replace
replace male = "Women" if male == "0"
replace male = "Men"   if male == "1"
tostring top20, replace
replace top20 = "Top 20%"    if top20 == "1"
replace top20 = "Bottom 80%" if top20 == "0"
rename top20 PGS_group

*graph PGS-13/PGS-18
clear all
cls
use "...\absurvival_topbottom", clear
mdesc
keep if prs == "c13" | prs == "c18"
drop npart events d2-ubd3
reshape long s lb ub, j(wpace) i(male prs top20)
foreach var of varlist s-ub {
	replace `var' = 100*(1-`var')
}
renames s lb ub \ f_est f_ub f_lb
order f_lb, before(f_ub)
tostring top20, replace
replace top20 = "Top 20%" 	 if top20 == "1"
replace top20 = "Bottom 80%" if top20 == "0"
tostring wpace, replace
replace wpace = "Slow"    if wpace == "1"
replace wpace = "Average" if wpace == "2"
replace wpace = "Fast"    if wpace == "3"
tostring male, replace
replace male = "Men"   if male == "1"
replace male = "Women" if male == "0"
split prs, p("c")
replace prs = "PGS-" + upper(prs2)
drop prs1 prs2
replace prs = prs + " [metaGRS_CAD]" if prs == "PGS-18"
replace prs = prs + " [GPS_CAD]"     if prs == "PGS-13"
gen comb = male + ", " + prs
renames top20 wpace \ xval1 xval2
gen top20 = xval1 + "_" + xval2
gen seq = _n
bys male prs (seq): gen xvn = _n-1
sort seq
replace xvn = xvn + 1 if xval1 == "Top 20%"
sencode comb, replace
labmask xvn, values(top20)

twoway (scatter f_est xvn if xval2 == "Slow",        sort mcolor(blue)   mlcolor(none) msize(small)) 		///
	   (scatter f_est xvn if xval2 == "Average",     sort mcolor(orange) mlcolor(none) msize(small)) 		///
	   (scatter f_est xvn if xval2 == "Fast",        sort mcolor(green)  mlcolor(none) msize(small)) 		///
	   (rspike  f_lb f_ub xvn if xval2 == "Slow",    sort lcolor(blue)   lwidth(thin)) 						///
	   (rspike  f_lb f_ub xvn if xval2 == "Average", sort lcolor(orange) lwidth(thin)) 						///
	   (rspike  f_lb f_ub xvn if xval2 == "Fast",    sort lcolor(green)  lwidth(thin)) 						///
	   , ytitle("10-year risk (%)") ylabel(0(2)12, labsize(small) angle(hor) gmax) xtitle("") ymtick(#2)	///
	   xlabel(0 1 2 4 5 6, labsize(vsmall) angle(90) valuelabel) 											///
	   by(, legend(off)) xsize(6) ysize(6) by(, graphregion(fcolor(white))) 								///
	   by(comb, note("")) subtitle(, size(small) fcolor(none) lcolor(black) lwidth(vthin))					///
	   name("crosscomp_topbottom", replace)
graph save "crosscomp_topbottom" "...\crosscomp_topbottom.gph", replace
graph close _all

/////GRAPH CROSS-COMPARISON DIFFERENCE: TOP-BOTTOM 20/80%///
clear all
cls
use "...\absurvival_topbottom", clear
mdesc
keep if prs == "c13" | prs == "c18"
drop npart-ub3
tostring top20, replace
replace top20 = "Top 20%" 	 if top20 == "1"
replace top20 = "Bottom 80%" if top20 == "0"
replace prs = upper(prs)
split prs, p("C")
replace prs = "PGS-" + prs2
drop prs1 prs2
replace prs = prs + " [metaGRS_CAD]" if prs == "PGS-18"
replace prs = prs + " [GPS_CAD]"     if prs == "PGS-13"
tostring male, replace
replace male = "Men"   if male == "1"
replace male = "Women" if male == "0"
gen comb = male + ", " + prs
foreach var of varlist male prs top20 comb {
	sencode `var', replace
}
reshape long d lbd ubd, j(diff) i(male prs top20)
tostring diff, replace
replace diff = "Slow vs Average" if diff == "2"
replace diff = "Slow vs Fast"    if diff == "3"
foreach var of varlist d-ubd {
	replace `var' = 100*`var'
}
gen seq = _n
bys male prs (seq): gen xvn = _n-1
replace xvn = xvn + 1 if top20 == 2
sdecode top20, replace
labmask xvn, values(top20)

twoway (scatter d xvn if diff == "Slow vs Average",       sort mcolor(orange) mlcolor(none) msize(small)) 	///
	   (scatter d xvn if diff == "Slow vs Fast",          sort mcolor(green)  mlcolor(none) msize(small)) 	///
	   (rspike  lbd ubd xvn if diff == "Slow vs Average", sort lcolor(orange) lwidth(thin)) 				///
	   (rspike  lbd ubd xvn if diff == "Slow vs Fast",    sort lcolor(green)  lwidth(thin)) 				///
	   , ytitle("10-year risk difference (%)") ylabel(0(1)5, labsize(vsmall) angle(horizontal) gmin gmax) 	///
	   xtitle("") xlabel(0 1 3 4, labsize(vsmall) angle(90) valuelabel) 									///
	   by(, legend(off)) xsize(6) ysize(6) by(, graphregion(fcolor(white))) 								///
	   by(comb, note("") ixaxes iylabel ixlabel iytitle ixtitle) 											///
	   subtitle(, size(small) fcolor(none) lcolor(black) lwidth(vthin))										///
	   name("crosscomp_diff_topbottom", replace)
graph save "crosscomp_diff_topbottom" "...\crosscomp_diff_topbottom.gph", replace
graph close _all

/////GRAPH C-INDEX COMBINED///
clear all
cls
use "...\results_cindex.dta", clear
append using "Results\results_cindex_WP_PGS13.dta"
keep if reference == "All conventional risk factors"
sort sex adjvars
duplicates drop sex adjvars, force
replace adjvars = "NGRF" if adjvars == "Null"
replace adjvars = "WP" if adjvars == "Walking pace"
split adjvars, p("s")
replace adjvars = "PGS-" + adjvars2 if adjvars1 == ""
drop adjvars1 adjvars2
replace adjvars = "WP + " + "PGS-13" if adjvars == "Walking pace + s13"
sort sex adjvars
compress
bys sex (adjvars): gen seq1 = _n
replace seq1 = seq1 + 1 if seq1>=2
gen ord1 = 1 if seq1 == 1
replace ord1 = 2 if seq1>1 & seq1<=11
replace ord1 = 3 if seq1 == 12
replace ord1 = 4 if seq1 == 13
gsort -sex ord1 seq1
sencode sex, replace
replace adjvars = "+ " + adjvars if adjvars != "NGRF"
labmask seq1, values(adjvars)
sum lb95 ub95

twoway (scatter cindex seq1 if ord1 == 1, sort mcolor(black)    msymbol(circle))  (rspike lb95 ub95 seq1 if ord1 == 1, sort lcolor(black)    lwidth(vthin)) ///
	   (scatter cindex seq1 if ord1 == 2, sort mcolor(blue)     msymbol(circle))  (rspike lb95 ub95 seq1 if ord1 == 2, sort lcolor(blue)     lwidth(vthin)) ///
	   (scatter cindex seq1 if ord1 == 3, sort mcolor(black)    msymbol(square))  (rspike lb95 ub95 seq1 if ord1 == 3, sort lcolor(black)    lwidth(vthin)) ///
	   (scatter cindex seq1 if ord1 == 4, sort mcolor(cranberry) msymbol(diamond)) (rspike lb95 ub95 seq1 if ord1 == 4, sort lcolor(cranberry) lwidth(vthin)) ///
		, xtitle("") xlabel(1(1)13, grid valuelabel ang(v) labsize(small)) 																	 		 ///
		ytitle("C-index") ylabel(0.68(0.02)0.82, angle(horizontal) format(%7.2f) grid gmin gmax) ymtick(##2) 										 ///
		by(, legend(off)) by(, graphregion(fcolor(white))) by(sex, cols(2) note("")) subtitle(, size(small) color(black) fcolor(none) lcolor(black)) ///
		name("cindex_comb", replace) xsize(5.5) ysize(3.5)
graph save "cindex_comb" "...\cindex_comb.gph", replace
graph close _all

/////GRAPH C-INDEX SINGLES///
clear all
cls
use "...\results_cindex.dta", clear

keep if reference != "All conventional risk factors"
drop if m == "All conventional risk factors"
drop if adjvars == "Null"
drop reference
compress
gsort -sex -m cindex
split adjvars, p("s")
replace adjvars = "PGS-" + adjvars2 if m == "Single genetic factors"
drop adjvars1 adjvars2
sum lb95 ub95

foreach sex in Women Men {
	preserve
	keep if sex == "`sex'"
	gen seq = _n
	replace seq = seq + 1 if m == "Single genetic factors"
	gen ord1 = 1 if m == "Single nogenetic factors"
	replace ord1 = 2 if m == "Single genetic factors"
	labmask seq, values(adjvars)
	twoway (scatter cindex seq if ord1 == 1, sort mcolor(black) msymbol(circle)) (rspike lb95 ub95 seq if ord1 == 1, sort lcolor(black) lwidth(vthin))   ///
		   (scatter cindex seq if ord1 == 2, sort mcolor(blue)  msymbol(circle)) (rspike lb95 ub95 seq if ord1 == 2, sort lcolor(blue)  lwidth(vthin))   ///
		    , xtitle("") xlabel(1(1)18, grid valuelabel ang(v) labsize(small)) 																 		 	 ///
		    ytitle("C-index") ylabel(0.50(0.02)0.7, angle(horizontal) format(%7.2f) grid gmin gmax) ymtick(##2)											 ///
		    by(, legend(off)) by(, graphregion(fcolor(white))) by(sex, cols(2) note("")) subtitle(, size(small) color(black) fcolor(none) lcolor(black)) ///
		    name("cs_`sex'", replace) xsize(5.5) ysize(3.5) nodraw
	restore
}
graph combine cs_Women cs_Men, ycommon cols(2) imargin(zero) graphregion(fcolor(white) ifcolor(white)) plotregion(margin(zero) fcolor(white) ifcolor(white)) xsize(6) ysize(3) name("cindex_singles", replace) nocopies
graph save "cindex_singles" "...\cindex_singles.gph", replace
graph close _all

******************************
***Graphs with PGS-18*********
******************************
clear all
cls
use "...\absurvival", clear
mdesc

foreach var of varlist s1-ubd3 {
	replace `var' = 100*`var'
}
split prs, p("s")
replace prs = "PGS-" + upper(prs2)
drop prs1 prs2

foreach var of varlist s1-ub3 {
	gen f_`var' = 100-`var'
}
drop s1-ub3
renames f_s1 f_lb1  f_ub1   \ f1 ub1 lb1				/*inversion of CI from survival to failure - split in rows for clarity*/ 
renames f_s2 f_lb2  f_ub2   \ f2 ub2 lb2
renames f_s3 f_lb3  f_ub3   \ f3 ub3 lb3
order f1 lb1 ub1 f2 lb2 ub2 f3 lb3 ub3 d2 lbd2 ubd2 d3 lbd3 ubd3, after(xval)
tostring sex, replace
replace sex = "Women" if sex == "0"
replace sex = "Men" if sex == "1"
sencode sex, replace

twoway (line f1 xval, 		sort lcolor(blue)      lwidth(medium)) 					///
	   (line f2 xval, 		sort lcolor(orange)    lwidth(medium)) 					///
	   (line f3 xval, 		sort lcolor(green)     lwidth(medium)) 					///
	   (rarea ub1 lb1 xval, sort fcolor(blue%30)   lwidth(none))					///
	   (rarea ub2 lb2 xval, sort fcolor(orange%30) lwidth(none))					///
	   (rarea ub3 lb3 xval, sort fcolor(green%30)  lwidth(none))					///
	   if (prs != "PGS-13"), scale(0.8)												///
	   ytitle("10-year risk (%)", size(small)) ymtick(##2)							///
	   ylabel(#10, labsize(small) angle(horizontal) format(%7.0f) grid gmin gmax) 	///
	   xtitle("Standardised score", size(small))									///
	   xlabel(#10, labsize(small) format(%7.0f) grid gmin gmax) xmtick(##2) 		///
	   by(, legend(off)) by(, graphregion(fcolor(white))) by(sex prs, rows(2) note("")) subtitle(, size(small) color(black) fcolor(none) lcolor(black)) xsize(4) ysize(2) ///
	   name("R1_abs_others_FigS4", replace)
graph save "R1_abs_others_FigS4" "...\Results_R1\R1_abs_others_FigS4.gph", replace	

twoway (line d2 xval, 		  sort lcolor(orange)    lwidth(medium)) 				///
	   (line d3 xval, 		  sort lcolor(green)     lwidth(medium)) 				///
	   (rarea lbd2 ubd2 xval, sort fcolor(orange%30) lwidth(none))  				///
	   (rarea lbd3 ubd3 xval, sort fcolor(green%30)  lwidth(none))  				///
	   if (prs != "PGS-13"), scale(0.8)												///
	   ytitle("10-year risk difference (%)", size(small)) ymtick(##2)				///
	   ylabel(#10, labsize(small) angle(horizontal) format(%7.0f) grid gmin gmax)  	///
	   xtitle("Standardised score", size(small))									///
	   xlabel(#10, labsize(small) format(%7.0f) grid gmin gmax) xmtick(##2) 		///
	   by(, legend(off)) by(, graphregion(fcolor(white))) by(sex prs, rows(2) note("")) subtitle(, size(small) color(black) fcolor(none) lcolor(black)) xsize(4) ysize(2)  ///
	   name("R1_absdif_others_FigS5", replace)
graph save "R1_absdif_others_FigS5" "...\Results_R1\R1_absdif_others_FigS5.gph", replace	
graph close _all 


************************************
******C-index by sex****************
************************************
clear all
cls
use "...\results_cindex.dta", clear
append using "...\results_cindex_WP_PGS13.dta"
keep if reference == "All conventional risk factors"
sort sex adjvars
duplicates drop sex adjvars, force
replace adjvars = "NGRF" if adjvars == "Null"
replace adjvars = "WP" if adjvars == "Walking pace"
split adjvars, p("s")
replace adjvars = "PGS-" + adjvars2 if adjvars1 == ""
drop adjvars1 adjvars2
replace adjvars = "WP + " + "PGS-13" if adjvars == "Walking pace + s13"
sort sex adjvars
compress
drop cindex secindex lb95 ub95 m
renames diff1 sediff1 d_lb95 d_ub95 \ d_ sed_ lbd_ ubd_
reshape wide d-ubd, i(adjvars) j(sex) string
order reference adjvars
gen mw = d_Men - d_Women
gen se_mw = sqrt(sed_Men^2 + sed_Women^2)
gen lb_mw = mw - 1.96*se_mw 
gen ub_mw = mw + 1.96*se_mw


*************************************************
***Formal interactions***************************
*************************************************
clear all
cls

tempfile betas
save `betas', emptyok replace

tempfile lrts
save `lrts', emptyok replace

cd "...\Updated_analysis2\"
use "db4_cad", clear

renames sa_pgs10 sa_pgs11 sa_pgs12 sa_pgs13 sa_pgs18 sa_pgs19 sa_pgs57 sa_pgs58 sa_pgs59 \ s10 s11 s12 s13 s18 s19 s57 s58 s59  /*simpler names for loop below*/
stset cad_time, f(cad==1) scale(365.24) id(id)

cls
forvalues sex = 0/1 {
	foreach prs in s10 s11 s12 s13 s18 s19 s57 s58 s59 {
		
		qui stpm2 c.`prs'##i.wp age0 tws i.smok msbp ldl drx frx if sex == `sex', df(4) scale(hazard)
		qui estimate store inter
		preserve
		qui parmest, fast eform
		qui keep if strpos(parm, "#")>1
		qui gen sex = `sex'
		qui gen prs = "`prs'"
		qui append using `betas'
		qui save `betas', replace
		restore
		
		qui stpm2 c.`prs' i.wp age0 tws i.smok msbp ldl drx frx if sex == `sex', df(4) scale(hazard)
		qui estimate store nointer
		
		qui lrtest (inter) (nointer)
		qui local p_int = `r(p)'
		preserve
		qui clear
		qui set obs 1
		qui gen sex = `sex'
		qui gen prs = "`prs'"
		qui gen pvalue_lrtest = `p_int'
		qui append using `lrts'
		qui save `lrts', replace
		restore
				
		di "Sex = `sex' | PGS = `prs' | p(lrt) = `p_int'"
	}
}		

use `betas', clear
merge m:1 sex prs using `lrts'
drop if p == .
drop eq z _merge
sort sex prs parm
split parm, p(#)
replace parm = "Average" if parm1 == "2.wp"
replace parm = "Brisk" if parm1 == "3.wp"
replace prs = subinstr(prs, "s", "PGS-",.)
drop parm1 parm2 p stderr
foreach var of varlist estimate min95 max95 {
	tostring `var', replace format(%04.2f) force
}
tostring pvalue_lrtest, replace format(%04.3f) force
tostring sex, replace
replace sex = "Women" if sex == "0"
replace sex = "Men" if sex == "1"
gen Hazard_Ratio = estimate + " (" + min95 + ", " + max95 + ")"
drop estimate-max95
order sex prs parm Hazard_Ratio
renames sex prs parm pvalue_lrtest \ Sex PGS Walking_pace p_valueLRT
bys Sex PGS (Walking_pace): gen seq = _n
replace p_valueLRT = "" if seq == 2
gsort -Sex PGS Walking_pace
drop seq 
compress
export delimited using "...\Table_S3.csv", datafmt replace


*************************************************
***Difference of difference**********************
*************************************************
clear all
cls

cd "...\Updated_analysis2\"
use "db4_cad", clear

tab wp, m
tab wp, m nolabel
tab sex, m
tab sex, m nolabel

foreach var of varlist wp smok {
	tab `var', gen(`var')
}
drop wp wp1 smok smok1

preserve
use "...\db4_cad_originalPRS", replace
keep id a_*
foreach var of varlist a_pgs10 a_pgs11 a_pgs12 a_pgs13 a_pgs18 a_pgs19 a_pgs57 a_pgs58 a_pgs59 {
    cumul `var', gen(c`var')
	gen cd`var' = c`var'>=0.8
}
keep id cda_*
renames cda_pgs10 cda_pgs11 cda_pgs12 cda_pgs13 cda_pgs18 cda_pgs19 cda_pgs57 cda_pgs58 cda_pgs59 \ c10 c11 c12 c13 c18 c19 c57 c58 c59
mdesc
tempfile cumdis
save `cumdis', replace
restore

merge 1:1 id using `cumdis', update
drop _merge date0 cadm_e cadi_e op_e cad_date unr_eur sa_pgs10-sa_pgs59

stset cad_time, f(cad==1) scale(365.24) id(id)
gen timevar = 10 in 1

gen int_c13wp2 = c13*wp2
gen int_c13wp3 = c13*wp3

preserve
clear
tempfile results
save `results', emptyok replace
restore

mata:
	function resmata(at) {
		return((at[4]-at[3]) - (at[2]-at[1]))	/* see below what these [at] are */
	}
end

forvalues sex = 0/1 {
		preserve
		qui stpm2 c13 wp2 wp3 int_c13wp2 int_c13wp3 age0 tws smok2 smok3 msbp ldl drx frx if sex == `sex', df(4) scale(hazard)
		
		standsurv if sex == `sex', atvars(surv1 surv2 surv3 surv4) 		/*
			*/	at1(wp2 0 wp3 0 c13 0 int_c13wp2 0 int_c13wp3 0)  		/* Low-Moderate genetic risk, slow
			*/	at2(wp2 0 wp3 1 c13 0 int_c13wp2 0 int_c13wp3 0) 		/* Low-Moderate genetic risk, fast
			*/	at3(wp2 0 wp3 0 c13 1 int_c13wp2 0 int_c13wp3 0)  		/* High genetic risk, slow
			*/	at4(wp2 0 wp3 1 c13 1 int_c13wp2 0 int_c13wp3 1) 		/* High genetic risk, fast
			*/	timevar(timevar) ci contrast(difference) userfunction(resmata)
			
		keep if timevar != .
		keep surv* _contrast4* _userfunc*
		gen sex = `sex'
		append using `results'
		save `results', replace
		restore
}

use `results', clear
gen high_slfa = 100*((1-surv3) - (1-surv4))
gen low_slfa  = 100*((1-surv1) - (1-surv2))
drop surv2* surv3*
renames _userfunc _userfunc_lci _userfunc_uci \ dd dd_lb dd_ub
foreach var of varlist dd dd_lb dd_ub {
	replace `var' = `var'*100
}
foreach var of varlist surv* {
	replace `var' = 100*(1-`var')
}
renames surv1_lci surv1_uci surv4_lci surv4_uci \ surv1_ub surv1_lb surv4_ub surv4_lb
foreach var of varlist _contrast* {
	replace `var' = 100*(-`var')
}
renames _contrast4_1 _contrast4_1_lci _contrast4_1_uci \ d41 d41_ub d41_lb
tostring sex, replace
replace sex = "Women" if sex == "0"
replace sex = "Men" if sex == "1" 
order sex surv1 surv1_lb surv1_ub surv4 surv4_lb surv4_ub d41 d41_lb d41_ub high_slfa low_slfa dd dd_lb dd_ub  
export delimited using "...\Diff_of_diff.csv", datafmt replace
