%macro csi_index(
	data_s=,
	data_w=,
	out=csi_out,
	wby=pri_site sex,
	oby=Year_intraval range,
	s=ASS,
	var=Var_ASS,
	wvar=total_count,
	count=TOT_CASES
);
	%local alpha eps ycap _wmiss _dup;
	%let alpha = 0.05;
	%let eps = 1e-6;
	%let ycap = 700;
	%let _wmiss=0;
	%let _dup=0;

	/*---------------------------------------
		1) Data Setup: _csi_base
	---------------------------------------*/
	data _csi_base;
		set &data_s(keep=&wby &oby &s &var &count);

		_s = &s;
		_var = &var;
		if _s < 0 then _s = .;
		if _var < 0 then _var = .;
		_avail = (nmiss(_s,_var)=0);

		length csi_gp $50;
		csi_gp = strip(pri_site); output;
		csi_gp = '000_ALL'; output;

		keep &wby &oby csi_gp _s _var _avail &count;
		run;
	/*---------------------------------------
		2) Aggregate Number of Cases: _csi_totcase
	---------------------------------------*/
	proc sort data=_csi_base; by &oby csi_gp;
	proc freq data=_csi_base noprint;
		by &oby;
		tables csi_gp / out=_csi_totcase(drop=percent rename=(count=TOT_CASES));
		weight &count;
		run;

	/*---------------------------------------
		3) Merge External Weights: _csi_avail_w
	---------------------------------------*/
	proc sort data=_csi_base; by &wby;
	proc sort data=&data_w; by &wby;
	data _csi_avail_w;
		merge _csi_base(in=ina keep=&wby &oby csi_gp _avail where=(_avail=1))
					&data_w (in=inw keep=&wby &wvar rename=(&wvar=_wgt));
		by &wby;
		if ina;
		if missing(_wgt) or _wgt<=0 then call symputx('_wmiss','1');
	run;
	proc sort data=_csi_avail_w; by &oby &wby csi_gp;
	data _csi_dup;set _csi_avail_w;
		by &oby &wby csi_gp;
		if not (first.csi_gp and last.csi_gp) then call symputx('_dup','1');
		run;
	%if %sysevalf(%superq(_wmiss)=1 or %superq(_dup)=1) %then %do;
		%if %superq(_wmiss)=1 %then
			%put ERROR: (age_std_survival) missing or incorrect weight (_wgt);
		%if %superq(_dup)=1 %then
			%put ERROR: (&data_s or &data_w) has duplicate rows in (&oby &wby csi_gp);
		%abort cancel;
	%end;

	/*---------------------------------------
		4) Weight Normalization: _csi_wgt
	---------------------------------------*/
	proc sort data=_csi_avail_w; by &oby csi_gp &wby;
	proc summary data=_csi_avail_w nway;
		class &oby csi_gp;
		var _wgt;
		output out=_csi_wgt_den(drop=_type_ _freq_) sum=w_den;
	data _csi_wgt;
		merge _csi_avail_w(in=a) _csi_wgt_den(in=b);
		by &oby csi_gp;
		if a;
		_wgt_norm = _wgt / w_den;
		keep &oby csi_gp &wby _wgt_norm;
		run;

	/*---------------------------------------
		5) Merge Weights Back: _csi_merged
	---------------------------------------*/
	proc sort data=_csi_base; by &oby csi_gp &wby;
	proc sort data=_csi_wgt; by &oby csi_gp &wby;
	data _csi_merged;
		merge _csi_base(in=inb)
					_csi_wgt (in=inw keep=&oby csi_gp &wby _wgt_norm);
		by &oby csi_gp &wby;
		if inb and inw;
	run;

	/*---------------------------------------
		6) Calculate CSI and Confidence Intervals: _csi_calc
	---------------------------------------*/
	data _csi_calc;
		set _csi_merged;
		by &oby csi_gp;

		retain _sumS _varS _z n_used;
		if _n_=1 then _z = probit(1-&alpha./2);

		if first.csi_gp then do;
			_sumS=0; _varS=0; n_used=0;
		end;

		n_used + 1;
		_sumS + _wgt_norm * _s;
		_varS + _wgt_norm*_wgt_norm * _var;

		if last.csi_gp then do;
			length CSI Var_CSI CSI_LCL CSI_UCL CSI_LCL_loglog CSI_UCL_loglog 8;

			CSI = _sumS;
			Var_CSI	= _varS;
			/* 95% CI (linear) */
			CSI_LCL = CSI - _z*sqrt(Var_CSI);
			CSI_UCL = CSI + _z*sqrt(Var_CSI);
			/* 95% CI (cloglog) */
			CSI_safe = min(max(CSI, &eps), 1-&eps);
			CSI_LCL_loglog = max(&eps, exp(-exp(min(log(-log(CSI_safe)) + _z*sqrt(Var_CSI)/(CSI_safe*abs(log(CSI_safe))), &ycap))));
			CSI_UCL_loglog = min(1-&eps, exp(-exp(max(log(-log(CSI_safe)) - _z*sqrt(Var_CSI)/(CSI_safe*abs(log(CSI_safe))), -&ycap))));

			/* Run these 3 lines if the results need to be restricted to the 0-1 range */
/*			CSI = min(max(CSI, 0), 1);*/
/*			CSI_LCL = min(max(CSI_LCL, 0), 1);*/
/*			CSI_UCL = min(max(CSI_UCL, 0), 1);*/

			keep &oby csi_gp n_used CSI Var_CSI CSI_LCL CSI_UCL CSI_LCL_loglog CSI_UCL_loglog;
			output;
			end;
		run;

	/*---------------------------------------
		7) Merge Total Cases Back and Output Final Dataset: &out
	---------------------------------------*/
	proc sort data=_csi_totcase; by &oby csi_gp;
	proc sort data=_csi_calc; by &oby csi_gp;
	data &out;
		merge _csi_totcase(in=int) _csi_calc;
		by &oby csi_gp;
		if int;
		if missing(n_used) then n_used=0;
		attrib _all_ label='';
		run;

	/*---------------------------------------
		8) Clean Up Temporary Datasets
	---------------------------------------*/
	proc datasets lib=work nolist;
		delete _csi_base _csi_totcase _csi_avail_w _csi_wgt_den _csi_dup _csi_wgt _csi_merged _csi_calc;
		quit;
%mend csi_index;



libname spc 'D:\path_to_project\output\spc_bootstrap_5000';

data rs_1yr; set spc.rs_1yr_final_:;
data rs_5yr; set spc.rs_5yr_final_:;

proc sort data=rs_1yr; by cohort_1yr pri_site bootstrap;
data dup_1yr;
	set rs_1yr;
	by cohort_1yr pri_site bootstrap;
	if not (first.bootstrap and last.bootstrap);
	run;*0;

proc sort data=rs_5yr; by cohort_5yr pri_site bootstrap;
data dup_5yr;
	set rs_5yr;
	by cohort_5yr pri_site bootstrap;
	if not (first.bootstrap and last.bootstrap);
	run;*0;



data rs_5000_1;
	format Year_intraval $9. Year_gap $1. pri_site $50. sex $1. bootstrap $5.;
	set rs_1yr;
	Year_intraval = put(cohort_1yr,z4.);
	Year_gap = '1';
	sex = substr(pri_site, length(pri_site),1);
	pri_site = substr(pri_site, 1, length(pri_site)-2);
	drop pop no_miss1 no_miss2 run_time cohort_1yr; 
data rs_5000_2;
	format Year_intraval $9. Year_gap $1. pri_site $50. sex $1. bootstrap $5.;
	set rs_5yr;
	Year_intraval = cohort_5yr;
	Year_gap = '5';
	sex = substr(pri_site, length(pri_site),1);
	pri_site = substr(pri_site, 1, length(pri_site)-2);
	rename cr_weight = cr;
	drop pop no_miss1 no_miss2 run_time cohort_5yr; 
	run;

data rs_5000;
	set rs_5000_1 rs_5000_2;
	ASS = cr;
	Var_ASS = 0;
	TOT_CASES = 1;
	drop cr;
	run;
proc sort data=rs_5000;
	by Year_intraval Year_gap pri_site sex bootstrap;
data dup;
	set rs_5000;
	by Year_intraval Year_gap pri_site sex bootstrap;
	if not (first.bootstrap and last.bootstrap);
	run;
proc freq data=rs_5000 noprint;table pri_site / out=A1;run;
proc freq data=rs_5000 noprint;table sex / out=A1;run;
proc freq data=rs_5000 noprint;table bootstrap / out=A1;run;
proc freq data=rs_5000 noprint;table Year_intraval / out=A1;run;




libname weight 'D:\path_to_project\data';
data rs_data;set weight.rs_data;run;

data rs_data2;
	set rs_data;
	if sex in ('M1' 'M2' 'S1' 'S2');
	if behavior in ('M3' 'S3');
	*pri_site = substr(pri_site,1,length(pri_site)-2);
	sex = substr(sex,2,1);
	gpage = gp;
	keep Year_intraval Year_gap pri_site sex behavior gpage no range l d w p_star cp cr se_cp;
	run;

data weight_1;
	set rs_data2;
	if Year_gap='1';
	if '100' <= year_intraval <= '109';
	run;

proc sort data=weight_1 nodupkey;by pri_site sex gpage Year_intraval;
proc freq data=weight_1 noprint;
	by pri_site sex;
	table gpage / out = weight_as;
	weight no;
	run;
data weight_as;set weight_as;percent = percent / 100;run;
proc transpose data=weight_as out=weight_as_table(drop=_NAME_ _LABEL_);
	by pri_site;
	id sex gpage;
	var percent;
	run;


proc sort data=weight_1 nodupkey;by pri_site sex gpage Year_intraval;
proc freq data=weight_1 noprint;
	by pri_site;
	table sex / out = weight_csi;
	weight no;
	run;
data weight_csi;set weight_csi;percent = percent / 100;run;

proc transpose data=weight_csi out=weight_csi_table(drop=_NAME_ _LABEL_);
	by pri_site;
	id sex;
	var percent;
	run;



%csi_index(
	data_s=rs_5000,
	data_w=weight_csi,
	out=csi_out,
	wby=pri_site sex,
	oby= Year_intraval Year_gap bootstrap,
	s=ASS,
	var=Var_ASS,
	wvar = COUNT,
	count  = TOT_CASES
);



proc sort data=csi_out;by csi_gp Year_intraval Year_gap bootstrap;
proc univariate data=csi_out noprint;
	by csi_gp Year_intraval Year_gap;
	var csi;
	output out=csi_ci_1 n=bootstrap_no mean=mean pctlpre=CI_  pctlpts=2.5,50,97.5; 
	run;
data csi_ci;
	format csi_gp;
	set csi_ci_1;
	rename csi_gp = pri_site;
	run;

proc sort data=csi_ci;by Year_gap pri_site Year_intraval;run;	
proc export data=csi_ci
	outfile = 'D:\path_to_project\output\Cancer_Survival_Index(ASRS-SPC-CSI).xls'
	dbms = excel
	replace;
	run;
