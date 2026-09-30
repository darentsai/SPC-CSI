ods _all_ close;
options NOMPRINT;

%macro bootstrap(in_data, bootstrap_no);

	%do value=1 %to &bootstrap_no;
		%let run_time = %sysfunc(datetime());
		%put bootstrap_no = &value--- %sysfunc(putn(&run_time., E8601DT26.3));

		%if &value. = 1 %then %do;
			data rs_1yr_final; length bootstrap $5. cohort_1yr 8 pop 8 cr 8 no_miss1 3 no_miss2 3 run_time 8; format run_time e8601dt19.0;stop;
			data rs_5yr_final; length bootstrap $5. cohort_5yr $9. cr_weight 8 run_time 8; format run_time e8601dt19.0;stop; run;
			%end;
		%else %do;
			proc datasets noprint;save rs_data2 filelist &in_data rs_1yr_final rs_5yr_final;quit;
			%end;

		proc sql noprint;select min(cohort) into : first_cohort from &in_data;quit;
		data spc1;
			set &in_data;
			c = cohort - %eval(&first_cohort-1);
			p = period - %eval(&first_cohort);
			s = period - cohort;
			keep c p s pop ass var_ass;
			run;

		%let first_c=; %let last_c=; %let ctr_c=; %let ref_c=; %let c_list=; %let c_list2=;
		%let first_p=; %let last_p=; %let ctr_p=; %let ref_p=; %let p_list=; %let p_list2=;
		%let first_s=; %let last_s=; %let ctr_s=; %let ref_s=; %let s_list=; %let s_list2=;
		proc sql noprint;
			select min(c), max(c), int(median(c)) into : first_c trimmed, :last_c trimmed, :ref_c trimmed from (select distinct c from spc1);
			select min(p), max(p), int(median(p)) into : first_p trimmed, :last_p trimmed, :ref_p trimmed from (select distinct p from spc1);
			select min(s), max(s), int(median(s)) into : first_s trimmed, :last_s trimmed, :ref_s trimmed from (select distinct s from spc1);
			select 'c'||put(c,3. -L) into : c_list separated by ' ' from (select distinct c from spc1) order by c;
			select 'p'||put(p,3. -L) into : p_list separated by ' ' from (select distinct p from spc1) order by p;
			select 's'||put(s,3. -L) into : s_list separated by ' ' from (select distinct s from spc1) order by s;
			select 'c'||put(c,3. -L) into : c_list2 separated by ' ' from (select distinct c from spc1 where c ^= %eval(&ref_c)) order by c;
			select 'p'||put(p,3. -L) into : p_list2 separated by ' ' from (select distinct p from spc1 where p ^= %eval(&first_p) and p ^= %eval(&ref_p)) order by p;
			select 's'||put(s,3. -L) into : s_list2 separated by ' ' from (select distinct s from spc1 where s ^= %eval(&ref_s)) order by s;
			quit;
		%let ctr_s = %sysevalf((&first_s + &last_s)/2);
		%let ctr_c = %sysevalf((&first_c + &last_c)/2);
		%let ctr_p = %sysevalf((&first_p + &last_p)/2);
/*		%put first_cohort = &first_cohort;*/
/*		%put cohort : first-ctr/ref-last = &first_c - &ctr_c/&ref_c - &last_c;*/
/*		%put period : first-ctr/ref-last = &first_p - &ctr_p/&ref_p - &last_p;*/
/*		%put survship : first-ctr/ref-last = &first_s - &ctr_s/&ref_s - &last_s;*/
/*		%put c_list = &c_list;*/
/*		%put c_list2 = &c_list2;*/
/*		%put p_list = &p_list;*/
/*		%put p_list2 = &p_list2;*/
/*		%put s_list = &s_list;*/
/*		%put s_list2 = &s_list2;*/

		data spc2;
			set spc1;
			if _N_=1 then call streaminit(12345 + (&value - 1));
			if missing(ass) or ass < 0 or missing(var_ass) or var_ass < 0 then do;
				log_log =.; w=.;
				end;
			else do;
				_ass = min(max(ass, 1e-6), 1-1e-6);
				y = log(-log(_ass));
				if var_ass = 0 then _pop_raw = 0;
					else _pop_raw = _ass*(1-_ass)/var_ass;
				_pop = min(max(_pop_raw, 1), 1e6);
				_val_raw = _pop*(1-_ass);
				_val_cc  = min(max(_val_raw, 0.5), _pop - 0.5);
				_os = 1 - (_val_cc / _pop);
				_var_os = _os * (1 - _os) / _pop;
				var_y   = _var_os / ( _os * log(_os) )**2;
				log_log = rand('NORMAL', y, sqrt(var_y));
				w = min(1 / var_y, 1e6);
				end;
			run;
		proc sql noprint;
			select count(*) into :no_miss1 from spc1 where missing(ass) or ass < 0 or missing(var_ass) or var_ass < 0;
			select count(*) into :no_miss2 from spc2 where missing(w) or w <= 0 or missing(log_log);
			quit;

		data spc3;
			set spc2(keep=c p s log_log w);
			array  dummy_c {%eval(&last_c)} c1-c%eval(&last_c);
			array  dummy_p {%eval(&last_p)} p1-p%eval(&last_p);
			array  dummy_s {%eval(&last_s)} s1-s%eval(&last_s);
			do i = 1 to %eval(&last_c);
				dummy_c(i) = 0;
				if c = %eval(&ref_c) and dummy_c(i) = 0 then dummy_c(i) = -1;
				end;
			dummy_c(c) = 1;
			do i = 1 to %eval(&last_p);
				dummy_p(i) = 0;
				if p = %eval(&first_p) and dummy_p(i) = 0 then dummy_p(i) = -1;
				if p = %eval(&ref_p) and dummy_p(i) = 0 then dummy_p(i) = -1;
				end;
			dummy_p(p) = 2;
			do i = 1 to %eval(&last_s);
				dummy_s(i) = 0;
				if s = %eval(&ref_s) and dummy_s(i) = 0 then dummy_s(i) = -1;
				end;
			dummy_s(s) = 1;
			drop i c&ref_c p&first_p p&ref_p s&ref_s;
			run;

		ods noresults;
		proc genmod data=spc3;
			model log_log = &s_list2  &p_list2 &c_list2;
			weight w;
/*			output out=b1 p=sur_hat;*/
		ods output parameterestimates=parameters;
		run;

		proc iml;
			use parameters;
			read point 1 var {Estimate} into intercept;
			read point (2:%eval(&ref_s)) var {Estimate} into effect_s1;
			read point (%eval(&ref_s+1):%eval(&last_s)) var {Estimate} into effect_s2;
			read point (%eval(&last_s+1):%eval(&last_s+&ref_p-2)) var {Estimate} into effect_p1;
			read point (%eval(&last_s+&ref_p-1):%eval(&last_s+&last_p-2)) var {Estimate} into effect_p2;
			read point (%eval(&last_s+&last_p-1):%eval(&last_s+&last_p+&ref_c-3)) var {Estimate} into effect_c1;
			read point (%eval(&last_s+&last_p+&ref_c-2):%eval(&last_s+&last_p+&last_c-3)) var {Estimate} into effect_c2;

			s=j(%eval(&last_s-1),1,1);
			effect_s = effect_s1`||effect_s2`;
			sI=-(effect_s*s);
			survival=(effect_s1 // sI // effect_s2);
			p=j(%eval(&last_p-2),1,1);
			effect_p = effect_p1`||effect_p2`;
			PJ1=-(effect_p*p);
			PJ2=-(effect_p*p)/2;
			effect_p1_=effect_p1*2;
			effect_p2_=effect_p2*2;
			period=(PJ1//effect_p1_//PJ2//effect_p2_);
			c=j(%eval(&last_c-1),1,1);
			effect_c = effect_c1`||effect_c2`;
			CI=-(effect_c*c);
			cohort=(effect_c1//CI//effect_c2);

			s=(%eval(&first_s):%eval(&last_s))`;
			L_s=s-&ctr_s;
			S_s=inv(L_s`*L_s)*(L_s`*s);
			s_curve=survival-L_s*S_s;

			p=(%eval(&first_p):%eval(&last_p))`;
			L_p=p-&ctr_p;
			S_p=inv(L_p`*L_p)*(L_p`*period);
			p_curve=period-L_p*S_p;

			c=(%eval(&first_c):%eval(&last_c))`;
			L_c=c-&ctr_c;
			S_c=inv(L_c`*L_c)*(L_c`*cohort);
			c_curve=cohort-L_c*S_c;

			_year = (%eval(&first_c):%eval(&last_c+&last_s))`;
			_deviation_period = p || p_curve;
			_deviation_cohort = c || c_curve;
			_deviation_spc = intercept || S_s || S_p || S_c || s_curve`;

			create _year from _year[colname={'year'}]; append from _year;
			create _deviation_period from _deviation_period[colname={'year' 'period_deviation'}]; append from _deviation_period;
			create _deviation_cohort from _deviation_cohort[colname={'year' 'cohort_deviation'}]; append from _deviation_cohort;
			create _deviation_spc from _deviation_spc[colname={intercept S_s S_p S_c &s_list}]; append from _deviation_spc;
			run;
		quit;

		data spline_cohort; merge _year _deviation_cohort; by year;
		proc glmselect data=spline_cohort outdesign(addinputvars fullmodel)=SplineBasis noprint;
			effect spl = spline(year / details naturalcubic basis=tpf(noint)
													knotmethod=percentiles(5));
			model cohort_deviation = spl / selection=none;
			output out=spline_cohort2 predicted=spline_cohort;
			quit;
		proc transpose data=spline_cohort2 out=spline_cohort3 prefix=c; var spline_cohort; id year;run;

		data spline_period; merge _year _deviation_period; by year;
		proc glmselect data=spline_period outdesign(addinputvars fullmodel)=SplineBasis noprint;
			effect spl = spline(year / details naturalcubic basis=tpf(noint)
													knotmethod=percentiles(5));
			model period_deviation = spl / selection=none;
			output out=spline_period2 predicted=spline_period;
			quit;
		proc transpose data=spline_period2 out=spline_period3 prefix=p; var spline_period; id year; run;

		data rs_effect(rename=(i=c));
			merge _deviation_spc spline_period3 spline_cohort3;
			array arrC c&first_c-c%eval(&last_c+&last_s);
			array arrP p&first_c-p%eval(&last_c+&last_s);
			do i = %eval(&first_c) to %eval(&last_c);
				%do j = %eval(&first_s) %to %eval(&last_s);
					effect&j = -( intercept + S_s*(&j-&ctr_s) + S_p*(i-&ctr_p+&j-1) + S_c*(i-&ctr_c) +s&j + arrP(i+&j-1) + arrC(i));
				%end;
				output;
				end;
			keep i effect&first_s - effect&last_s;
			run;

		proc sort data=spc1;by c p;
		proc sort data=spc1 out=rs_pop(keep=c pop) nodupkey;by c;run;

		data rs_1yr;
			format bootstrap $5. cohort_1yr pop;
			merge rs_effect rs_pop;
			by c;
			if pop = . then delete;
			bootstrap = put(%eval(&value),z5.);
			cohort_1yr = c + %eval(&first_cohort) - 1;
			s1=exp(-exp(-effect1));
			s2=exp(-exp(-effect2));
			s3=exp(-exp(-effect3));
			s4=exp(-exp(-effect4));
			s5=exp(-exp(-effect5));
			cr = /*s1*s2*s3*s4**/s5;
			length no_miss1 no_miss2 3;
			no_miss1 = %eval(&no_miss1); no_miss2 = %eval(&no_miss2);
			format run_time e8601dt19.0; run_time = &run_time.;
			keep bootstrap cohort_1yr pop cr run_time no_miss1 no_miss2;
			run;

		data rs_5yr_1;
			format bootstrap $5. cohort_5yr;
			set rs_1yr;
			do i = %eval(&first_cohort) to %eval(&first_cohort+&last_c-&last_s);
				cohort_5yr = put(i,z4.)||'_'||put(i+%eval(&last_s)-1,z4.);
				if i <= cohort_1yr and cohort_1yr < i + %eval(&last_s) then output;
				end;
			keep bootstrap cohort_5yr pop cr run_time;
			run;
		proc sort data=rs_5yr_1;by bootstrap cohort_5yr run_time;run;
		proc means data=rs_5yr_1 noprint;
			by bootstrap cohort_5yr run_time;
			var cr;
			output out=rs_5yr(drop=_type_ _freq_) mean=cr_weight;
			weight pop;
		proc append base=rs_1yr_final data=rs_1yr force; run;
		proc append base=rs_5yr_final data=rs_5yr force; run;
		proc datasets noprint;save rs_data2 filelist &in_data rs_1yr_final rs_5yr_final;quit;

		%if %sysfunc(mod(&value,100))=0 %then %do;dm 'log; clear;';%end;
			%put bootstrap_no = &value;
	%end;
%mend bootstrap;



libname input 'D:\path_to_project\data';
libname output 'D:\path_to_project\output\spc_bootstrap_5000';
	data chk;
		set input.asrs_data;
		if Year_gap = '1';
		if surv_gp = 'RS';
		if sex in ('1' '2');
		if range <= '05';
		pri_site = strip(pri_site)||'_'||sex;
		*keep pri_site Year_intraval range tot_cases y vary;
		run;
	proc freq data=chk noprint;table pri_site / out=A1;run;
	proc sort data=chk;by pri_site Year_intraval range;
	data dup;
		set chk;
		by pri_site Year_intraval range;
		if first.range=0 | last.range=0;
		run;*0;
	proc freq data=chk noprint;table pri_site*Year_intraval / out=A1;
	proc freq data=A1 noprint;table pri_site / out=A2;run;
	data A1;set A1;if count=5 | count=4 & Year_intraval='107' | count=3 & Year_intraval='108' | count=2 & Year_intraval='109' then delete;run;*1;

	proc sql;
	create table chk2 as
		select distinct base.pri_site, base.Year_intraval, base.range, 
								coalesce(surv.tot_cases, 0) as pop, surv.ass, surv.var_ass
		from
		( select distinct * from
			( select distinct pri_site from chk ) as a,
			( select distinct Year_intraval, range from chk ) as b
		) as base
		left join chk as surv
		on base.pri_site = surv.pri_site and	
				base.Year_intraval = surv.Year_intraval and	
				base.range = surv.range
		order by base.pri_site, base.Year_intraval, base.range;
	quit;
	proc freq data=chk2 noprint;table pri_site*Year_intraval / out=A1;
	proc freq data=A1 noprint;table pri_site / out=A2;run;
	data A1;set A1;if count=5 | count=4 & Year_intraval='107' | count=3 & Year_intraval='108' | count=2 & Year_intraval='109' then delete;run;*0;



%macro ca28;
	data rs_data;
		set input.asrs_data;
		if Year_gap = '1';
		if surv_gp = 'RS';
		if sex in ('1' '2');
		if range <= '05';
		pri_site = strip(pri_site)||'_'||sex;
		keep pri_site Year_intraval range tot_cases ass var_ass;
		run;
	proc sql;
	create table rs_data2 as
		select distinct base.pri_site, base.Year_intraval, base.range, 
								coalesce(surv.tot_cases, 0) as pop, surv.ass, surv.var_ass
		from
		( select distinct * from
			( select distinct pri_site from rs_data ) as a,
			( select distinct Year_intraval, range from rs_data ) as b
		) as base
		left join rs_data as surv
		on base.pri_site = surv.pri_site and	
				base.Year_intraval = surv.Year_intraval and	
				base.range = surv.range
		order by base.pri_site, base.Year_intraval, base.range;
	quit;
	proc freq data=rs_data2 noprint;table pri_site/out=filelist;run;
	data filelist;set filelist;/*if _N_=1 | _N_=2*/;
	data filelist;set filelist;call symputx('num_files',_n_);run;

	%do n=1 %to &num_files;
/*	%do n=15 %to 15;*/
		data _null_;
			set  filelist;
			if _n_=&n;
			call symputx('pri_site',pri_site);
			run;
		%put no=&n;
		%put pri_site=&pri_site;

		data cancer;
			length cohort period pop y vary 8.;
			set rs_data2;
			if pri_site = "&pri_site";
			cohort = substr(Year_intraval,1,3)+1911;
			period = substr(Year_intraval,1,3)+1911+range;
			keep cohort period pop ass var_ass;
			run;

		%bootstrap(cancer, 5000);

		data output.rs_1yr_final_&n;format pri_site $50.;set rs_1yr_final; pri_site = "&pri_site";
		data output.rs_5yr_final_&n;format pri_site $50.;set rs_5yr_final; pri_site = "&pri_site";run;
		proc datasets noprint;delete cancer rs_1yr_final rs_5yr_final;quit;dm 'log; clear;';
	%end;

	%put &num_files;
%mend ca28;
%ca28;
