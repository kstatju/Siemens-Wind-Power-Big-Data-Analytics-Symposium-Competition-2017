proc import 
datafile="C:\Users\md351969\Desktop\Siemens\Data\Pattern_data.csv" 
dbms=csv 
out=work.Pattern_data 
replace;
run;

*Fisher's exact test was used when expected cell count was less than 5, MC=Monte carlo option is used to do the test more quickly;

%macro association(var1=,var2=);	

Proc Freq data=Pattern_data;
tables &var1 * &var2/Chisq;			
Exact Fisher/MC;				
run;

%mend  association;

%association (var1=Park_Name, var2=FactorA);
%association (var1=Park_Name, var2=FactorB);
%association (var1=Park_Name, var2=FactorC);
%association (var1=Park_Name, var2=FactorD);
%association (var1=StationID, var2=FactorA);
%association (var1=StationID, var2=FactorB);
%association (var1=StationID, var2=FactorC);
%association (var1=StationID, var2=FactorD);
%association (var1=FactorA, var2=FactorB);
%association (var1=FactorA, var2=FactorC);
%association (var1=FactorA, var2=FactorD);
%association (var1=FactorB, var2=FactorC);
%association (var1=FactorB, var2=FactorD);
%association (var1=FactorC, var2=FactorD);
%association (var1=ManualStop_during_visit, var2=code);


* Creating new variable called Day, described in the report;
data visit_data;
set visit_data;

if  0 <= visit_start_Time <= 7
 THEN Day='Night';

else if 7 < visit_start_Time <= 12
 THEN Day='Morning';

else if 12 < visit_start_Time <= 18
 THEN Day='Afternoon';

ELSE Day='Evening';
run;

*association between Parks and the new variable called DAY;
Proc Freq data=V_data2;
tables Park_Name*Day/Chisq norow nocol nopercent;
Exact Fisher/MC;
run;

Proc Freq data=V_data2;
tables Park_Name*Day*EventWarningStop/Chisq norow nocol nopercent;
Exact Fisher/MC;
run;


*Average and SD of Time difference variable ;

proc sql;
select Park_Name, avg(Time_dif) as Mean, STD(Time_dif) as Standard_Deviation
from Pattern_data
group by Park_Name;
quit;
