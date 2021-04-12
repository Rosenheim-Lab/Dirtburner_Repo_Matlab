function [D,Dn,Tmax, f,CompAge, Compd13C]=DecompDB3(filename,Int)
%************Decompose Dirt Burner Data************************************
%
%FOR USE ON RUNS FROM 25MAY2012 ONWARD!!!
%
%DecompDB(filename,Int) plots a thermograph from data in text file
%'filename' (6 columns: date number, pCO2,temperature,He, He-side, O2-side)
%using temperature intervals in Int.
%
%Outputs
%
%[D]=DecompDB(filename,Int) plots a thermograph from data in teh text file
%'filename' (6 columns: timestamp, pCO2(umol/mol), temperature, 
%He main flow (mL/min), He sidearm flow (mL/min), O2 sidearm flow(mL/min),
%using the temperature intervals in Int (use [] for no temperature
%intervals when a shape run is performed).  D contains the data in the
%file 'filename'.
%
%[D,Dn]=DecompDB(filename,Int) plots a thermograph from data in teh text file
%'filename' (6 columns: timestamp, pCO2(umol/mol), temperature, 
%He main flow (mL/min), He sidearm flow (mL/min), O2 sidearm flow(mL/min),
%using the temperature intervals in Int (use [] for no temperature
%intervals when a shape run is performed).  D contains the data in the
%file 'filename'.  Dn contains zeroed and then normalized data in two
%columns, temperature and zeroed, normalized pCO2.
%
%[D,Dn,Tmax]=DecompDB(filename,Int) plots a thermograph from data in teh text file
%'filename' (6 columns: timestamp, pCO2(umol/mol), temperature, 
%He main flow (mL/min), He sidearm flow (mL/min), O2 sidearm flow(mL/min),
%using the temperature intervals in Int (use [] for no temperature
%intervals when a shape run is performed).  D contains the data in the
%file 'filename'.  Dn contains zeroed and then normalized data in two
%columns, temperature and zeroed, normalized pCO2.  Tmax is the
%temperature at which maximum CO2 is observed.
%
%***THE FOLLOWING OUTPUTS ARE PLANNED IN FUTURE VERSIONS AND HAVE NOT BEEN
%ACTIVATED!****************************************************************
%
%[D,Dn,Tmax,f]=DecompDB(filename,Int) plots a thermograph from data in teh text file
%'filename' (6 columns: timestamp, pCO2(umol/mol), temperature, 
%He main flow (mL/min), He sidearm flow (mL/min), O2 sidearm flow(mL/min),
%using the temperature intervals in Int (use [] for no temperature
%intervals when a shape run is performed).  D contains the data in the
%file 'filename'.  Dn contains zeroed and then normalized data in two
%columns, temperature and zeroed, normalized pCO2. Tmax is the
%temperature at which maximum CO2 is observed.  Also returns f - 
%fraction matrix explaining how much OC was in each Gaussian component.
%
%[D,Dn,Tmax,f, CompAge]=DecompDB(filename, Int) plots a thermograph from data in the 
%text file 'filename' (6 columns: timestamp, pCO2(umol/mol), 
%temperature, He main flow (mL/min), He sidearm flow (mL/min), O2 sidearm 
%flow(mL/min), using the temperature intervals in Int (use [] for no 
%temperature intervals when a shape run is performed).  D contains the data
%in the file 'filename'.  Dn contains zeroed and then normalized data in 
%two columns, temperature and zeroed, normalized pCO2. Tmax is the
%temperature at which maximum CO2 is observed.  Also returns f - 
%fraction matrix explaining how much OC was in each Gaussian component and
%CompAge - the age of each component.  If commanded to return CompAge, you 
%will be asked for a .txt file containing the results in 6 columns:  
%Interval number, umoles of CO2, Fm, sigma Fm, d13C, and sigma d13C.
%
%[D,Dn,Tmax,f, CompAge Compd13C]=DecompDB(filename, Int) plots a thermograph from 
%data in the text file 'filename' (6 columns: timestamp, pCO2
%(umol/mol), temperature, He main flow (mL/min), He sidearm flow (mL/min), 
%O2 sidearm flow(mL/min), using the temperature intervals in Int (use [] 
%for no temperature intervals when a shape run is performed).  D contains 
%the data in the file 'filename'.  Dn contains zeroed and then normalized 
%data in two columns, temperature and zeroed, normalized pCO2. Tmax is the
%temperature at which maximum CO2 is observed.  Also returns
%f - fraction matrix explaining how much OC was in each Gaussian component 
%and CompAge - the age of each component.  If commanded to return CompAge, 
%you will be asked for a .txt file containing the results in 6 columns:  
%Interval number, umoles of CO2, Fm, sigma Fm, d13C, and sigma d13C.
%Compd13C - the d13C composition of each component.  
%
%All input files 'filename' should be fixed in excel to remove blank
%entries in the first row - THIS SHOULD NOT A PROBLEM AFTER MAY 24, 2012.  
%
%Brad E. Rosenheim, Tulane University, November, 2011
%modified May, 2012, Brad E. Rosenheim

D=load(filename);
%D = dlmread(filename, '\t');
D0=D(:,2)-min(D(:,2));
Dn=[D(:,3) D0./sum(D0)*10000/length(D0);];  %normalized to measurement
    %length of 10000 points for comparison to plots of data taken at
    %different rates

%Plot figure 1 - zeroed thermograph
f1=figure;
set(f1,'Name',['Thermograph for ',filename,', zeroed'])
plot(D(:,3),D0,'k','LineWidth',2) %plot zeroed pCO2
xlabel('Temperature, ^{o}C')
ylabel('pCO_{2} (\mumol^{.}mol^{-1})')
%Draw Temperature interval lines
for n=1:length(Int)
    line([Int(n) Int(n)],[min(D0) max(D0)],'Color','r','LineStyle','--');
end

[C I]=max(D(:,2));
Tmax=(D(I,3));
text(1.1*Tmax,C,['T_{max} = ',num2str(roundn(Tmax,0)),'^{o}C']);

%Plot Figure 2, normalized and zeroed thermograph
f2=figure;
set(f2,'Name',['Normalized Thermograph for ',filename])
plot(Dn(:,1),Dn(:,2),'b','LineWidth',2) %plot normalized pCO2
xlabel('Temperature, ^{o}C')
ylabel('normalized pCO_{2} (mol^{-1})')
%Draw Temperature interval lines
for n=1:length(Int)
    line([Int(n) Int(n)],[min(Dn(:,2)) max(Dn(:,2))],'Color','k','LineStyle','--');
end
text(1.1*Tmax,0.9*max(ylim),['T_{max} = ',num2str(roundn(Tmax,0)),'^{o}C']);

