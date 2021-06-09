function [ ages_BC_rounded ] = round_14C_ages( ages_BC )
%round_14C_ages rounds ages calculated from fraction modern radiocarbon after 
%Stuiver and Pollach, 1977 .  
%   This function works on calculated ages, and it can work with the output
%   BlankCorrect14C with slight modification.
%
%	Note that this function will work on any set of positive radiocarbon ages
%	no matter whether they are blank corrected or not. However, this function
%	should not be used prior to calculated blank corrected fractions modern due
%	to the possibility for rounding errors. 
%      
%Syntax
%   [ages_BC_rounded] = round_14C_ages(ages_BC)
%       ages_BC can be any positive ages calculated from fractions modern 14C. Data
%	should be in two columns:
%		Ages | Age Uncertainties 
%       
%   [ages_BC_rounded] = round_14C_ages([-8033*ln(FmBC(:,1); 8033*FmBC(:,2)./FmBC(:,1)])
%       With this syntax on the right side of the equation, this function can be used 
%	directly with the output of BlankCorrect14C. This is the best time to round 
%	ages. 
%
%   
%       
%
%Written by B.E. Rosenheim, USF-CMS, 9Jun2021


%% Initialize matrix to be filled in two loops
rounded_ages = ones(size(ages_BC(1)), 2)

%% Run a for loop to fill rounded_ages
for ind=1:length(ages_BC(:,1))
	if ages_BC(ind, 1)<1000
            rounded_ages(ind, 1) = roundn(ages_BC(ind,1)/5, 0)*5
        elseif ages_BC<10000
            rounded_ages(ind, 1) = roundn(ages_BC(ind,1)/10, 0)*10
        elseif age<20001
            rounded_ages(ind, 1) = roundn(ages_BC(ind,1)/50, 0)*50
        else
            rounded_ages(ind, 1) = roundn(ages_BC(ind,1)/100, 0)*100
	end
end

for ind=a:length(ages_BC(:,1))
        if ages_BC(ind,2)<100
            rounded_ages(ind,2) = roundn(ages_BC(ind,2)/5, 0)*5
        elseif ages_BC(ind,2)<1001
            rounded_ages(ind,2) = roundn(ages_BC(ind,2)/10, 0)*10
        else
            rounded_ages(ind,2) = roundn(ages_BC(ind,2)/100, 0)*100
	end
end

