function [ D_BC ] = D14C_BC( BlankCorrect14C_output,sigmas )
%Using the output of BlankCorrect14C, this routine calculated ages and age
%limits.
%
%   This function accepts the output of BlankCorrect14C (column of blank-
%   corrected Fm's and a column of sigma values on these.  It rounds
%   according to conventions established by Stuiver and Pollach, 1977.  It
%   defaults to 2-sigma ages, but this can be changed using the last input
%   variable.  For age limits, the routine keys on 1. negative Fm's from
%   the blank correction routine and 2. values of Fm-(sigmas)*sigFm that
%   are less than 0.  These ages incur 2-sigma age limits, whereby the
%   Fm+2sigma (or 2 sigma in the case of a negative fraction modern) are 
%   output as ages.
%
%Syntax:
%   [D_BC]=D14C_BC(BlankCorrect14C_output, sigmas)
%       BlankCorrect14C is a 2-column matrix, with one column being
%       blank-corrected fractions modern and the second column being sigma
%       values for those fractions modern, with error propagated into them
%       from blank correction.  Sigmas is a scalar representing to what
%       precision you want to report ages; default is 2 sigma.  Output is a
%       two column matrix of the same length as the input matrix with age
%       in column one and error in column two.
%
%Written by B.E. Rosenheim, USF-CMS, 18Mar2015


%Set defaults
if nargin<2
    sigmas=2; %2-sigma default
end

D=BlankCorrect14C_output;
%Calculate ages and errors
for j=1:length(D(:,1))
    if D(j,1)<=0
        D_BC(j,1)=(2*D(j,2)-1)*1000;
        D_BC(j,2)=-999;
    elseif D(j,1)<2*D(j,2)&&D(j,1)>0
        D_BC(j,1)=((D(j,1)+2*D(j,2))-1)*1000;
        D_BC(j,2)=-999;
    else
        D_BC(j,1)=(D(j,1)-1)*1000;
        D_BC(j,2)=D(j,2)*1000;
    end
    
end

        
        
