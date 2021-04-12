function [ FmBC ] = BlankCorrect14C( Fmdat,m,d,sigm,sigd )
%BlankCorrectC14 corrects fractions modern by a mass of dead and modern
%carbon as per Santos et al., 2007 NIM-B.  
%   This function accepts single-value Fm's or vectors of Fm, as well as
%   the masses of the live (m) and dead (d) blanks and their uncertainties.
%    These can be calculated after Santos et al., 2010 (Radiocarbon).
%    Defaults are no blanks and zero uncertainties - i.e. no correction.
%    Blank masses and uncertainties must be entered.  If only masses are
%    entered, the default uncertainty on the mass is 50% the mass of each
%    blank as per Santos et al., 2007 NIM-B.  
%Syntax
%   [FmBC] = BlankCorrect14C(Fmdat,m,d,sigm,sigd)
%       Fmdat is a row vector for each sample of fraction modern | analytical uncertainty |
%       sample size (micromoles C) | number of runs.  It can also be a matrix formed from
%       several of these vectors (different sample, different row).  FmBC will be a vector or a matrix 
%       consisting of corrected Fm and corrected uncertainties.  m, d, 
%       sigm, and sigd are all properties of the blank contamination.
%         At the writing of this function, these properties are determined 
%       in Fernandez et al., 2014, Analytical Chemistry.  These are
%       functional inputs in order to be able to use different amounts of
%       blank contamination with changes in protocol and equipment (such as 
%       during the move of the equipment to USF).
%       
%   [FmBC] = BlankCorrect14C(Fmdat)
%       Default blank values are 0 micrograms and 0 error - in other words
%       the values will simply be the input values.
%
%   [FmBC] = BlankCorrect14C(Fmdat,m,d)
%       Default blank mass uncertainties if blank masses are entered are
%       50% of the blank mass as per Santos et al., 2007 NIM-B
%       
%
%Written by B.E. Rosenheim, USF-CMS, 18Nov2014
%Revised by B.E. Rosenheim, USF-CMS, 18Mar2015
%   -added functionality for age limits based on Michelle Guitard's SMO
%   data
%   -fixed bug with error propagation calculation
%   -allowed for variable blank corrections entered as matrices

%% Defaults
if nargin==4
    sigd=0.5*d;
end
if nargin==3
    sigm=0.5*m;
    sigd=0.5*d;
end
if nargin==2
    d=0;
    sigm=0.5*m;
    sigd=0.5*d;
end
if nargin==1
    m=0;
    d=0;
    sigm=0.5*m;
    sigd=0.5*d;
end

%% Calculation of blank-corrected Fm and associated uncertainty

runs = Fmdat(:,4);
m = m*runs;
sigm = sigm*runs;
d = d*runs;
sigd = sigd*runs;

mass=Fmdat(:,3).*12;
if length(m)~=length(d) || length(m)~=length(sigm) || length(m)~=length(sigd) || length(d)~=length(sigm) || length(d)~=length(sigd) ||length(sigm)~=length(sigd)
    error('Please ensure that all blank mass and uncertainty inputs are of equal size! Restart function.')
end
%Turn scalar inputs for blanks into vectors for if loop below.  
if length(m)==1
    m=ones(length(Fmdat(:,1)),1)*m;
    d=ones(length(Fmdat(:,1)),1)*d;
    sigm=ones(length(Fmdat(:,1)),1)*sigm;
    sigd=ones(length(Fmdat(:,1)),1)*sigd;
end

FmBC=zeros(length(Fmdat(:,1)),2);
for j=1:length(Fmdat(:,1))
    FmBC(j,1)=(Fmdat(j,1)-m(j)/mass(j)+Fmdat(j,1)*d(j)/mass(j)+Fmdat(j,1)*m(j)/mass(j));
    FmBC(j,2)=sqrt(Fmdat(j,2)^2*((d(j)+m(j)+mass(j))/mass(j))^2+sigd(j)^2*((Fmdat(j,1))/...
        mass(j))^2+sigm(j)^2*((Fmdat(j,1)-1)/mass(j))^2);
end



end

