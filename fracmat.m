function [BinSumCO2,f]=fracmat(T,nCO2,Model,TsumCO2)
%
%function [BinSumCO2,f]=fracmat(nCO2,Model,TSumCO2)
%
%INPUT VARIABLES
%
%nCO2 is an [nx1]row vector of temperatures at which samples were binned.
%
%Model is an [mxn] matrix containing the component values of pCO2 at each
%temperature increment.
%
%TSumCO2 is an [mx1] matrix containing the sum, at each temperature
%increment, of total CO2 at that increment.  (This is the value closest to
%the actual CO2 measurement taken during run time).
%
%mark is the marker matrix from the integration loop of the main frame
%
%OUTPUT VARIABLES
%
%BinSumCO2 is the sum of total CO2 from each temperature bin
%
%f is the coeffient matrix containing the ratio of each component's
%proportional contribution to BinSumCO2
%
%Brad E. Rosenheim
%Jan2007
%last modified 18Jan2007

%Construct marker matrix
mark=ones(length(nCO2)+1,1);
mark=mark(:);
for z=1:length(nCO2)
    mark(z+1)=max(find(T<nCO2(z)));
end
%mark(z+1)=length(T);

%For every measurement bin,
for m=1:length(nCO2)
    if m==1
        %Calculate a denominator for the component fractions.
        BinSumCO2(m)=sum(TsumCO2(mark(m):mark(m+1)));
        %For every component within every bin,
        for n=1:length(Model(1,:))
            %compute the fraction of the binCO2 that the component
            %constitutes.
            f(m,n)=(sum(Model(mark(m):mark(m+1),n)))/BinSumCO2(m);
        end
    else
        %Calculate a denominator for the component fractions (after first bin).
        BinSumCO2(m)=sum(TsumCO2((mark(m)+1):(mark(m+1))));
        %For every component within every bin (after first bin),
        for n=1:length(Model(1,:))
            %compute the fraction of the binCO2 that the component
            %constitutes (after first bin).
            f(m,n)=(sum(Model((mark(m)+1):mark(m+1),n)))/BinSumCO2(m);
        end
    end
end
for w=1:length(f(:,1))
    Check(w)=sum(f(w,:));
    if roundn(Check(w),-2)~=1
        Check(w)=1;
    else
        Check(w)=0;
    end
end
if sum(Check)>0
    warning('Fractions are inaccurate.  Some temp. increments do not sum to 1.')
end