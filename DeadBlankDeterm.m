function [dC, sigdC, linvals_d]=DeadBlankDeterm(datafile,sigmas,analystindexfile)
%
%This function calculates and graphically displays the amount of blank
%likely from a series of radiocarbon determinations of different masses of
%a substance of known radiocarbon content.  The input file is a .txt file
%that consists of:
%Sample Mass (micrograms C) | Fraction modern | Fraction modern
%uncertainty|Date Analyzed|Date Reported|Analyst|Run number|Split
%The other input concerns the true value of the measured entity in fraction
%modern (1.29 for C3, 1.0398 for Ox-I and, 0 for graphite). 
%Ox-I - 1.0398
%
%Written: B.E. Rosenheim, USF CMS
%Revised: B.E. Rosenheim, USF CMS 12Oct2016
%   -plot blank data showing the date of analysis, date of 14C
%   determination and the analyst as color-coded points

%Load data
D=load(datafile);

%Plot the blank determinations on log-log space
figure(201)
title('Dead Blank Analysis, Whole')
loglog(D(:,1)/1000,D(:,2),'ok','MarkerFaceColor','g');

%Now plot contours of equal blank contamination along different masses of
%blank.
mM=ones(length(D(:,1)),1);  %Vector of isotope dilution estimations of blank
for n=1:length(D(:,1))
    mM(n)=D(n,1)-(D(n,2)*D(n,1))/1.0398; %Equation 8 of Alvarez et al,
    %supplemental materials
end
%It is necessary to avoid negative masses here - we assign them a mass of
%0.1 microgram.
ctr=1;
for n=1:length(mM)
    if mM(n)>=0
        mMnonZ(ctr)=mM(n);
        ctr=ctr+1;
    end
end
mC=mean(mMnonZ);
dC=mC;
sigmC=std(mMnonZ);
sigdC=sigmC;
%Find the order (power) of the x-axis limits in order to plot correctly on
%a logarithmic scale.  If the difference in order is too small, increase
%the order by one on the high end and decrease by one on the low end to
%give at least three orders of magnitude.
Xes=order(xlim);
if Xes(1)==Xes(2);
    Xes(1)=Xes(1)-1;
    Xes(2)=Xes(2)+1;
end
xint=logspace(Xes(1), Xes(2),50);
minline=floor(min(mM-sigmC*sigmas));
if minline<1
    minline=1;
end
maxline=ceil(max(mM+sigmC*sigmas));
rightline=maxline-minline+2; %Add one for masses, and one for the one that is skipped by subtraction
lines=ones(length(xint),rightline); 
negs=ylim;
negs=[0.95*negs(1) 1.05*negs(2)];
for m=minline:maxline
    linvals_d(m)=m;
    for n=1:length(xint)
        
        lines(n,m+1)=((xint(n)-m/1000)*1.0398)/xint(n);
        if lines(n,m+1)<negs(1)
            lines(n,m+1)=NaN;
        end
        
    end
end
linvals_d=linvals_d(1:m-1);
lines(:,1)=xint(:);
%Construct matrix for right answer and error bounds
corrlines=ones(length(xint),4);
for z=1:length(xint)
    corrlines(z,3)=((xint(z)-mC/1000)*1.0398)/xint(z);
    if corrlines(z,3)<negs(1)
            corrlines(z,3)=NaN;
    end
    corrlines(z,2)=((xint(z)-(mC-sigmC*sigmas)/1000)*1.0398)/xint(z);
    if corrlines(z,2)<negs(1)
            corrlines(z,2)=NaN;
    end
    corrlines(z,4)=((xint(z)-(mC+sigmC*sigmas)/1000)*1.0398)/xint(z);
    if corrlines(z,4)<negs(1)
            corrlines(z,4)=NaN;
    end
end
corrlines(:,1)=xint(:);
%show plots
hold on
for k=1:length(lines(1,:))-2
    loglog(lines(:,1),lines(:,k+1),'k--')
    %text(xint(1)*5,(minline-(k-1))/xint(1)*5,[num2str(minline+(k-1)),' /mug']);
    %Text line doesn't throw an error, but it doesn't show up.  I need to
    %think it through better.
end
%plot the "right" answer and error bounds
loglog(corrlines(:,1),corrlines(:,2),'-.r')
loglog(corrlines(:,1),corrlines(:,4),'-.r')
loglog(corrlines(:,1),corrlines(:,3),'-r','LineWidth',2)

xlabel('Sample Size (mg C)')
ylabel('Fraction Modern')

%Plot by date analyzed 
figure(202)
title('Dead Blank by Analysis Date');
scatter(D(:,1)/1000,D(:,2),30,D(:,5),'filled','MarkerEdgeColor','k')
set(gca,'XScale','log','YScale','log');
c1=colorbar;
yy=get(c1,'ylim');
cylabs=linspace(min(yy),max(yy),5);
mdatenum=cylabs+693960; %Convert to matlab date numbers assuming 1900 start date from Excel (macs use 1904)
labdates=datestr(mdatenum,'yyyy-mmm');
set(c1,'Ticks',cylabs,'TickLabels',labdates)
hold on
for k=1:length(lines(1,:))-2
    loglog(lines(:,1),lines(:,k+1),'k--')
    %text(xint(1)*5,(minline-(k-1))/xint(1)*5,[num2str(minline+(k-1)),' /mug']);
    %Text line doesn't throw an error, but it doesn't show up.  I need to
    %think it through better.
end
%plot the "right" answer and error bounds
loglog(corrlines(:,1),corrlines(:,2),'-.r')
loglog(corrlines(:,1),corrlines(:,4),'-.r')
loglog(corrlines(:,1),corrlines(:,3),'-r','LineWidth',2)

xlabel('Sample Size (mg C)')
ylabel('Fraction Modern')

%Plot by analyst
figure(203)
title('Modern Blank by Analyst');

set(gca,'XScale','log','YScale','log');

namecode=D(:,6);    %Use analyst numbers as indices to find abbreviations
%Renumber analysts based on input file (data analyst file has all of them).
%Keep abbreviations but assign new indices.
fileID=fopen(analystindexfile);
C=textscan(fileID,'%s');
fclose(fileID);
abbcounter=1;   %abbreviation counter
cctemp={[0.2123 0.2138 0.6270] [0.0165 0.4266 0.8786] 'k' [0.1453 0.7098 0.6646] [0.5709 0.7485 0.4494] [0.9139 0.7258 0.3063] [0.9763 0.9831 0.0538]};
sytemp={'v' 'o' 'p' 's' '^' '<' 'h' '>'}; %Note: ensure that symbols template and color code template are different lengths
for n=1:length(D(:,1))
    if n ==1
        k1=rem(abbcounter,length(sytemp));
        k2=rem(abbcounter,length(cctemp));
        legtext(abbcounter)=C{1}(namecode(n));
        anaind(n)=abbcounter;
        if k1==0
            symbols{abbcounter}=sytemp{end};
        else
            symbols{abbcounter}=sytemp{k1};
        end
        if k2==0
            colorcodes{abbcounter}=cctemp{end};
        else
            colorcodes{abbcounter}=cctemp{k2};
        end
    else
        if namecode(n)~=namecode(n-1)
            abbcounter=abbcounter+1;
            legtext(abbcounter)=C{1}(namecode(n));
        end
        anaind(n)=abbcounter;
        k1=rem(abbcounter,length(sytemp));
        k2=rem(abbcounter,length(cctemp));
        if k1==0
            symbols{abbcounter}=sytemp{end};
        else
            symbols{abbcounter}=sytemp{k1};
        end
        if k2==0
            colorcodes{abbcounter}=cctemp{end};
        else
            colorcodes{abbcounter}=cctemp{k2};
        end
    end
end

hold on
for j=1:abbcounter
    k3=find(anaind==j);
    hdl=plot(D(k3,1)/1000,D(k3,2));
    set(hdl,'MarkerFaceColor',colorcodes{j},'Marker',symbols{j},'MarkerEdgeColor','k','MarkerSize', 10,'LineStyle','none');
end


%PlotbyNameCode([D(:,1)/1000 D(:,2)],namecode,ss)


legend(legtext)


for k=1:length(lines(1,:))-2
    loglog(lines(:,1),lines(:,k+1),'k--')
    %text(xint(1)*5,(minline-(k-1))/xint(1)*5,[num2str(minline+(k-1)),' /mug']);
    %Text line doesn't throw an error, but it doesn't show up.  I need to
    %think it through better.
end
%plot the "right" answer and error bounds
loglog(corrlines(:,1),corrlines(:,2),'-.r')
loglog(corrlines(:,1),corrlines(:,4),'-.r')
loglog(corrlines(:,1),corrlines(:,3),'-r','LineWidth',2)

xlabel('Sample Size (mg C)')
ylabel('Fraction Modern')