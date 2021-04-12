function [Tmax,Data,Results,DataN]=PlotThermo(datafile,results, xlimit, ylim14C, ylim13C,ylimCO2)

%This function takes a Dirtburner data file and a dirtburner results file
%and plots a thermograph with two panels.  Both panels have Temperature on
%the x-axis, and the top panel has two y-axes (pCO2 on left and radiocarbon
%age on right) and the bottom panel has d13C on the left (only) axis.  In 
%the top panel, pCO2 data from the data file are plotted against 
%temperature from the data file, and results are plotted as ages depicted 
%by points (with error bars) and bars that depict both the sampled 
%pyrolysis temperature interval and the age (height).  The bottom plot
%plots d13C from the results file against temperature from the data file.
%Both the data file and the results files are .txt files with no column
%labels (organization of these files is discussed below).  The function
%returns the Tmax value (temperature of maximum pCO2), the data in case
%additional plotting is desired, and the results in case additional
%calculation or plotting is required.
%
%Inputs:
%datafile:  Standard Dirtburner .txt data file containing the following
%(unlabeled) columns:
%Temperature | Date | Time | pCO2
%results:  Modified Dirtburner .txt results file containing the following
%(unlabeled) columnns:
%Interval number (integer) |  micromoles of CO2 | Upper Temperature
%Interval Limit | Fraction Modern | 1sigma Fraction Modern | d13C | 1sigma
%d13C
%xlimit: Temperature limits for plots, in form of [ll ul] where ll is lower
%limit and ul is upper limit.  Default = [0 1000]
%ylim14C:  Radiocarbon Age Limits for plots, [ll ul], default = [0 40000]
%ylim13C:  Stable Isotope Ratio Limits for plots, [ll ul], default [-30
%-20]
%ylimCO2:  pCO2 limits for left y-axis in plots, no default
%
%Outputs:
%Tmax:  scalar quantity that is the temperature at which the highest
%amounts of pCO2 were detected in the IRGA
%Data:  the loaded values of the input data file (see above for format of
%columns)
%Results:  the loaded values of the input results file (see above for
%format of columns)
%DataN:  The normalized and zeroed Data, 2 columns Temperature|pCO2

%Load Data and Results files
Data=load(datafile);
D=sortrows(Data,3);
T=D(:,3);
y=D(:,2);
R=load(results);
D0=D(:,2)-min(D(:,2));
Dn=[D(:,3) D0./sum(D0)*10000/length(D0);];

%Check results file for column with temperature intervals
if size(R(1,:))<7
    xmark=[];
    fprintf(['Results file does not seem to have temperature intervals in 3rd column',...
        'as judged by total number of columns in file (#C<8).  Empty matrix [] used instead.'])
    Fmcolumn=3;
    yFm=R(:,3); errFm=R(:,4);   
else
    xmark(1)=min(T);
    xmark(2:length(R(:,1)))=R(1:length(R(:,1))-1,3);
    xmark(length(R(:,1))+1)=max(T);
    Fmcolumn=4;
    yFm=R(:,4);errFm=R(:,5);
end


yages=yFm;  %This is fraction modern, as opposed to an age - formerly -8033ln(yFm)


for ii=1:length(errFm)
    err(ii)=(2*errFm(ii));
end

%Set defaults
if nargin<6
    ylimCO2=[0 roundn(1.1*max(D(:,2)),1)];
end
if nargin<5
    ylim13C=[-30 -20];ylimCO2=[0 roundn(1.1*max(D(:,2)),1)];
end
if nargin<4
    ylim14C=[0 1];ylim13C=[-30 -20];ylimCO2=[0 roundn(1.1*max(D(:,2)),1)];
end
if nargin<3
    xlimit=[0 1000];ylim14C=[0 40000];ylim13C=[-30 -20];ylimCO2=[0 roundn(1.1*max(D(:,2)),1)];
end
    

f=figure;
%First plot the thermograph with ages
s1=subplot(10,1,1:7);
ppos=(get(s1,'Position'));
sqfc=ppos(3)/ppos(4);
title([datafile,' Thermograph'])
plot(T,y,'k','LineWidth',2.5);
[C I]=max(y);
Tmax=(T(I,1));
text(1.1*Tmax,C,['T_{max} = ',num2str(roundn(Tmax,0)),'^{o}C']);
ylabel('pCO_{2} (\mumol/mol)','FontSize',14);
h1=gca;
set(gca,'YLim',ylimCO2,'XLim',xlimit);
h2=axes('Position',get(h1,'Position'),'YLim',ylim14C);
set(h2,'YAxisLocation','right','Color','none','XTickLabel',[]);
set(h2, 'XLim', xlimit,'Layer','Top','YLim',ylim14C);
%Loop to prevent negative ages from being passed into unevenbarerrcirc
%function:
for k=1:length(yages)
    if yages(k)<=0
        base_D14C(k)=yages(k);
        top_D14C(k) = 0;
    else
        base_D14C(k)=0;
        top_D14C(k)=yages(k);
    end
end
unevenbarerrcirc(xmark,top_D14C,base_D14C,err,sqfc,[0 0 0]);
ylabel('Fraction modern','FontSize',14);

%Calculate and plot weighted bulk Fm from dirtburner results, if all
%results present
if sum(isnan(R(:,Fmcolumn)))==0     %If there are not any NaN's
    MolSum=sum(R(:,2));
    MolFrac=R(:,2)./MolSum;
    for n=1:length(yFm)
        D14Cweight(n)=(yFm(n)-1)*1000*MolFrac(n);
        Fmweighted(n)=yFm(n)*MolFrac(n);
        E2(n)=errFm(n)^2;
        d13Cweight(n)=R(n,Fmcolumn+2)*MolFrac(n);
    end
    WBD14C=sum(D14Cweight);
    WBd13C=sum(d13Cweight);
    WBFm=WBD14C/1000+1;
    WBFmsigma=sqrt(sum(E2));
    WBAgeerr=(-8033*log(WBFm-2*WBFmsigma)-(-8033*log(WBFm+2*WBFmsigma)))/2;
    WBd13Cerr=sqrt(sum(R(:,Fmcolumn+3).^2));
    line([min(xlimit) max(xlimit)],[(WBFm) (WBFm)],'Color','r')
    line([min(xlimit) max(xlimit)],[(WBFm+2*WBFmsigma) (WBFm+2*WBFmsigma)],'Color','r','LineStyle','--')
    line([min(xlimit) max(xlimit)],[(WBFm-2*WBFmsigma) (WBFm-2*WBFmsigma)],'Color','r','LineStyle','--')
    text(min(T),1.05*(WBFm),['Weig. Bulk F_{m} = ',...
        num2str(roundn((WBFm),-3)),' +/- ',...
        num2str(roundn(WBFmsigma,-4))],'Color','r');
else
    fprintf('Incomplete results file with missing Fm.  Calculated bulk age neither performed nor plotted.')
end



%Plot d13C in remaining panel
s2=subplot(10,1,9:10);
mp=zeros(length(xmark)-1);
for h=2:(length(xmark))
    mp(h)=(xmark(h-1)+xmark(h))/2;  %find rectangle midpoints
end
errorbar(nonzeros(mp),R(:,Fmcolumn+2),R(:,Fmcolumn+3),'ko')
%Check for complete d13C dataset and plot calculated bulk if complete.
if sum(isnan(R(:,Fmcolumn+3)))==0     %If there are not any NaN's
    MolSum=sum(R(:,2));
    MolFrac=R(:,2)./MolSum;
    for n=1:length(R(:,Fmcolumn+3))
        d13Cweight(n)=R(n,Fmcolumn+2)*MolFrac(n);
    end
    line([min(xlimit) max(xlimit)],[WBd13C WBd13C],'Color','b')
    line([min(xlimit) max(xlimit)],[WBd13C-2*WBd13Cerr WBd13C-2*WBd13Cerr],'Color','b','LineStyle','--')
    line([min(xlimit) max(xlimit)],[WBd13C+2*WBd13Cerr WBd13C+2*WBd13Cerr],'Color','b','LineStyle','--')
    text(300,-27,['Weig. Bulk \delta^{13}C = ', num2str(roundn(WBd13C,-1)),' +/- ',num2str(roundn(2*WBd13Cerr,-2))],'Color','b')
else
    fprintf('Incomplete results file with missing d13C.  Calculated bulk d13C neither performed nor plotted.')
end
xlabel('Temperature (^{o}C)','FontSize',14);
h3=gca;
set(h3,'XLim',xlimit,'YLim',ylim13C);

%Calculate outputs

Data=D0;
DataN=Dn;
Results=R;