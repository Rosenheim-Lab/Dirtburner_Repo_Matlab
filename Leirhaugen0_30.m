%Svalbard Leirhaugen 0-30cm
%
%This command shows the gaussian decomposition of sediment from a core from
%Svalbard, originally treated by Gesine Mollenhauer. The first figure  
%shows the stacked thermographs (bold red) and gaussian components solved 
%in iteratively using nlleasq routine.  The second figure shows individual 
%components.  A third part of the exercise involves calculating a linear
%mixing model in D14C/inverse cumulative yield space.  This is used as a 
%further constraint on the young component age.
%
%Inputs:  
%   Leirhaugen0-30.txt - 2 column matrix, |Temperature|pCO2(zeroed|
%   Leirhaugen0-30_Results.txt - 6-column matrix, |Interval Number |
%   micromoles CO2 | Fm | sigFm | d13C | sigd13C
%
%The data were analyzed by Brad Rosenheim, April 2011.
%
%B.E. Rosenheim
%26Apr2011
%Last Modified 26Apr2011



%load the data files for data and for results.
load Leirhaugen0_30_AT.txt;
D=sortrows(Leirhaugen0_30_AT,1);
T=D(:,1);
y=D(:,2);

load Leirhaugen0_30_DeminResults.txt;
R=Leirhaugen0_30_DeminResults;  %These are currently dummy results (26Apr2011)


%Create data array with pin matrices containing initial guess information
%for data set in D
pin=[500 320 50 400 400 40 900 450 60 400 505 35 150 605 40];


%Create vectors of temperatures at which valve was switched for each result.
nCO2=[307 382 440 max(D(:,1))];


%Take an initial guess from PIN for the Gaussian components that would 
%best approximate the thermograph data in D
ModelIg= pin(1)*exp(-((T-pin(2))/pin(3)).^2)+pin(4)*exp(-((T-pin(5))/pin(6)).^2)+pin(7)*exp(-((T-pin(8))/pin(9)).^2)+pin(10)*exp(-((T-pin(11))/pin(12)).^2);
[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]=nlleasqr(T,y,pin,'Gauss5',0.00001,100);
%Calculate the model reproduction of the curve and the root mean square
%of the residual between model and data
Modelout= p(1)*exp(-((T-p(2))/p(3)).^2)+p(4)*exp(-((T-p(5))/p(6)).^2)+p(7)*exp(-((T-p(8))/p(9)).^2)+p(10)*exp(-((T-p(11))/p(12)).^2)+p(13)*exp(-((T-p(14))/p(15)).^2);
RMS=sqrt(sum((Modelout-y).^2)/length(y));
%Compute individual components of the model
for h=1:length(pin)/3
    Models(:,h)=p(3*h-2)*exp(-((T-p(3*h-1))/p(3*h)).^2);
end

%One Gaussian is a baseline Gaussian, and is added into the other 4:
Models(1:290,1)=Models(1:290,1)+Models(1:290,5);
Models(291:380,1)=Models(291:380,1)+0.5*Models(291:380,5);
Models(291:380,2)=Models(291:380,2)+0.5*Models(291:380,5);
Models(381:450,1)=Models(381:450,1)+0.05*Models(381:450,5);
Models(381:450,2)=Models(381:450,2)+0.75*Models(381:450,5);
Models(381:450,3)=Models(381:450,3)+0.02*Models(381:450,5);
Models(381:450,4)=Models(381:450,4)+0.18*Models(381:450,5);
Models(451:510,2)=Models(451:510,2)+0.25*Models(451:510,5);
Models(451:510,3)=Models(451:510,3)+0.35*Models(451:510,5);
Models(451:510,4)=Models(451:510,4)+0.45*Models(451:510,5);
Models(510:end,4)=Models(510:end,4)+Models(510:end,5);
MM=Models(:,1:4);
clear Models
Models=MM;

       
%Create figure to plot thermograph with model fit and components
figure(603)
subplot(2,1,1)
%Thermograph and model fit
plot(T,y,'r','LineWidth',2.5)
hold on
plot(T,ModelIg,'c')
plot(T,Modelout,'k')
ylabel('pCO_{2} (\mumol/mol)','FontSize',14)
legend('Leirhaugen 0-30cm Thermograph','Initial Guess','5 Gaussian Model')
subplot(2,1,2)
%Thermograph with components
plot(T,y,'r','LineWidth',2.5)
hold on
for kk=1:length(Models(1,:))
    plot(T,Models(:,kk),'Color',[length(Models(1,:))/(kk-1+length(Models(1,:))) 0.5 kk/length(Models(1,:))])
end

xlabel('Temperature (^{o}C)','FontSize',14)
ylabel('pCO_{2} (\mumol/mol)','FontSize',14)



%Now calculate component ages once the thermographs are satisfactorily
%fit with funcformTC18
%Create vectors of measurements (d13C and Fm)
A=R(:,3); SI=R(:,5);
    %Calculate Uncorrected Fractions modern:
    UA=A./(1-2*(25+SI)/1000);
    T=D(:,1);
    y=D(:,2);
    for k=1:length(T)
        %Calculate T Sum CO2 (ppm) at each temperature, create column 
        %vector of summed CO2 (ppm) at each T
        TsumCO2(k)=sum(Models(k,:));
    end
    %Calculate fractions of each component, as modeled, in each measurement.
    [BinSumCO2,f]=fracmat(T,nCO2,Models,TsumCO2);
    B=BinSumCO2;
   
    
    %Use algebraic solution of equation set for at least an initial guess,
    %if the system is overdetermined
   
    xx=(f'*f)^-1*(f'*UA(:)); %Calculate uncorrected component Fm's
    yy=(f'*f)^-1*(f'*SI(:));    %Calculate component d13C's
    xxCorr=xx.*(1-2*(25+yy)/1000);   %Correct xx for d13C
    %Here, I need to install a check for rank so that I don't get stuck
    %with an error message in program.
    CompAge=-8033*log(xxCorr);
    Compd13C=yy;
    Fm=exp(-CompAge./8033);
    D14C=(Fm-1)*1000;
    
    %Use algebraic solution to check constrained minimization routine;
    %If system is overdetermined, use constrained minimizing nonlinear
    %regression. First, create an anonymous function to express the quantity
    %that needs to be minimized.  This will be the residual between the
    %measured ages and those calculated from the regression solution using the
    %fractions in f.
    
    %igAge=ones(size(xx));
    %igSIR=-27*ones(size(yy));
    %lbAge=zeros(size(CompAge));
    %a=igAge;
    %b=igSIR;
    %AgeMin{j} = @(f,a,A) A(:)-(f{j}*a);
    %SIRMin{j} = @(f,b,SI) SI(:)-(f{j}*b);
    %[a,fvala] = fmincon('minKS1802',igAge,f,A,[],[],lbAge,[1 1 1]);
    %[b,fvalb] = fmincon('minKS1802',igSIR,f,SIRKC49,[],[],[],[]);

    
    %Calculate linear mixing model to further constrain the young component
    %age.

Err=ones(length(A));
for k=1:length(A)
    D14C(k)=(A(k)-1)*1000;              %Calculate D14C, not age-corrected
    D14CErr(k)=1000*R(k,4);
    Ncum=0;
end
for n=1:length(A)
    Ncum=Ncum+R(n,2);     %Cumulative micromoles yield
    NmInv(n,1)=1/Ncum;                %Inverse cumulative yield
    Err(n,n)=1/D14CErr(n);
end



%Calculate linear model between D14C and NmInv.
AA=ones(length(NmInv),2);
AA(:,2)=NmInv(:);
BB=D14C;
X=(AA'*AA)^-1*AA'*BB;
Model=X(1)+X(2)*AA(:,2);

%Calculate error weighted model between D14C and NmInv
AAe=Err*AA;
BBe=Err*BB;
Xe=(AAe'*AAe)^-1*AAe'*BBe;
Modele=Xe(1)+Xe(2)*AA(:,2);


figure(703)
errorbar(NmInv,D14C,D14CErr,'ko');
hold on
%plot(NmInv,Model,'k');     %Plot linear Model
plot(NmInv,Modele,'k--');   %Plot weighted linear model
xlabel('Inverse Cumulative Yield','FontSize',14)
ylabel('\Delta^{14}C (permille modern)','FontSize',14)
T1=['\Delta^{14}C = ',num2str(roundn(Xe(1)),3),'+ ',num2str(roundn(Xe(2)),4),'(ICY)'];
T2=['r^{2} = ',num2str(roundn(corr2(D14C,NmInv)),2)];
TT=char(T1,T2);
text(0.01,-750,TT)




