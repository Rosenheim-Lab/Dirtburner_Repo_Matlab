function [p_out, iter, model_out, RMS] = fit_Gauss(run_file_path, number_of_Gaussians, p_initial)

%Function to fit Gaussian curves to Ramped PyrOx run data. 
%
%This function fits the desired number of Gaussian curves to run data from the RPO system. It has
%dependencies on the following functions:
% dfdp
% nlleasq
%
%	Inputs:
%		run_file_path - (string) text file where the run data can be located. The data must be
%			in the following columns: Temperature in column 1, pCO2 in column 2. If not, 
%			a dummy file can be made with the data in the proper columns.
%		number_of_Gaussians - (int) number of Gaussians you believe to fit into the curve. 
%		p_initial - (vector) values for Gaussian variables Height, Center, Width. These values
%			should total to 3*number_of_Gaussians because they need to be specified for each
%			curve in the order above. Note that the width is the width of 1 standard deviation
%			thus, visually, you should provide guess that is roughly 1/3 of what you think
%			you see.
%
%Created by B.E. Rosenheim, University of South Florida 7 Oct 2021



%Test to ensure that there are enough p_initial values for the number of Gaussians specified
if ~isinf(number_of_Gaussians) & floor(number_of_Gaussians) ~= number_of_Gaussians
	error('Input "number_of_Gaussians" must be an integer.')
end

if length(p_initial) ~= number_of_Gaussians*3
	error('There must be 3x as many p_initial values as there are Gaussians in your model!')
end
	

%load the data files for data and for results.
D=sortrows(load(run_file_path),1);
T=D(:,1);
y=D(:,2);

pin = p_initial;

%Model initial guess:
for x=1:length(pin)/3
	Init_Models(:,x) = pin(x*3-2)*exp(-((T-pin(x*3-1))/pin(x*3)).^2);
end

%Take sum of each row into a column vector that is the sum of inital guess Guassians
ModelIg = sum(Init_Models, 2);

%Calculate the best fit of the number of specified Gaussians:
[f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2]=nlleasqr(T,y,pin,'Gauss_func',0.00001,100);

%Calculate the model reproduction of the curve and the root mean square
%of the residual between model and data
for x=1:length(p)/3
	Out_Models(:,x) = p(x*3-2)*exp(-((T-p(x*3-1))/p(x*3)).^2);
end

Modelout= sum(Out_Models, 2)
RMS=sqrt(sum((Modelout-y).^2)/length(y));

%Compute individual components of the model
for h=1:length(pin)/3
    Models(:,h)=p(3*h-2)*exp(-((T-p(3*h-1))/p(3*h)).^2);
end


       
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


%Output variables
model_out = Models;
p_out = p;
