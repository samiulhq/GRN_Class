clc;
close all;
clear all;

geneExpression = GeneRegulatorModel(12,[-10 0 0 0.001 0 0 -10 0 0 -10 0 0],[10 10 1 0.4 10 1 10 10 1 10 10 1]);

[par, likelihood]=geneExpression.estimateParam();
 geneExpression.optimumParam=par;
% geneExpression.plotresults('COL4');
% geneExpression.plotresults('ETF9');
% geneExpression.plotresults('ASIL2');
% geneExpression.plotresults('MYB55');

theta=zeros(1000,12);
L=zeros(1,1000);
L(1)=likelihood;
theta(1,:)=par;

sd=geneExpression.microArrayData.GeneData.Standard_Deviation; 
 %  T=geneExpression.microArrayData.GeneData.TimeCourse;
   sd=sd(1:4,:);
%   T=T(1:4,:);
  [T, tc]= geneExpression.simulate(par);
  ind =[1 3 6 12 24 48 72]/0.5 - 1;
  T=T(ind,:)';
for it=2:1000
   
   boot_strap=normrnd(T,sd);
   geneExpression.TimeC=boot_strap;
   [theta(it,:), L(it)]=geneExpression.estimateParam();
end

%save('result1','theta',' L')
%plotting

%bootstarping 


