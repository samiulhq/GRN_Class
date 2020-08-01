%profile likelihood

clc;
close all;
clear all;

%paramNum=[2 3 5 6 8 9 11 12];

default_step= 0.2;
min_step=10e-6;
max_step=0.2;

prof_like=cell(2,12);
GData=load('GeneData');
qPCRData=load('qpcr_regulators_plus_Fe');
geneExpression = GRM_class(12,[-10 0 0 0.001 0 0 -10 0 0 -10 0 0],[10 10 1 0.4 10 1 10 10 1 10 10 1],GData,qPCRData);
[optimP,minLike]=geneExpression.estimateParam();
%parpool
poolobj = gcp;
addAttachedFiles(poolobj,{'profile_calc.m','GRM_class.m','check_step_size'});


for j=1:4
    j
[likelihood,theta]=profile_calc(GData,optimP,minLike,j,default_step,qPCRData);
prof_like(j,:)={likelihood,theta};
end

save('updated_PL_calc.mat','prof_like')

