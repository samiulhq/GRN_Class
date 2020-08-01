%profile likelihood

clc;
close all;
clear all;
%parpool
paramNum=[2 3 5 6 8 9 11 12];
minVal=-5;
dVal=0.1;
maxVal=20;
val=minVal:dVal:maxVal;
likelihood=zeros(length(paramNum),length(val));
par=zeros(length(paramNum),length(val),12);

for j=1:length(paramNum)
   j 
   
   if mod(j,2)==1
       val=minVal:dVal:maxVal;
   else
       val=minVal:dVal:maxVal;
   end
   
parfor i=1:length(val)
i
geneExpression = GeneRegulatorModel(12,[-10 0 0 0.001 0 0 -10 0 0 -10 0 0],[10 10 1 0.4 10 1 10 10 1 10 10 1],[paramNum(j)],val(i));

[par(j,i,:), likelihood(j,i)]=geneExpression.estimateParam();

end

end

save 'ProfileLikelihood_2.mat'

for i=1:8
    subplot(2,4,i)
    if mod(i,2)==1
       val=minVal:dVal:maxVal;
   else
       val=minVal:dVal:maxVal;
    end
   plot(val,likelihood(i,:));
end


