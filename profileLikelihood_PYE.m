%profile likelihood

clc;
close all;
clear all;
%parpool
maxVal=10;
dx=0.1;

paramNum=[1 2 3 4];
val=0:dx:maxVal;
likelihood=zeros(length(paramNum),length(val));
par=zeros(length(paramNum),length(val),length(paramNum));

geneExpression = GeneRegulatorModelTargetTest(4,[-10 0 0 0.001],[10 10 1 0.4]);

[p,l]=geneExpression.estimateParam();
geneExpression.optimumParam=p;

for j=1:length(paramNum)
   j 
   
   if mod(j,2)==1
       val=0:dx:maxVal;
   else
       val=0:dx:maxVal;
   end
   
   len=length(val);
parfor i=1:len
i
geneExpression = GeneRegulatorModelTargetTest(4,[-10 0 0 0.001],[10 10 1 0.4],[paramNum(j)],val(i));

[par(j,i,:), likelihood(j,i)]=geneExpression.estimateParam();

end

end



for i=1:4
    subplot(2,2,i)
    if mod(i,2)==1
       val=0:dx:maxVal;
   else
       val=0:dx:maxVal;
    end
   plot(val,likelihood(i,:));
end

%save 'ProfileLikelihood_PYE_mArray.mat'
save 'ProfileLikelihood_bHLH115_mArray.mat'
