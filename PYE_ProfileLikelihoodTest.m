%load ProfileLikelihood_PYE_mArray.mat
load 'ProfileLikelihood_bHLH115_mArray.mat'
close all;


geneExpression = GeneRegulatorModelTargetTest(4,[-10 0 0 0.001],[10 10 1 0.4]);

[oppar, lh]=geneExpression.estimateParam();

thresh=lh+chi2inv(0.8,1);
paramNum=[1 2 3 4];
parname={'q_1','a_1','b_1','\lambda'}
for i=1:4
    
    subplot(2,2,i)
    plot(val,likelihood(i,:));
    hold on;
    plot(val,ones(1,length(val))*thresh,'--');
    hold on
    plot(oppar(paramNum(i)),lh,'Marker','x','MarkerSize',5);
    xlabel(parname{i});
    ylabel('CostFunction');
%     if i<7
%     xlim([0 0.5])
%     elseif i==7
%     xlim([0 2]);
%     else
%     xlim([0 8]);
%     end
  % if i==1 
 %      xlim([0 1]);
  % if i==2
%       xlim([2 4]);
   %end
   
   %ylim([thresh-0.5 thresh+0.1]);
    
end

%plt=Plot();
export_fig profileLikelihood_mArray_bHLH115.png -m5