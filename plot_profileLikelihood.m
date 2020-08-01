close all;
%clear all;
load ProfileLikelihood_mArray.mat
%load ProfileLikelihood_mArray.mat
%thresh=min(min(likelihood))+0.004;
%thresh=min(min(likelihood))+chi2inv(0.80,1);

gm=GeneRegulatorModel2(12,[-10 0 0 0.001 0 0 -10 0 0 -10 0 0],[10 10 1 0.4 10 1 10 10 1 10 10 1]);
%gm=GeneRegulatorModel(12,[-10 0 0 0.001 0 0 -10 0 0 -10 0 0],[10 10 1 0.4 10 1 10 10 1 10 10 1]);
[oppar, lh]=gm.estimateParam();
paramNum=[2 3 5 6 8 9 11 12];
parname={'a_1','b_1','a_2','b_2','a_3','b_3','a_4','b_4'}
figure;

for i=1:8
    thresh=lh+chi2inv(0.80,1);
    subplot(4,2,i)
    plot(val,likelihood(i,:));
    hold on;
    plot(val,ones(1,length(val))*thresh,'--');
    hold on
    plot(oppar(paramNum(i)),lh,'Marker','x','MarkerSize',5);
    xlabel(parname{i});
    ylabel('Cost');
%     if i<7
%     xlim([0 3])
%     elseif i==7
%     xlim([0 3]);
%     else
%     xlim([0 8]);
%     end
%    
%     if i==4 || i==6
%     xlim([0 10]);
%     end
    ylim([0 thresh+1]);
    
end

%plt=Plot();
export_fig profileLikelihood_mArray_only.png -m5