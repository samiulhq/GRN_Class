clc;
close all;
clear all;

load result1.mat
parname={'q_1','a_1','b_1','\lambda','a_2','b_2','q_3','a_3','b_3','q_4','a_4','b_4'};

for i=1:12
   subplot(3,4,i)
   histogram(theta(:,i),100);
  %[N,EDGE]=histcounts(theta(:,i),100);
 % N=N/sum(N);
  %E=(EDGE(1:end-1)+EDGE(2:end))/2;
  %plot(E,N);
   xlabel(parname{i});
   xlim([min(theta(:,i))-0.1 max(theta(:,i)+0.1)]); 
end
export_fig test.png -m5