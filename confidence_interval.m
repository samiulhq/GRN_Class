%determine confidence region for parameters

%load bootstarp data

%load result1.mat

%geneExpression = GeneRegulatorModel(12,[-10 0 0 0.001 0 0 -10 0 0 -10 0 0],[10 10 1 0.4 10 1 10 10 1 10 10 1]);

%[par, likelihood]=geneExpression.estimateParam();
%at first find range for 12 params in 95% confidence interval consdier all
%gaussian rng

range=zeros(2,12);
for i=1:12
    range(1,i)=mean(theta(:,i))- 1.96*std(theta(:,i));
    range(2,i)=mean(theta(:,i))+ 1.96*std(theta(:,i));
end
max=zeros(1,145);
min=100*ones(1,145);
%now search param space simulate only those satisfying 95% CI for all params
count=0;
figure;
geneNum=4;
for i=1:1000
    chk = theta(i,  :) >= range(1,:) & theta(i,:)<=range(2,:);
    if sum(chk) ==12 
    [T, tc]=geneExpression.simulate(theta(i,:));
   % plot(tc,T(:,4),'.','color','red');hold on;
   % count=count+1;
    a=T(:,geneNum)';  
    ind=a>max;
    max((find(ind==1)))=a((find(ind==1)));
    ind=a<min;
    min((f/                      ind(ind==1)))=a((find(ind==1)));
    else
    %    [T, tc]=geneExpression.simulate(theta(i,:));
    %    plot(tc,T(:,4),'.','color','black');hold on;
    end
end
load GeneData.mat
geneNames={'COL4','ETF9','ASIL2','MYB55'};
timeCourseTable=GeneData.TimeCourse;
sd=GeneData.Standard_Deviation;
h=errorbar([0 3 6 12 24 48 72],timeCourseTable(geneNum,:),sd(geneNum,:),'Marker','o','color','r','LineStyle','none');hold on
plot(tc,min,'color','red');hold on;
plot(tc,max,'color','red');hold on;
[T, tc]=geneExpression.simulate(par);
plot(tc,T(:,geneNum),'color','blue');
xlabel('Time');
ylabel(geneNames{geneNum});
plt=Plot();
plt.LineStyle{1}='none';
plt.Colors{2}=plt.Colors{3}
plt.export(strcat(geneNames{geneNum},'.png'));