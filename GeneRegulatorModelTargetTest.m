classdef GeneRegulatorModelTargetTest
    properties
        numOfParameters
        microArrayData 
        qPCRData
        initialParams
        lb
        ub
        optimumParam
        optimumlikelihood
        TimeC
       
    end
    
    methods
        
        function obj=set.optimumParam(obj,val)
            obj.optimumParam=val;
        end
        
        function obj=set.TimeC(obj,val)
            obj.TimeC=val;
        end
        function obj = GeneRegulatorModelTargetTest(varargin)
       lb=  varargin{2};
       ub= varargin{3};
          if length(varargin)==5 
               paramInd=varargin{4};
               paramVal=varargin{5};
               
              for i=1:length(paramInd)
                  lb(paramInd)=paramVal(i);
                  ub(paramInd)=paramVal(i);
                  
              end
               
           end
            obj.numOfParameters= varargin(1);
            obj.lb=lb;
            obj.ub=ub;
            obj.microArrayData=load('GeneData');
            obj.TimeC=obj.microArrayData.GeneData.TimeCourse;
            obj.qPCRData=load('qpcr_plus_Fe');
            
            obj.initialParams=(lb+ub)/2;
            
            
        end
        
        function ode=odefuncWithOutIron(obj,p)
            
            ode=@(t,x) [p(1)*(1-exp(-p(4)*t))+p(2)-p(3)*x(1); ...
                p(5)-p(6)*x(2);...
                p(7)*(1-exp(-p(4)*(t)))+p(8)-p(9)*x(3);...
                p(10)*(1-exp(-p(4)*(t)))+p(11)-p(12)*x(4)];
            
        end
        function ode=odefuncWithIron(obj,p)
            ode=@(t,x) [p(2)-p(3)*x(1); ...
                p(5)-p(6)*x(2);...
                p(8)-p(9)*x(3);...
                p(11)-p(12)*x(4)];
            
        end
        
        
        function likelihood=objective_func(obj,param,GeneData,qpcr)
                        
            
            TC=obj.TimeC;
            data1=TC(6,:);            
            
            p=param;
            f=@(t,x) [p(1)*(1-exp(-p(4)*(t)))+p(2)-p(3)*x(1)];
            
            [t,xa] = ode45(f,(0:0.5:72),[3.8912]);
            
            ind =[1 3 6 12 24 48 72]/0.5 - 1;
                        
            likelihood= (data1-xa(ind,1)')* (data1-xa(ind,1)')';                                
           
            
        end
        
        
        
        function [params, L]=estimateParam(obj)
            
           
            [r,fval,~ ,~]=fmincon(@(param)obj.objective_func(param,obj.microArrayData.GeneData,obj.qPCRData.qpcr_plus_Fe),obj.initialParams,[],[],[],[],obj.lb,obj.ub);
            params=r;
            L=fval;
            
            obj.optimumParam=r;
            obj.optimumlikelihood=L;
        end
        
        function [geneExpression,tc] = simulate(obj,p)
            
            f=@(t,x) [p(1)*(1-exp(-p(4)*(t)))+p(2)-p(3)*x(1)];
            
            [t,xa] = ode45(f,(0:0.5:72),[3.8912]);    %initial value for ode should be set intelligently %to do
            
            geneExpression=xa;
            tc=t;
            
            
        end
        
        
        function [geneExpression,tc] = simulate2(obj,p)
            f=@(t,x) [p(2)-p(3)*x(1)];
            
            [t,xa] = ode45(f,(0:0.5:72),[3.8912]);    %initial value for ode should be set intelligently %to do
            ind =[1 3 6 12 24 48 72]/0.5 - 1;
            geneExpression=xa;
            tc=t;
            
            
        end
        
        
        function plotresults(obj,geneName)
            geneNames={'COL4','ETF9','ASIL2','MYB55','MYB10','bHLH115','MYB72','bHLH39','BTS','PYE','bHLH101'};
            ind=strncmp(geneNames,geneName,length(geneName));
            num=find(ind>0);
            load qpcr_plus_Fe
            load qpcr_minus_Fe
          %  load(fname);
            timeCourseTable=obj.TimeC;
            
            figure;
                        
            sd=obj.microArrayData.GeneData.Standard_Deviation;
                        
            [result, tc]=obj.simulate(obj.optimumParam);
            result_without_ironsignal=obj.simulate2(obj.optimumParam);
                                   
            h=errorbar([0 3 6 12 24 48 72],timeCourseTable(num,:),sd(num,:),'Marker','o','color','r','LineStyle','none');
            
            hold on;
            plot(tc,result,'color','blue');
            
            plot([0 12 24 36],qpcr_plus_Fe(num,:),'x','color','black');
            plot(tc,result_without_ironsignal,'color','green');
            plot([0 12 24 36],qpcr_minus_Fe(num,:),'o','color','black');
            legend('measurement in -Fe','prediction in -Fe ','PCR measurement in +Fe','Prediction in +Fe','PCR measurement -Fe')
            xlabel('Time');
            ylabel(geneName);
            xlim([0 80]);
            plt=Plot();
            plt.LineStyle={'none','-','none','-','none'};
            plt.BoxDim=[10 6];
            
            
        end
    end
    
end