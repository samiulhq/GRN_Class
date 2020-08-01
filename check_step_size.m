function [flag, current_p, current_like ,step_size]=check_step_size(step_size,current_p,current_like,parNum,direction,GData,qPCRData)

if direction == 2
    if (step_size>0)
        step_size = (-1) * step_size;
    end
end

flag=1;
q=0.1;

proposed_par=current_p;
proposed_par(parNum)=proposed_par(parNum)+step_size;
gs=GRM_class(12,[-10 0 0 0.001 0 0 -10 0 0 -10 0 0],[10 10 1 0.4 10 1 10 10 1 10 10 1],GData,qPCRData,parNum,proposed_par(parNum));
thresh_alpha=chi2inv(0.95,1);
[proposed_par, proposed_like]=gs.estimateParam();

likelihood_diff = abs(current_like -proposed_like);
if likelihood_diff >q*thresh_alpha
    while (likelihood_diff >q*thresh_alpha)
        step_size=step_size/2;
        proposed_par=current_p;
        proposed_par(parNum)=proposed_par(parNum)+step_size;
        gs=GRM_class(12,[-10 0 0 0.001 0 0 -10 0 0 -10 0 0],[10 10 1 0.4 10 1 10 10 1 10 10 1],GData,qPCRData,parNum,proposed_par(parNum));
        [proposed_par, proposed_like]=gs.estimateParam();       
        likelihood_diff = abs(current_like -proposed_like);
       % if abs(step_size)<10e-6 
       %     break
       % end
    end
else
    while (likelihood_diff <q*thresh_alpha)
        step_size=step_size * 2;
        proposed_par=current_p;
        proposed_par(parNum)=proposed_par(parNum)+step_size;
         gs=GRM_class(12,[-10 0 0 0.001 0 0 -10 0 0 -10 0 0],[10 10 1 0.4 10 1 10 10 1 10 10 1],GData,qPCRData,parNum,proposed_par(parNum));
        [proposed_par, proposed_like]=gs.estimateParam();
       
        likelihood_diff = abs(current_like -proposed_like);
        %if abs(step_size)>0.2 
          %  break
        %end
    end
    step_size=step_size/2;
end

current_like=proposed_like;
current_p = proposed_par;
end