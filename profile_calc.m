
function [likelihood,theta]=profile_calc(GData,optimP,minLike,j,default_step,qPCRData)
 
theta=[];
likelihood=[];
theta(1,:)=optimP;
likelihood(1)=minLike;

step_size=default_step;

i=1;
%increasing theta
thresh = minLike + chi2inv(0.95,1);
current_p=theta(i,:);
current_like=likelihood(1);

direction=1;%increasing
steps=0;
while current_like<thresh
parNum=j;
%i
[flag, current_p, current_like, step_size]=check_step_size(step_size,current_p,current_like,parNum,direction,GData,qPCRData);
i=i+1;
likelihood(i)=current_like;
theta(i,:)=current_p;
if steps == 100
    break;
end
steps=steps+1;
end


current_p=theta(1,:);
current_like=likelihood(1);
direction=2; %decreasing
step_size=default_step;
steps=0;
while current_like<thresh
    
parnum=j;
%i

[flag, current_p, current_like, step_size]=check_step_size(step_size,current_p,current_like,parNum,direction,GData,qPCRData);
i=i+1;
likelihood(i)=current_like;
theta(i,:)=current_p;
if steps == 100
    break;
end
steps=steps+1;
end

[l b]=sort(theta(:,j));
theta=theta(b,:);
likelihood=likelihood(b);


end
