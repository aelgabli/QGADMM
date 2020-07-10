clear
close all


%EbNo = 10;

transmissionTime=1E-3;%100E-6;
[Xdata_28] = load('data99/data.txt'); 
[ydata_28] = load('data99/y.txt'); 

%[Xdata_28] = load('data28/data.txt'); % UCLA Housing dataset
%[ydata_28] = load('data28/y.txt'); 

%[Xdata_28] = load('data29/data.txt'); % Body Fat dataset
%[ydata_28] = load('data29/y.txt'); 


for i=1:size(Xdata_28,2)
    maxValue=max(Xdata_28(:,i));
    for j=1:size(Xdata_28,1)
       Xdata_28(j,i)= Xdata_28(j,i)/maxValue;
    end
end

maxValue=max(ydata_28);
for i=1:length(ydata_28)
    ydata_28(i)=ydata_28(i)/maxValue;
end

num_feature=size(Xdata_28,2);
total_sample=size(Xdata_28,1);
bitsToSend=2; 


%num_feature=size(Xdata_28,2);
%total_sample=4000;%size(Xdata_28,1);
%no_workers = 100;

% no_workers = 40;
% dataSamples_per_worker=400;%floor(total_sample/no_workers);

no_workers = 50;
dataSamples_per_worker=floor(total_sample/no_workers);


total_sample =no_workers*dataSamples_per_worker;

num_iter=200000;
%tot_virtual_workers=24;
X=cell(no_workers);
y=cell(no_workers);
num_workers=no_workers;
diagmatrix=diag(ones(total_sample,1));
% [lambda]=eig(Xdata'*Xdata);
Hmax=zeros(num_workers,1);

for n=1:no_workers
        first = (n-1)*dataSamples_per_worker+1;
        last = first+dataSamples_per_worker-1;
        X{n}=Xdata_28(first:last,1:num_feature);
        y{n}=ydata_28(first:last);
        Hmax(n)=max(eig(X{n}'*X{n}));
        
end


num_feature=size(X{1},2);

Hmax_sum=sum(Hmax);
hfun=Hmax_sum./Hmax;
nonprob=Hmax/Hmax_sum;

Hmin=zeros(num_workers,1);
Hcond=zeros(num_workers,1);
for i=1:num_workers
   Hmin(i)=min(abs(eig(X{i}'*X{i}))); 
   Hcond(i)=Hmax(i)/Hmin(i);
end



X_fede=[];
y_fede=[];
for i=1:no_workers
  X_fede=[X_fede;X{i}];
  y_fede=[y_fede;y{i}];
end

triggerslot=10;
Hmaxall=max(eig(X_fede'*X_fede)); 
Hminall=min(eig(X_fede'*X_fede)); 
condnum=Hmaxall/Hminall;
[cdff,cdfx] = ecdf(Hmax*num_workers/Hmaxall);
comm_save=0;
for i=1:triggerslot
    comm_save=comm_save+(1/i-1/(i+1))*cdff(find(cdfx>=min(max(cdfx),2*sqrt(1/(triggerslot*i))),1));
end


lambda=0.000;

%% Optimal solution

XX=X_fede;
YY=y_fede;
size = num_feature;

% obj0 = opt_sol(XX,YY, size);%(0.5*sum_square(XX*z1 - YY));
obj0 = opt_sol_closedForm(XX,YY);%(0.5*sum_square(XX*z1 - YY));
%obj1 = opt_sol(XX,YY,num_feature);%(0.5*sum_square(XX*z1 - YY));

opt_obj = obj0*ones(1,num_iter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

acc = 1E-4;
alpha=1;


num_workers=no_workers;
%%----------------------------------------------------------------
heterconst=mean(exp(Hmax/Hmaxall));
heterconst2=mean(Hmax/Hmaxall);
rate=1/(1+sum(Hmin)/(4*sum(Hmax)));
%% parameter initialization
%triggerslot=100;
theta_agent_1=zeros(num_feature,num_iter);
theta_agent_2=zeros(num_feature,num_iter);
theta_agent_3=zeros(num_feature,num_iter);

theta=zeros(num_feature,num_iter);
grads=ones(num_feature,num_workers);
%stepsize=1/(num_workers*max(Hmax));
stepsize=1/Hmaxall;
thrd=10/(stepsize^2*num_workers^2)/triggerslot;
comm_count=ones(num_workers,1);

theta2=zeros(num_feature,num_iter);
grads2=zeros(num_feature,1);
q_grad_w=1/num_workers*ones(num_feature,num_workers);
q_grad_laq=1/num_workers*ones(num_feature,num_workers);
grad_laq=1/num_workers*ones(num_feature,num_workers);
grads2=sum(q_grad_w,2);%(num_feature,1);
stepsize2=stepsize;
%stepsize2=1E-4;

theta3=zeros(num_feature,num_iter);
grads3=ones(num_feature,num_workers);
stepsize3=stepsize2/num_workers; % cyclic access learning

theta4=zeros(num_feature,num_iter);
grads4=ones(num_feature,num_workers);
stepsize4=stepsize/num_workers; % nonuniform-random access learning


thrd5=1/(stepsize^2*num_workers^2)/triggerslot;
theta5=zeros(num_feature,1);
grads5=ones(num_feature,num_workers);
stepsize5=stepsize;
%stepsize5=1E-4;
comm_count5=ones(num_workers,1);

%thrd6=2/(stepsize*num_workers);
theta6=zeros(num_feature,1);
grads6=ones(num_feature,1);
stepsize6=0.5*stepsize;
comm_count6=ones(num_workers,1);


theta7=zeros(num_feature,1);
grads7=ones(num_feature,num_workers);
stepsize7=stepsize;
comm_count7=ones(num_workers,1);

lambda=0.000;

%%%%%%%%%%%%%%%%%%%%%%%%%%

%[path0, pathCost0, grid]=findPath2(num_workers);
% non-quantized algorithms
Rate=32*num_feature/transmissionTime;%3.84E6;%needs to transmit 32*d bits in trabsmsiiionTime seconds
sysBand=2E6;%20E6;
Band=sysBand/(num_workers/2);
%Band=1E6;
[path, pathCost, d_square, P_central_uplink, P_central_downlink, center]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink=sum(P_central_uplink);
cost_star_topology_downlink=max(P_central_downlink);
cost_star_topology=cost_star_topology_uplink+cost_star_topology_downlink;
cost_decentralized=sum(pathCost);

%quantized algorithms
Rate=(32+num_feature*bitsToSend)/transmissionTime; %1E6; needs to transmit 32+6*3 bits in (32+2*num_feature*3)/transmissionTime;
%Band=1E6;
[~, pathCost_Q, ~, P_central_uplink_Q, P_central_downlink_Q, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q=sum(P_central_uplink_Q);
cost_star_topology_downlink_Q=max(P_central_downlink_Q);
cost_star_topology_Q=cost_star_topology_uplink_Q+cost_star_topology_downlink_Q;
cost_decentralized_Q=sum(pathCost_Q);

%A-DIANA 
Rate=(32+2*num_feature*bitsToSend)/transmissionTime;% 1.36E6;%needs to transmit 32+2*6*3 bits in (32+2*num_feature*3)/transmissionTime;
%Band=1E6;
[~, pathCost_Q_diana, ~, P_central_uplink_Q_diana, P_central_downlink_Q_diana, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q_diana=sum(P_central_uplink_Q_diana);
cost_star_topology_downlink_diana=max(P_central_downlink);
cost_star_topology_Q_diana=cost_star_topology_uplink_Q_diana+cost_star_topology_downlink_diana;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  ADIANA
comm_error2=[];
comm_grad2=[];
flagg=0;
s1=num_feature;
s2=dataSamples_per_worker;
q_grad_w=1/num_workers*ones(num_feature,num_workers);
theta2=zeros(num_feature,num_iter);
W=theta2;
y=theta2;
z=theta2;
h=zeros(num_feature,num_workers);
h_central=zeros(num_feature,1);
s=2^bitsToSend;
omega=2+2*sqrt(num_feature)/s;
alpha_s=1/(omega+1);
eta=num_workers/(64*omega*(4*(omega+1)+1)^2);
if(eta > 0.5)
    eta=0.5;
end
zeta_1=sqrt(eta/2);
if(zeta_1 > 0.25)
    zeta_1=0.25;
end

zeta_2=0.5;
gamma=eta/(2*(zeta_1+eta));
beta=1-gamma;
temp1=sqrt(no_workers/(32*omega))-1;
if(temp1 < 1)
    temp1=1;
end
prob=temp1/(2*(1+omega));
if (prob > 1)
    prob =1;
end

eta=1E-4;
gamma=9*eta;
for iter=1:num_iter
    
    % central server computation
    theta2(:,iter)=zeta_1*z(:,iter)+zeta_2*W(:,iter)+(1-zeta_1-zeta_2)*y(:,iter);
    %theta2(:,iter)=zeta_2*W(:,iter)+(1-zeta_1-zeta_2)*y(:,iter);
    
    %if iter>1
    agg_gradient=zeros(num_feature,1);
    agg_h=zeros(num_feature,1);
        for i=1:num_workers
            first = (i-1)*s2+1;
            last = first+s2-1;
            grad_w(:,i)=X_fede(first:last,1:s1)'*X_fede(first:last,1:s1)*theta2(:,iter)-X_fede(first:last,1:s1)'*y_fede(first:last);
            [q_G(:,i),number_of_bits_toSend]=ADIANA_Q(grad_w(:,i)-h(:,i),zeros(num_feature,1),bitsToSend);
            agg_gradient=agg_gradient+q_G(:,i);
            
            grad_w2(:,i)=X_fede(first:last,1:s1)'*X_fede(first:last,1:s1)*W(:,iter)-X_fede(first:last,1:s1)'*y_fede(first:last);
            [q_G_W2(:,i),~]=ADIANA_Q(grad_w2(:,i)-h(:,i),zeros(num_feature,1),bitsToSend);
            agg_h=agg_h+q_G_W2(:,i);
            
            %q_grad_w(:,i)=q_grad_w(:,i)+q_delta(:,i);
            %grads2=sum(q_grad_w(:,i));%q_delta(:,i);
        end
        
        agg_gradient = agg_gradient + h_central;%1/num_workers*
        
        h_central=h_central+alpha_s*agg_h;%*1/num_workers;
    %end

        obj_diana(iter)=0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2;
        loss_diana(iter)=abs(0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2-obj0);%+0.5*num_workers*lambda*norm(theta2(:,iter),2)^2;
        toPrintOut_diana=loss_diana(iter)
        
        if(loss_diana(iter) < acc)
            iter_diana=iter;
            break;
        end
        %theta2(:,iter+1)=theta2(:,iter)-stepsize2*sum(q_G,2);
        y(:,iter+1)=theta2(:,iter)-eta*agg_gradient;
        z(:,iter+1)= beta*z(:,iter)+(1-beta)*theta2(:,iter)+gamma/eta*(y(:,iter+1)-theta2(:,iter));
        randNum=rand;
        if(randNum <=prob)
            W(:,iter+1)=y(:,iter);
        else
            W(:,iter+1)=W(:,iter);
        end
        %theta2(:,iter+1)=theta2(:,iter)-stepsize2*sum(q_grad_w,2);


        
      for i=1:num_workers
           h(:,i) = h(:,i) +alpha_s*q_G_W2(:,i);
      end
        
    
end

for iter=1:iter_diana
    cumulative_com_rounds_diana(iter)=iter*num_workers+iter; 
    cumulative_com_bits_diana(iter)=iter*num_workers*(32+s1*2*bitsToSend)+2*iter*32*s1; 
    %errorPer_GD(iter) = abs(loss_GD(iter)/opt_obj(iter)*100);
end

%num_iter=iter_diana;

for iter=1:iter_diana
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_diana(iter)=cost_star_topology_Q_diana*transmissionTime;
    else
        cumulative_EnergyCost_diana(iter)= cumulative_EnergyCost_diana(iter-1)+cost_star_topology_Q_diana*transmissionTime;
    end        
end



%%  QGD
comm_error2=[];
comm_grad2=[];
flagg=0;
s1=num_feature;
s2=dataSamples_per_worker;

for iter=1:num_iter

    % central server computation
    if iter>1
        for i=1:num_workers
            first = (i-1)*s2+1;
            last = first+s2-1;
            grad_w(:,i)=X_fede(first:last,1:s1)'*X_fede(first:last,1:s1)*theta2(:,iter)-X_fede(first:last,1:s1)'*y_fede(first:last);
            [q_delta(:,i),number_of_bits_toSend]=GD_Q(grad_w(:,i),q_grad_w(:,i),bitsToSend);
            q_grad_w(:,i)=q_grad_w(:,i)+q_delta(:,i);
            %grads2=sum(q_grad_w(:,i));%q_delta(:,i);
        end
    end
        %grad_error2(iter)=norm(sum(grads2,2),2);
        %obj_GD(iter)=0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2;
        obj_GD_2(iter)=0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2;
        loss_GD_2(iter)=abs(0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2-obj0);%+0.5*num_workers*lambda*norm(theta2(:,iter),2)^2;
        %theta2(:,iter+1)=theta2(:,iter)-stepsize2*grads2;
        theta2(:,iter+1)=theta2(:,iter)-stepsize2*sum(q_grad_w,2);
        %comm_error2=[comm_error2;iter*num_workers,loss_GD(iter)]; 
        %comm_grad2=[comm_grad2;iter*num_workers,grad_error2(iter)];
        toPrintOut_qgd=loss_GD_2(iter)
        
        if(loss_GD_2(iter) < acc)
            iter_GD_2=iter;
            break;
        end
        
    
end

for iter=1:iter_GD_2
    cumulative_com_rounds_GD_2(iter)=iter*num_workers+iter; 
    cumulative_com_bits_GD_2(iter)=iter*num_workers*(32+s1*bitsToSend)+iter*32*s1; 
    %errorPer_GD(iter) = abs(loss_GD(iter)/opt_obj(iter)*100);
end

%num_iter=iter_GD_2;

for iter=1:iter_GD_2
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgd(iter)=cost_star_topology_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgd(iter)= cumulative_EnergyCost_qgd(iter-1)+cost_star_topology_Q*transmissionTime;
    end        
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%  stochastic-GD
% comm_error2=[];
% comm_grad2=[];
% flagg=0;
% s1=num_feature;
% s2=dataSamples_per_worker;
% q_grad_w=1/num_workers*ones(num_feature,num_workers);
% theta2=zeros(num_feature,num_iter);
% 
% for iter=1:num_iter
% 
%     % central server computation
%     if iter>1
%         for i=1:num_workers
%             first = (i-1)*s2+1;
%             last = first+s2-1;
%             grad_w(:,i)=X_fede(first:last,1:s1)'*X_fede(first:last,1:s1)*theta2(:,iter)-X_fede(first:last,1:s1)'*y_fede(first:last);
%             [q_delta(:,i),number_of_bits_toSend]=GD_stochasticQ(grad_w(:,i),q_grad_w(:,i),bitsToSend);
%             q_grad_w(:,i)=q_grad_w(:,i)+q_delta(:,i);
%             %grads2=sum(q_grad_w(:,i));%q_delta(:,i);
%         end
%     end
%         %grad_error2(iter)=norm(sum(grads2,2),2);
%         %obj_GD(iter)=0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2;
%         obj_GD_s(iter)=0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2;
%         loss_GD_s(iter)=abs(0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2-obj0);%+0.5*num_workers*lambda*norm(theta2(:,iter),2)^2;
%         %theta2(:,iter+1)=theta2(:,iter)-stepsize2*grads2;
%         theta2(:,iter+1)=theta2(:,iter)-stepsize2*sum(q_grad_w,2);
%         %comm_error2=[comm_error2;iter*num_workers,loss_GD(iter)]; 
%         %comm_grad2=[comm_grad2;iter*num_workers,grad_error2(iter)];
%         toPrintOut_qgd=loss_GD_s(iter)
%         
%         if(loss_GD_s(iter) < acc)
%             iter_GD_s=iter;
%             break;
%         end
%         
%     
% end
% 
% for iter=1:iter_GD_s
%     cumulative_com_rounds_GD_s(iter)=iter*num_workers+iter; 
%     cumulative_com_bits_GD_s(iter)=iter*num_workers*(32+s1*bitsToSend)+iter*32*s1; 
%     %errorPer_GD(iter) = abs(loss_GD(iter)/opt_obj(iter)*100);
% end
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho=24;


    [obj_GADMM, loss_GADMM, Iter_gadmm, feas_gadmm, comp_time_gadmm] = group_ADMM_closedForm_2...
        (X_fede,y_fede, rho, no_workers, num_feature, dataSamples_per_worker,...
        num_iter, obj0, acc, alpha);
    rho
    
    for iter=1:Iter_gadmm%num_iter%Iter_gadmm
        cumulative_com_GADMM(iter)=iter*no_workers*num_feature*32;
    end
    
    for iter=1:Iter_gadmm%num_iter%Iter_gadmm
        cumulative_com_rounds_GADMM(iter)=iter*no_workers;
    end
    
    %num_iter=Iter_2GADMM;

for iter=1:Iter_gadmm
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_gadmm(iter)=cost_decentralized*transmissionTime;
    else
        cumulative_EnergyCost_gadmm(iter)= cumulative_EnergyCost_gadmm(iter-1)+cost_decentralized*transmissionTime;
    end        
end
    

    [obj_2GADMM, loss_2GADMM, Iter_2GADMM, cumulative_com_2GADMM, feas_2gadmm, comp_time_qgadmm] = quantized_group_ADMM...
        (X_fede,y_fede, rho, no_workers, num_feature, dataSamples_per_worker,...
        num_iter, obj0, acc,alpha,bitsToSend);
    
    for iter=1:Iter_2GADMM%num_iter%Iter_gadmm
        cumulative_com_rounds_2GADMM(iter)=iter*no_workers;
    end
    
    
%num_iter=Iter_2GADMM;

for iter=1:Iter_2GADMM
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgadmm(iter)=cost_decentralized_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgadmm(iter)= cumulative_EnergyCost_qgadmm(iter-1)+cost_decentralized_Q*transmissionTime;
    end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GD
comm_error2=[];
comm_grad2=[];
flagg=0;
s1=num_feature;
s2=dataSamples_per_worker;
%bitsToSend=3;
for iter=1:num_iter

    % central server computation
    if iter>1
        for i=1:num_workers
            first = (i-1)*s2+1;
            last = first+s2-1;
            grad_w(:,i)=X_fede(first:last,1:s1)'*X_fede(first:last,1:s1)*theta2(:,iter)-X_fede(first:last,1:s1)'*y_fede(first:last);
            %[q_delta(:,i),number_of_bits_toSend]=laq_quantizer(grad_w(:,i),q_grad_w(:,i),bitsToSend);
            %q_grad_w(:,i)=q_grad_w(:,i)+q_delta(:,i);
            %grads2=sum(q_grad_w(:,i));%q_delta(:,i);
        end
    end
        %grad_error2(iter)=norm(sum(grads2,2),2);
        %obj_GD(iter)=0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2;
        obj_GD(iter)=0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2;
        loss_GD(iter)=abs(0.5*norm(X_fede*theta2(:,iter)-y_fede,2)^2-obj0);%+0.5*num_workers*lambda*norm(theta2(:,iter),2)^2;
        %theta2(:,iter+1)=theta2(:,iter)-stepsize2*grads2;
        theta2(:,iter+1)=theta2(:,iter)-stepsize2*sum(grad_w,2);
        %comm_error2=[comm_error2;iter*num_workers,loss_GD(iter)]; 
        %comm_grad2=[comm_grad2;iter*num_workers,grad_error2(iter)];
        
        if(loss_GD(iter) < acc)
            iter_GD=iter;
            break;
        end
        
    
end

for iter=1:iter_GD
    cumulative_com_rounds_GD(iter)=iter*num_workers+iter; 
    cumulative_com_bits_GD(iter)=iter*num_workers*s1*32+iter*32*s1; 
    %errorPer_GD(iter) = abs(loss_GD(iter)/opt_obj(iter)*100);
end

for iter=1:iter_GD
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_gd(iter)=cost_star_topology*transmissionTime;
    else
        cumulative_EnergyCost_gd(iter)= cumulative_EnergyCost_gd(iter-1)+cost_star_topology*transmissionTime;
    end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
subplot(1,3,1);
semilogy(cumulative_com_rounds_GD, loss_GD,'b-','LineWidth',3);
hold on
semilogy(cumulative_com_rounds_GD_2, loss_GD_2,'m--','LineWidth',3);
hold on
semilogy(cumulative_com_rounds_diana, loss_diana,'c--','LineWidth',3);
%semilogy(cumulative_com_rounds_GD_2, loss_GD_2,'c-','LineWidth',3);
semilogy(cumulative_com_rounds_GADMM, loss_GADMM,'k-','LineWidth',3);
hold on
semilogy(cumulative_com_rounds_2GADMM, loss_2GADMM,'r--','LineWidth',3);
hold on
%semilogy(cumulative_com_rounds_LAG, loss_LAG_WK,'b:','LineWidth',3);
%semilogy(loss_1GADMM,'m-','LineWidth',3);
%semilogy(loss_17GADMM,'k--','LineWidth',3);
xlabel({'Number of communication rounds';'(a)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Loss','fontsize',16,'fontname','Times New Roman')
legend('GD','QGD-2bits','ADIANA','GADMM', 'Q-GADMM-2bits');
%ylim([0 10^2])
%xlim([20 1800])
%xlim([50 6000])
%xlim([3000 6000])

set(gca,'fontsize',14,'fontweight','bold');

subplot(1,3,2);
semilogy(cumulative_com_bits_GD, loss_GD,'b-','LineWidth',3);
hold on
semilogy(cumulative_com_bits_GD_2, loss_GD_2,'m--','LineWidth',3);
hold on
semilogy(cumulative_com_bits_diana, loss_diana,'c--','LineWidth',3);
%semilogy(cumulative_com_bits_GD_2, loss_GD_2,'c-','LineWidth',3);
hold on
semilogy(cumulative_com_GADMM, loss_GADMM,'k-','LineWidth',3);
hold on
semilogy(cumulative_com_2GADMM, loss_2GADMM,'r--','LineWidth',3);
%semilogy(cumulative_com_bits_lag, loss_LAG_WK,'b:','LineWidth',3);
%semilogy(cumulative_com_1GADMM,'m-','LineWidth',3);
%semilogy(cumulative_com_17GADMM,'k--','LineWidth',3);
xlabel({'Total number of transmitted bits';'(b)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Loss','fontsize',16,'fontname','Times New Roman')
legend('GD','QGD-2bits','ADIANA','GADMM', 'Q-GADMM-2bits');
set(gca,'fontsize',14,'fontweight','bold');
%ylim([10^-4 10^3])
%xlim([10 6000])



subplot(1,3,3);
%semilogy(cumulative_EnergyCost_gd, loss_GD,'b-','LineWidth',3);
%hold on
semilogy(cumulative_EnergyCost_qgd, loss_GD_2,'m--','LineWidth',3);
hold on
semilogy(cumulative_EnergyCost_diana, loss_diana,'c--','LineWidth',3);
%semilogy(cumulative_com_bits_GD_2, loss_GD_2,'c-','LineWidth',3);
hold on
semilogy(cumulative_EnergyCost_gadmm, loss_GADMM,'k-','LineWidth',3);
hold on
semilogy(cumulative_EnergyCost_qgadmm, loss_2GADMM,'r--','LineWidth',3);
xlabel({'Sum energy';'(c)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Loss','fontsize',16,'fontname','Times New Roman')
legend('QGD-2bits','ADIANA','GADMM', 'Q-GADMM-2bits');
set(gca,'fontsize',14,'fontweight','bold');
%ylim([10^-4 10^3])
%xlim([10 6000])

figure(17);
semilogy(cumulative_com_rounds_GADMM, loss_GADMM,'k-','LineWidth',3);
hold on
semilogy(cumulative_com_rounds_2GADMM, loss_2GADMM,'r-','LineWidth',3);
legend('GADMM', 'Q-GADMM-2bits');
set(gca,'fontsize',14,'fontweight','bold');

figure(18);
semilogy(comp_time_gadmm, loss_GADMM,'k-','LineWidth',3);
hold on
semilogy(comp_time_qgadmm, loss_2GADMM,'r-','LineWidth',3);
xlabel({'Computation time'},'fontsize',16,'fontname','Times New Roman')
ylabel('Loss','fontsize',16,'fontname','Times New Roman')
legend('GADMM', 'Q-GADMM-2bits');
set(gca,'fontsize',14,'fontweight','bold');



% save diana_results_2bits_rho24.mat cumulative_EnergyCost_diana iter_diana...
%     cumulative_EnergyCost_qgd iter_GD_2...
%     cumulative_EnergyCost_gadmm Iter_gadmm...
%     cumulative_EnergyCost_qgadmm Iter_2GADMM

% save diana_results.mat cumulative_EnergyCost_diana iter_diana...
%     cumulative_EnergyCost_qgd iter_GD_2...
%     cumulative_EnergyCost_gadmm Iter_gadmm...
%     cumulative_EnergyCost_qgadmm Iter_2GADMM