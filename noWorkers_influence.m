clear
close all


%EbNo = 10;
figure(3);


no_workers_set=[25, 50, 75, 100];
for rr=1:4
    
transmissionTime=1E-3;
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
        
%num_feature=size(Xdata_28,2);
%total_sample=4000;%size(Xdata_28,1);
%no_workers = 100;

% no_workers = 40;
% dataSamples_per_worker=400;%floor(total_sample/no_workers);



no_workers = no_workers_set(rr);
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



X_fede=[];
y_fede=[];
for i=1:no_workers
  X_fede=[X_fede;X{i}];
  y_fede=[y_fede;y{i}];
end


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

lambda=0.000;
bitsToSend=3;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%[path0, pathCost0, grid]=findPath2(num_workers);
% non-quantized algorithms
Rate=32*num_feature/transmissionTime;
sysBand=1E6;
Band=sysBand/(num_workers/2);
%Band=1E6;
[path, pathCost, d_square, P_central_uplink, P_central_downlink, center]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink=sum(P_central_uplink);
cost_star_topology_downlink=max(P_central_downlink);
cost_star_topology=cost_star_topology_uplink+cost_star_topology_downlink;
cost_decentralized=sum(pathCost);

%quantized algorithms
Rate=(32+num_feature*3)/transmissionTime;
%Band=1E6;
[~, pathCost_Q, ~, P_central_uplink_Q, P_central_downlink_Q, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q=sum(P_central_uplink_Q);
cost_star_topology_downlink_Q=max(P_central_downlink_Q);
cost_star_topology_Q=cost_star_topology_uplink_Q+cost_star_topology_downlink_Q;
cost_decentralized_Q=sum(pathCost_Q);


rho=20;%rho_set(rr);


    [obj_GADMM, loss_GADMM, Iter_gadmm, feas_gadmm] = group_ADMM_closedForm_2...
        (X_fede,y_fede, rho, no_workers, num_feature, dataSamples_per_worker,...
        num_iter, obj0, acc, alpha);
    rho
    
    for iter=1:Iter_gadmm%num_iter%Iter_gadmm
        cumulative_com_GADMM(iter)=iter*no_workers*num_feature*32;
    end

for iter=1:Iter_gadmm
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_gadmm(iter)=cost_decentralized*transmissionTime;
    else
        cumulative_EnergyCost_gadmm(iter)= cumulative_EnergyCost_gadmm(iter-1)+cost_decentralized*transmissionTime;
    end        
end
 


    [obj_2GADMM, loss_2GADMM, Iter_2GADMM, cumulative_com_2GADMM, feas_2gadmm] = quantized_group_ADMM...
        (X_fede,y_fede, rho, no_workers, num_feature, dataSamples_per_worker,...
        num_iter, obj0, acc,alpha,bitsToSend);
    

    

for iter=1:Iter_2GADMM
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgadmm(iter)=cost_decentralized_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgadmm(iter)= cumulative_EnergyCost_qgadmm(iter-1)+cost_decentralized_Q*transmissionTime;
    end        
end


semilogy(cumulative_com_GADMM, loss_GADMM,'LineWidth',3);
hold on
semilogy(cumulative_com_2GADMM, loss_2GADMM,'LineWidth',3);


clearvars -except rr no_workers_set 
end

xlabel({'Total number of transmitted bits'},'fontsize',16,'fontname','Times New Roman')
ylabel('Loss','fontsize',16,'fontname','Times New Roman')
legend('GADMM-25 workers', 'Q-GADMM-3bits-25 workers',...
    'GADMM-50 workers', 'Q-GADMM-3bits-50 workers'...
    ,'GADMM-75 workers', 'Q-GADMM-3bits-75 workers'...
    ,'GADMM-100 workers', 'Q-GADMM-3bits-100 workers');
set(gca,'fontsize',14,'fontweight','bold');
xlim([0 5.5E8])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

