clear
close all


for rr=1:100
[Xdata_28] = load('data99/data.txt'); 
[ydata_28] = load('data99/y.txt'); 
transmissionTime=1E-3;
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
no_workers = 50;
dataSamples_per_worker=floor(total_sample/no_workers);


total_sample =no_workers*dataSamples_per_worker;

num_iter=200000;
%tot_virtual_workers=24;
X=cell(no_workers);
y=cell(no_workers);

diagmatrix=diag(ones(total_sample,1));
% [lambda]=eig(Xdata'*Xdata);
%Hmax=zeros(num_workers,1);

for n=1:no_workers
        first = (n-1)*dataSamples_per_worker+1;
        last = first+dataSamples_per_worker-1;
        X{n}=Xdata_28(first:last,1:num_feature);
        y{n}=ydata_28(first:last);
        %Hmax(n)=max(eig(X{n}'*X{n}));
        
end


num_feature=size(X{1},2);







num_workers=no_workers;

bitsToSend=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%[path0, pathCost0, grid]=findPath2(num_workers);
Rate=32*num_feature/transmissionTime;%3.84E6;
sysBand=10E6;%20E6;
Band=sysBand/(num_workers/2);
[path, pathCost, d_square, P_central_uplink, P_central_downlink, center]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink=sum(P_central_uplink);
cost_star_topology_downlink=max(P_central_downlink);
cost_star_topology=cost_star_topology_uplink+cost_star_topology_downlink;
cost_decentralized=sum(pathCost);


Rate=(32+num_feature*bitsToSend)/transmissionTime;%1E6;
%Band=1E6;
[~, pathCost_Q, ~, P_central_uplink_Q, P_central_downlink_Q, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q=sum(P_central_uplink_Q);
cost_star_topology_downlink_Q=max(P_central_downlink);
cost_star_topology_Q=cost_star_topology_uplink_Q+cost_star_topology_downlink_Q;
cost_decentralized_Q=sum(pathCost_Q);

Rate=(32+2*num_feature*bitsToSend)/transmissionTime;%1.36E6;
%Band=1E6;
[~, pathCost_Q_diana, ~, P_central_uplink_Q_diana, P_central_downlink_Q_diana, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q_diana=sum(P_central_uplink_Q_diana);
cost_star_topology_downlink_diana=max(P_central_downlink);
cost_star_topology_Q_diana=cost_star_topology_uplink_Q_diana+cost_star_topology_downlink_diana;

%load diana_results.mat iter_diana iter_GD_2 Iter_gadmm Iter_2GADMM
load diana_results_2bits_rho24.mat iter_diana iter_GD_2 Iter_gadmm Iter_2GADMM

for iter=1:iter_diana
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_diana(iter)=cost_star_topology_Q_diana*transmissionTime;
    else
        cumulative_EnergyCost_diana(iter)= cumulative_EnergyCost_diana(iter-1)+cost_star_topology_Q_diana*transmissionTime;
    end        
end

for iter=1:iter_GD_2
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgd(iter)=cost_star_topology_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgd(iter)= cumulative_EnergyCost_qgd(iter-1)+cost_star_topology_Q*transmissionTime;
    end        
end

for iter=1:Iter_gadmm
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_gadmm(iter)=cost_decentralized*transmissionTime;
    else
        cumulative_EnergyCost_gadmm(iter)= cumulative_EnergyCost_gadmm(iter-1)+cost_decentralized*transmissionTime;
    end        
end

for iter=1:Iter_2GADMM
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgadmm(iter)=cost_decentralized_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgadmm(iter)= cumulative_EnergyCost_qgadmm(iter-1)+cost_decentralized_Q*transmissionTime;
    end        
end

final_ebnergyCost(rr,1)=cumulative_EnergyCost_qgd(iter_GD_2);
final_ebnergyCost(rr,2)=cumulative_EnergyCost_diana(iter_diana);
final_ebnergyCost(rr,3)=cumulative_EnergyCost_gadmm(Iter_gadmm);
final_ebnergyCost(rr,4)=cumulative_EnergyCost_qgadmm(Iter_2GADMM);

clearvars -except rr final_ebnergyCost transmissionTime

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(4);
subplot(1,3,1)
%semilogy(cumulative_com_bits_GD, loss_GD,'b-','LineWidth',3);
%hold on
cdfplot(final_ebnergyCost(:,1));
hold on
cdfplot(final_ebnergyCost(:,2));
hold on
cdfplot(final_ebnergyCost(:,3));
hold on
cdfplot(final_ebnergyCost(:,4));
xlabel({'Sum energy';'(a)'},'fontsize',16,'fontname','Times New Roman')
ylabel('CDF(sum energy)','fontsize',16,'fontname','Times New Roman')
legend('QGD-3bits','ADIANA','GADMM', 'Q-GADMM-3bits');
set(gca,'fontsize',14,'fontweight','bold');
%ylim([10^-4 10^3])
%xlim([10 6000])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear
%close all


for rr=1:100
[Xdata_28] = load('data99/data.txt'); 
[ydata_28] = load('data99/y.txt'); 
transmissionTime=1E-3;
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
no_workers = 50;
dataSamples_per_worker=floor(total_sample/no_workers);


total_sample =no_workers*dataSamples_per_worker;

num_iter=200000;
%tot_virtual_workers=24;
X=cell(no_workers);
y=cell(no_workers);

diagmatrix=diag(ones(total_sample,1));
% [lambda]=eig(Xdata'*Xdata);
%Hmax=zeros(num_workers,1);

for n=1:no_workers
        first = (n-1)*dataSamples_per_worker+1;
        last = first+dataSamples_per_worker-1;
        X{n}=Xdata_28(first:last,1:num_feature);
        y{n}=ydata_28(first:last);
        %Hmax(n)=max(eig(X{n}'*X{n}));
        
end


num_feature=size(X{1},2);







num_workers=no_workers;

bitsToSend=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%[path0, pathCost0, grid]=findPath2(num_workers);
Rate=32*num_feature/transmissionTime;%3.84E6;
sysBand=2E6;%20E6;
Band=sysBand/(num_workers/2);
[path, pathCost, d_square, P_central_uplink, P_central_downlink, center]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink=sum(P_central_uplink);
cost_star_topology_downlink=max(P_central_downlink);
cost_star_topology=cost_star_topology_uplink+cost_star_topology_downlink;
cost_decentralized=sum(pathCost);


Rate=(32+num_feature*3)/transmissionTime;%1E6;
%Band=1E6;
[~, pathCost_Q, ~, P_central_uplink_Q, P_central_downlink_Q, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q=sum(P_central_uplink_Q);
cost_star_topology_downlink_Q=max(P_central_downlink);
cost_star_topology_Q=cost_star_topology_uplink_Q+cost_star_topology_downlink_Q;
cost_decentralized_Q=sum(pathCost_Q);

Rate=(32+2*num_feature*3)/transmissionTime;%1.36E6;
%Band=1E6;
[~, pathCost_Q_diana, ~, P_central_uplink_Q_diana, P_central_downlink_Q_diana, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q_diana=sum(P_central_uplink_Q_diana);
cost_star_topology_downlink_diana=max(P_central_downlink);
cost_star_topology_Q_diana=cost_star_topology_uplink_Q_diana+cost_star_topology_downlink_diana;

%load diana_results.mat iter_diana iter_GD_2 Iter_gadmm Iter_2GADMM
load diana_results_2bits_rho24.mat iter_diana iter_GD_2 Iter_gadmm Iter_2GADMM

for iter=1:iter_diana
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_diana(iter)=cost_star_topology_Q_diana*transmissionTime;
    else
        cumulative_EnergyCost_diana(iter)= cumulative_EnergyCost_diana(iter-1)+cost_star_topology_Q_diana*transmissionTime;
    end        
end

for iter=1:iter_GD_2
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgd(iter)=cost_star_topology_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgd(iter)= cumulative_EnergyCost_qgd(iter-1)+cost_star_topology_Q*transmissionTime;
    end        
end

for iter=1:Iter_gadmm
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_gadmm(iter)=cost_decentralized*transmissionTime;
    else
        cumulative_EnergyCost_gadmm(iter)= cumulative_EnergyCost_gadmm(iter-1)+cost_decentralized*transmissionTime;
    end        
end

for iter=1:Iter_2GADMM
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgadmm(iter)=cost_decentralized_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgadmm(iter)= cumulative_EnergyCost_qgadmm(iter-1)+cost_decentralized_Q*transmissionTime;
    end        
end

final_ebnergyCost(rr,1)=cumulative_EnergyCost_qgd(iter_GD_2);
final_ebnergyCost(rr,2)=cumulative_EnergyCost_diana(iter_diana);
final_ebnergyCost(rr,3)=cumulative_EnergyCost_gadmm(Iter_gadmm);
final_ebnergyCost(rr,4)=cumulative_EnergyCost_qgadmm(Iter_2GADMM);

clearvars -except rr final_ebnergyCost transmissionTime

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(4);
subplot(1,3,2)
%semilogy(cumulative_com_bits_GD, loss_GD,'b-','LineWidth',3);
%hold on
cdfplot(final_ebnergyCost(:,1));
hold on
cdfplot(final_ebnergyCost(:,2));
hold on
cdfplot(final_ebnergyCost(:,3));
hold on
cdfplot(final_ebnergyCost(:,4));
xlabel({'Sum energy';'(b)'},'fontsize',16,'fontname','Times New Roman')
ylabel('CDF(sum energy)','fontsize',16,'fontname','Times New Roman')
legend('QGD-3bits','ADIANA','GADMM', 'Q-GADMM-3bits');
set(gca,'fontsize',14,'fontweight','bold');
%ylim([10^-4 10^3])
%xlim([10 6000])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




clear
%close all


for rr=1:100
[Xdata_28] = load('data99/data.txt'); 
[ydata_28] = load('data99/y.txt'); 
transmissionTime=1E-3;
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
no_workers = 50;
dataSamples_per_worker=floor(total_sample/no_workers);


total_sample =no_workers*dataSamples_per_worker;

num_iter=200000;
%tot_virtual_workers=24;
X=cell(no_workers);
y=cell(no_workers);

diagmatrix=diag(ones(total_sample,1));
% [lambda]=eig(Xdata'*Xdata);
%Hmax=zeros(num_workers,1);

for n=1:no_workers
        first = (n-1)*dataSamples_per_worker+1;
        last = first+dataSamples_per_worker-1;
        X{n}=Xdata_28(first:last,1:num_feature);
        y{n}=ydata_28(first:last);
        %Hmax(n)=max(eig(X{n}'*X{n}));
        
end


num_feature=size(X{1},2);







num_workers=no_workers;

bitsToSend=2;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%[path0, pathCost0, grid]=findPath2(num_workers);
Rate=32*num_feature/transmissionTime;%3.84E6;
sysBand=1E6;%20E6;
Band=sysBand/(num_workers/2);
[path, pathCost, d_square, P_central_uplink, P_central_downlink, center]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink=sum(P_central_uplink);
cost_star_topology_downlink=max(P_central_downlink);
cost_star_topology=cost_star_topology_uplink+cost_star_topology_downlink;
cost_decentralized=sum(pathCost);


Rate=(32+num_feature*3)/transmissionTime;%1E6;
%Band=1E6;
[~, pathCost_Q, ~, P_central_uplink_Q, P_central_downlink_Q, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q=sum(P_central_uplink_Q);
cost_star_topology_downlink_Q=max(P_central_downlink);
cost_star_topology_Q=cost_star_topology_uplink_Q+cost_star_topology_downlink_Q;
cost_decentralized_Q=sum(pathCost_Q);

Rate=(32+2*num_feature*3)/transmissionTime;%1.36E6;
%Band=1E6;
[~, pathCost_Q_diana, ~, P_central_uplink_Q_diana, P_central_downlink_Q_diana, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q_diana=sum(P_central_uplink_Q_diana);
cost_star_topology_downlink_diana=max(P_central_downlink);
cost_star_topology_Q_diana=cost_star_topology_uplink_Q_diana+cost_star_topology_downlink_diana;

%load diana_results.mat iter_diana iter_GD_2 Iter_gadmm Iter_2GADMM
load diana_results_2bits_rho24.mat iter_diana iter_GD_2 Iter_gadmm Iter_2GADMM

for iter=1:iter_diana
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_diana(iter)=cost_star_topology_Q_diana*transmissionTime;
    else
        cumulative_EnergyCost_diana(iter)= cumulative_EnergyCost_diana(iter-1)+cost_star_topology_Q_diana*transmissionTime;
    end        
end

for iter=1:iter_GD_2
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgd(iter)=cost_star_topology_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgd(iter)= cumulative_EnergyCost_qgd(iter-1)+cost_star_topology_Q*transmissionTime;
    end        
end

for iter=1:Iter_gadmm
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_gadmm(iter)=cost_decentralized*transmissionTime;
    else
        cumulative_EnergyCost_gadmm(iter)= cumulative_EnergyCost_gadmm(iter-1)+cost_decentralized*transmissionTime;
    end        
end

for iter=1:Iter_2GADMM
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgadmm(iter)=cost_decentralized_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgadmm(iter)= cumulative_EnergyCost_qgadmm(iter-1)+cost_decentralized_Q*transmissionTime;
    end        
end

final_ebnergyCost(rr,1)=cumulative_EnergyCost_qgd(iter_GD_2);
final_ebnergyCost(rr,2)=cumulative_EnergyCost_diana(iter_diana);
final_ebnergyCost(rr,3)=cumulative_EnergyCost_gadmm(Iter_gadmm);
final_ebnergyCost(rr,4)=cumulative_EnergyCost_qgadmm(Iter_2GADMM);

clearvars -except rr final_ebnergyCost transmissionTime

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(4);
subplot(1,3,3)
%semilogy(cumulative_com_bits_GD, loss_GD,'b-','LineWidth',3);
%hold on
cdfplot(final_ebnergyCost(:,1));
hold on
cdfplot(final_ebnergyCost(:,2));
hold on
cdfplot(final_ebnergyCost(:,3));
hold on
cdfplot(final_ebnergyCost(:,4));
xlabel({'Sum energy';'(c)'},'fontsize',16,'fontname','Times New Roman')
ylabel('CDF(sum energy)','fontsize',16,'fontname','Times New Roman')
legend('QGD-3bits','ADIANA','GADMM', 'Q-GADMM-3bits');
set(gca,'fontsize',14,'fontweight','bold');
%ylim([10^-4 10^3])
%xlim([10 6000])





