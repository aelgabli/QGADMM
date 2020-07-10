clear
close all


%EbNo = 10;

transmissionTime=100E-3;

no_workers = 10;


s1=109184;

num_workers=no_workers;

bitsToSend=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%[
Rate=34938880;%need to transmit 109184*32 bits in 100ms
sysBand=40E6;
Band=sysBand/(num_workers/2);
%Band=1E6;
[path, pathCost, d_square, P_central_uplink, P_central_downlink, center]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink=sum(P_central_uplink);
cost_star_topology_downlink=max(P_central_downlink);
cost_star_topology=cost_star_topology_uplink+cost_star_topology_downlink;
cost_decentralized=sum(pathCost);

%quantized algorithms
Rate=8735040;%need to transmit 32+109184*8 bits in 100ms
%Band=1E6;
[~, pathCost_Q, ~, P_central_uplink_Q, P_central_downlink_Q, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q=sum(P_central_uplink_Q);
cost_star_topology_downlink_Q=max(P_central_downlink_Q);
cost_star_topology_Q=cost_star_topology_uplink_Q+cost_star_topology_downlink_Q;
cost_decentralized_Q=sum(pathCost_Q);




acc_GD = readNPY('acc_test_gd.npy');
acc_GD_2 = readNPY('acc_test_qgd.npy');

acc_GADMM = readNPY('gadmm_dnn_max_test.npy');
acc_2GADMM = readNPY('qgadmm_dnn_max_test.npy');



%%  QGD

maxIter=399;
for iter=1:maxIter
    cumulative_com_rounds_GD_2(iter)=iter*num_workers+iter; 
    cumulative_com_bits_GD_2(iter)=iter*num_workers*(32+s1*bitsToSend)+iter*32*s1; 
    %errorPer_GD(iter) = abs(loss_GD(iter)/opt_obj(iter)*100);
end

%num_iter=iter_GD_2;

for iter=1:maxIter
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgd(iter)=cost_star_topology_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgd(iter)= cumulative_EnergyCost_qgd(iter-1)+cost_star_topology_Q*transmissionTime;
    end        
end


    
    for iter=1:maxIter%num_iter%Iter_gadmm
        cumulative_com_GADMM(iter)=iter*no_workers*s1*32;
    end
    
    for iter=1:maxIter%num_iter%Iter_gadmm
        cumulative_com_rounds_GADMM(iter)=iter*no_workers;
    end
    


for iter=1:maxIter
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_gadmm(iter)=cost_decentralized*transmissionTime;
    else
        cumulative_EnergyCost_gadmm(iter)= cumulative_EnergyCost_gadmm(iter-1)+cost_decentralized*transmissionTime;
    end        
end
    

    
    
    for iter=1:maxIter%num_iter%Iter_gadmm
        cumulative_com_2GADMM(iter)=iter*num_workers*(32+s1*bitsToSend);%iter*no_workers*(32+num_feature*);
    end
    
    for iter=1:maxIter%num_iter%Iter_gadmm
        cumulative_com_rounds_2GADMM(iter)=iter*no_workers;
    end
    

for iter=1:maxIter
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgadmm(iter)=cost_decentralized_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgadmm(iter)= cumulative_EnergyCost_qgadmm(iter-1)+cost_decentralized_Q*transmissionTime;
    end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GD
for iter=1:maxIter
    cumulative_com_rounds_GD(iter)=iter*num_workers+iter; 
    cumulative_com_bits_GD(iter)=iter*num_workers*s1*32+iter*32*s1; 
    %errorPer_GD(iter) = abs(loss_GD(iter)/opt_obj(iter)*100);
end

for iter=1:maxIter
    
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
plot(cumulative_com_rounds_GD, acc_GD(1:399),'b-','LineWidth',3);
hold on
plot(cumulative_com_rounds_GD_2, acc_GD_2(1:399),'m--','LineWidth',3);
hold on
%semilogy(cumulative_com_rounds_diana, acc_diana,'c--','LineWidth',3);

plot(cumulative_com_rounds_GADMM, acc_GADMM(1:399),'k-','LineWidth',3);
hold on
plot(cumulative_com_rounds_2GADMM, acc_2GADMM(1:399),'r--','LineWidth',3);
hold on
%semilogy(cumulative_com_rounds_LAG, loss_LAG_WK,'b:','LineWidth',3);
%semilogy(loss_1GADMM,'m-','LineWidth',3);
%semilogy(loss_17GADMM,'k--','LineWidth',3);
xlabel({'Number of communication rounds';'(a)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Test accuracy','fontsize',16,'fontname','Times New Roman')
legend('SGD','QSGD-8bits','GADMM', 'QSGADMM-8bits');%,'ADIANA'
%ylim([9E-2 10^0])
%xlim([20 1800])
%xlim([50 6000])
xlim([0 3500])

set(gca,'fontsize',14,'fontweight','bold');

subplot(1,3,2);
plot(cumulative_com_bits_GD, acc_GD(1:399),'b-','LineWidth',3);
hold on
plot(cumulative_com_bits_GD_2, acc_GD_2(1:399),'m--','LineWidth',3);
hold on
%semilogy(cumulative_com_bits_diana, loss_diana,'c--','LineWidth',3);

hold on
plot(cumulative_com_GADMM, acc_GADMM(1:399),'k-','LineWidth',3);
hold on
plot(cumulative_com_2GADMM, acc_2GADMM(1:399),'r--','LineWidth',3);
%semilogy(cumulative_com_bits_lag, loss_LAG_WK,'b:','LineWidth',3);
%semilogy(cumulative_com_1GADMM,'m-','LineWidth',3);
%semilogy(cumulative_com_17GADMM,'k--','LineWidth',3);
xlabel({'Total number of transmitted bits';'(b)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Test accuracy','fontsize',16,'fontname','Times New Roman')
legend('SGD','QSGD-8bits','GADMM', 'QSGADMM-8bits');%,'ADIANA'
set(gca,'fontsize',14,'fontweight','bold');
%ylim([9E-2 10^0])
xlim([0 2.5E9])



subplot(1,3,3);
%semilogy(cumulative_EnergyCost_gd, loss_GD,'b-','LineWidth',3);
%hold on
plot(cumulative_EnergyCost_qgd, acc_GD_2(1:399),'m--','LineWidth',3);
hold on
%semilogy(cumulative_EnergyCost_diana, loss_diana,'c--','LineWidth',3);

hold on
plot(cumulative_EnergyCost_gadmm, acc_GADMM(1:399),'k-','LineWidth',3);
hold on
plot(cumulative_EnergyCost_qgadmm, acc_2GADMM(1:399),'r--','LineWidth',3);
xlabel({'Sum energy';'(c)'},'fontsize',16,'fontname','Times New Roman')
ylabel('Test accuarcy','fontsize',16,'fontname','Times New Roman')
legend('QSGD-8bits','SGADMM', 'QSGADMM-8bits');%,'ADIANA'
set(gca,'fontsize',14,'fontweight','bold');
%ylim([10^-4 10^3])
xlim([0 7E7])

