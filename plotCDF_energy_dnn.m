clear
close all

for rr=1:100
transmissionTime=100E-3;

no_workers = 10;

target_acc=0.9;

s1=109184;

num_workers=no_workers;

bitsToSend=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%[
Rate=34938880;%needs to transmit 32*6 bits in 100 microseconds
sysBand=400E6;
Band=sysBand/(num_workers/2);
%Band=1E6;
[path, pathCost, d_square, P_central_uplink, P_central_downlink, center]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink=sum(P_central_uplink);
cost_star_topology_downlink=max(P_central_downlink);
cost_star_topology=cost_star_topology_uplink+cost_star_topology_downlink;
cost_decentralized=sum(pathCost);

%quantized algorithms
Rate=8735040;%needs to transmit 32+6*3 bits in 100 microseconds
%Band=1E6;
[~, pathCost_Q, ~, P_central_uplink_Q, P_central_downlink_Q, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q=sum(P_central_uplink_Q);
cost_star_topology_downlink_Q=max(P_central_downlink_Q);
cost_star_topology_Q=cost_star_topology_uplink_Q+cost_star_topology_downlink_Q;
cost_decentralized_Q=sum(pathCost_Q);

% %A-DIANA 
% Rate=1.36E6;%needs to transmit 32+2*6*3 bits in 100 microseconds
% %Band=1E6;
% [~, pathCost_Q_diana, ~, P_central_uplink_Q_diana, P_central_downlink_Q_diana, ~]=findPath2(num_workers, Rate, Band);
% 
% cost_star_topology_uplink_Q_diana=sum(P_central_uplink_Q_diana);
% cost_star_topology_downlink_diana=max(P_central_downlink);
% cost_star_topology_Q_diana=cost_star_topology_uplink_Q_diana+cost_star_topology_downlink_diana;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%


acc_GD = readNPY('acc_test_gd.npy');
acc_GD_2 = readNPY('acc_test_qgd.npy');

acc_GADMM = readNPY('gadmm_dnn_max_test.npy');
acc_2GADMM = readNPY('qgadmm_dnn_max_test.npy');



maxIter=399;
for iter=1:maxIter
    cumulative_com_rounds_GD_2(iter)=iter*num_workers+iter; 
    cumulative_com_bits_GD_2(iter)=iter*num_workers*(32+s1*bitsToSend)+iter*32*s1;


end

%num_iter=iter_GD_2;

for iter=1:maxIter
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgd(iter)=cost_star_topology_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgd(iter)= cumulative_EnergyCost_qgd(iter-1)+cost_star_topology_Q*transmissionTime;
    end
    
    if(acc_GD_2 > target_acc)
        break;
    end
end

cumulative_EnergyCost_final(rr,1)= cumulative_EnergyCost_qgd(iter);
    
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
    
    if(acc_GADMM > target_acc)
        break;
    end
    
end
    
cumulative_EnergyCost_final(rr,2)=cumulative_EnergyCost_gadmm(iter);
    
    
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
    if(acc_2GADMM > target_acc)
        break;
    end
end

cumulative_EnergyCost_final(rr,3)=cumulative_EnergyCost_qgadmm(iter);

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

rr

clearvars -except rr cumulative_EnergyCost_final cumulative_EnergyCost_final...
    cumulative_EnergyCost_final

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(5);
subplot(1,3,1);
cdfplot(cumulative_EnergyCost_final(:,1));
hold on
cdfplot(cumulative_EnergyCost_final(:,2));
hold on
cdfplot(cumulative_EnergyCost_final(:,3));
xlabel({'Sum energy';'(a)'},'fontsize',16,'fontname','Times New Roman')
ylabel('CDF(sum energy)','fontsize',16,'fontname','Times New Roman')
legend('SQGD-8bits', 'SGADMM', 'Q-SGADMM-8bits');
set(gca,'fontsize',14,'fontweight','bold');
%xlim([0 1E12])




clear
%close all

for rr=1:100
transmissionTime=100E-3;

no_workers = 10;

target_acc=0.9;

s1=109184;

num_workers=no_workers;

bitsToSend=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%[
Rate=34938880;%needs to transmit 32*6 bits in 100 microseconds
sysBand=100E6;
Band=sysBand/(num_workers/2);
%Band=1E6;
[path, pathCost, d_square, P_central_uplink, P_central_downlink, center]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink=sum(P_central_uplink);
cost_star_topology_downlink=max(P_central_downlink);
cost_star_topology=cost_star_topology_uplink+cost_star_topology_downlink;
cost_decentralized=sum(pathCost);

%quantized algorithms
Rate=8735040;%needs to transmit 32+6*3 bits in 100 microseconds
%Band=1E6;
[~, pathCost_Q, ~, P_central_uplink_Q, P_central_downlink_Q, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q=sum(P_central_uplink_Q);
cost_star_topology_downlink_Q=max(P_central_downlink_Q);
cost_star_topology_Q=cost_star_topology_uplink_Q+cost_star_topology_downlink_Q;
cost_decentralized_Q=sum(pathCost_Q);

% %A-DIANA 
% Rate=1.36E6;%needs to transmit 32+2*6*3 bits in 100 microseconds
% %Band=1E6;
% [~, pathCost_Q_diana, ~, P_central_uplink_Q_diana, P_central_downlink_Q_diana, ~]=findPath2(num_workers, Rate, Band);
% 
% cost_star_topology_uplink_Q_diana=sum(P_central_uplink_Q_diana);
% cost_star_topology_downlink_diana=max(P_central_downlink);
% cost_star_topology_Q_diana=cost_star_topology_uplink_Q_diana+cost_star_topology_downlink_diana;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%


acc_GD = readNPY('acc_test_gd.npy');
acc_GD_2 = readNPY('acc_test_qgd.npy');

acc_GADMM = readNPY('gadmm_dnn_max_test.npy');
acc_2GADMM = readNPY('qgadmm_dnn_max_test.npy');



maxIter=399;
for iter=1:maxIter
    cumulative_com_rounds_GD_2(iter)=iter*num_workers+iter; 
    cumulative_com_bits_GD_2(iter)=iter*num_workers*(32+s1*bitsToSend)+iter*32*s1;


end

%num_iter=iter_GD_2;

for iter=1:maxIter
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgd(iter)=cost_star_topology_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgd(iter)= cumulative_EnergyCost_qgd(iter-1)+cost_star_topology_Q*transmissionTime;
    end
    
    if(acc_GD_2 > target_acc)
        break;
    end
end

cumulative_EnergyCost_final(rr,1)= cumulative_EnergyCost_qgd(iter);
    
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
    
    if(acc_GADMM > target_acc)
        break;
    end
    
end
    
cumulative_EnergyCost_final(rr,2)=cumulative_EnergyCost_gadmm(iter);
    
    
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
    if(acc_2GADMM > target_acc)
        break;
    end
end

cumulative_EnergyCost_final(rr,3)=cumulative_EnergyCost_qgadmm(iter);

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

rr

clearvars -except rr cumulative_EnergyCost_final cumulative_EnergyCost_final...
    cumulative_EnergyCost_final

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(5);
subplot(1,3,2);
cdfplot(cumulative_EnergyCost_final(:,1));
%hold on
%cdfplot(cumulative_EnergyCost_final(:,2));
hold on
cdfplot(cumulative_EnergyCost_final(:,3));
xlabel({'Sum energy';'(b)'},'fontsize',16,'fontname','Times New Roman')
ylabel('CDF(sum energy)','fontsize',16,'fontname','Times New Roman')
legend('SQGD-8bits', 'Q-SGADMM-8bits');
set(gca,'fontsize',14,'fontweight','bold');
%xlim([0 1E12])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
%close all

for rr=1:100
transmissionTime=100E-3;

no_workers = 10;

target_acc=0.9;

s1=109184;

num_workers=no_workers;

bitsToSend=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%

%[
Rate=34938880;%needs to transmit 32*6 bits in 100 microseconds
sysBand=40E6;
Band=sysBand/(num_workers/2);
%Band=1E6;
[path, pathCost, d_square, P_central_uplink, P_central_downlink, center]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink=sum(P_central_uplink);
cost_star_topology_downlink=max(P_central_downlink);
cost_star_topology=cost_star_topology_uplink+cost_star_topology_downlink;
cost_decentralized=sum(pathCost);

%quantized algorithms
Rate=8735040;%needs to transmit 32+6*3 bits in 100 microseconds
%Band=1E6;
[~, pathCost_Q, ~, P_central_uplink_Q, P_central_downlink_Q, ~]=findPath2(num_workers, Rate, Band);

cost_star_topology_uplink_Q=sum(P_central_uplink_Q);
cost_star_topology_downlink_Q=max(P_central_downlink_Q);
cost_star_topology_Q=cost_star_topology_uplink_Q+cost_star_topology_downlink_Q;
cost_decentralized_Q=sum(pathCost_Q);

% %A-DIANA 
% Rate=1.36E6;%needs to transmit 32+2*6*3 bits in 100 microseconds
% %Band=1E6;
% [~, pathCost_Q_diana, ~, P_central_uplink_Q_diana, P_central_downlink_Q_diana, ~]=findPath2(num_workers, Rate, Band);
% 
% cost_star_topology_uplink_Q_diana=sum(P_central_uplink_Q_diana);
% cost_star_topology_downlink_diana=max(P_central_downlink);
% cost_star_topology_Q_diana=cost_star_topology_uplink_Q_diana+cost_star_topology_downlink_diana;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%


acc_GD = readNPY('acc_test_gd.npy');
acc_GD_2 = readNPY('acc_test_qgd.npy');

acc_GADMM = readNPY('gadmm_dnn_max_test.npy');
acc_2GADMM = readNPY('qgadmm_dnn_max_test.npy');



maxIter=399;
for iter=1:maxIter
    cumulative_com_rounds_GD_2(iter)=iter*num_workers+iter; 
    cumulative_com_bits_GD_2(iter)=iter*num_workers*(32+s1*bitsToSend)+iter*32*s1;


end

%num_iter=iter_GD_2;

for iter=1:maxIter
    
    %cumulative_com_ADMM_rho5(iter)=iter*num_workers;   
    if(iter==1)
        cumulative_EnergyCost_qgd(iter)=cost_star_topology_Q*transmissionTime;
    else
        cumulative_EnergyCost_qgd(iter)= cumulative_EnergyCost_qgd(iter-1)+cost_star_topology_Q*transmissionTime;
    end
    
    if(acc_GD_2 > target_acc)
        break;
    end
end

cumulative_EnergyCost_final(rr,1)= cumulative_EnergyCost_qgd(iter);
    
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
    
    if(acc_GADMM > target_acc)
        break;
    end
    
end
    
cumulative_EnergyCost_final(rr,2)=cumulative_EnergyCost_gadmm(iter);
    
    
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
    if(acc_2GADMM > target_acc)
        break;
    end
end

cumulative_EnergyCost_final(rr,3)=cumulative_EnergyCost_qgadmm(iter);

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

rr

clearvars -except rr cumulative_EnergyCost_final cumulative_EnergyCost_final...
    cumulative_EnergyCost_final

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(5);
subplot(1,3,3);
cdfplot(cumulative_EnergyCost_final(:,1));
%hold on
%cdfplot(cumulative_EnergyCost_final(:,2));
hold on
cdfplot(cumulative_EnergyCost_final(:,3));
xlabel({'Sum energy';'(c)'},'fontsize',16,'fontname','Times New Roman')
ylabel('CDF(sum energy)','fontsize',16,'fontname','Times New Roman')
legend('SQGD-8bits', 'Q-SGADMM-8bits');
set(gca,'fontsize',14,'fontweight','bold');
%xlim([0 1E12])







