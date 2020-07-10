clear
close all


%EbNo = 10;
%for rr=1:100
transmissionTime=100E-3;

no_workers = 10;


s1=109184;

bitsToSend=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%

gadmm_10=readNPY('gadmm_dnn_max_test_rho10.npy');
gadmm_20=readNPY('gadmm_dnn_max_test_rho20.npy');
gadmm_30=readNPY('gadmm_dnn_max_test_rho30.npy');

qgadmm_10=readNPY('qgadmm_dnn_max_test_rho10.npy');
qgadmm_20=readNPY('qgadmm_dnn_max_test_rho20.npy');
qgadmm_30=readNPY('qgadmm_dnn_max_test_rho30.npy');




maxIter=399;

    
    for iter=1:maxIter%num_iter%Iter_gadmm
        cumulative_com_GADMM(iter)=iter*no_workers*s1*32;
    end
    
        
    for iter=1:maxIter%num_iter%Iter_gadmm
        cumulative_com_2GADMM(iter)=iter*no_workers*(32+s1*bitsToSend);%iter*no_workers*(32+num_feature*);
    end
    

    



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
plot(cumulative_com_GADMM, gadmm_10(1:399),'LineWidth',3);
hold on
plot(cumulative_com_2GADMM, qgadmm_10(1:399),'LineWidth',3);

plot(cumulative_com_GADMM, gadmm_20(1:399),'LineWidth',3);
hold on
plot(cumulative_com_2GADMM, qgadmm_20(1:399),'LineWidth',3);

plot(cumulative_com_GADMM, gadmm_30(1:399),'LineWidth',3);
hold on
plot(cumulative_com_2GADMM, qgadmm_30(1:399),'LineWidth',3);

xlabel({'Total number of transmitted bits'},'fontsize',16,'fontname','Times New Roman')
ylabel('Test accuracy','fontsize',16,'fontname','Times New Roman')
legend('GADMM-\rho=10', 'QSGADMM-8bits-\rho=10','GADMM-\rho=20', 'QSGADMM-8bits-\rho=20','GADMM-\rho=30', 'QSGADMM-8bits-\rho=30');%,'ADIANA'
%ylim([9E-2 10^0])
xlim([0 3E9])






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





clear



noWorkers_set=[10, 20, 30];



for rr=1:3
    
    
transmissionTime=100E-3;
s1=109184;
bitsToSend=8;
%%%%%%%%%%%%%%%%%%%%%%%%%%

gadmm_10=readNPY('gadmm_dnn_max_test_users10.npy');
gadmm_20=readNPY('gadmm_dnn_max_test_users20.npy');
gadmm_30=readNPY('gadmm_dnn_max_test_users30.npy');

qgadmm_10=readNPY('qgadmm_dnn_max_test_rho20.npy');
qgadmm_20=readNPY('qgadmm_dnn_max_test_users20.npy');
qgadmm_30=readNPY('qgadmm_dnn_max_test_users30.npy');    

    
no_workers = noWorkers_set(rr);

maxIter=399;

    
    for iter=1:maxIter%num_iter%Iter_gadmm
        cumulative_com_GADMM(iter)=iter*no_workers*s1*32;
    end
    
        
    for iter=1:maxIter%num_iter%Iter_gadmm
        cumulative_com_2GADMM(iter)=iter*no_workers*(32+s1*bitsToSend);%iter*no_workers*(32+num_feature*);
    end
    

  cumulative_com_gadmm_final(rr,:)=cumulative_com_GADMM;
  cumulative_com_qgadmm_final(rr,:)=cumulative_com_2GADMM;
  
  clearvars -except rr gadmm_final cumulative_com_gadmm_final cumulative_com_qgadmm_final noWorkers_set
  
end

gadmm_10=readNPY('gadmm_dnn_max_test_users10.npy');
gadmm_20=readNPY('gadmm_dnn_max_test_users20.npy');
gadmm_30=readNPY('gadmm_dnn_max_test_users30.npy');

qgadmm_10=readNPY('qgadmm_dnn_max_test_rho20.npy');
qgadmm_20=readNPY('qgadmm_dnn_max_test_users20.npy');
qgadmm_30=readNPY('qgadmm_dnn_max_test_users30.npy'); 

figure(4);
plot(cumulative_com_gadmm_final(1,:), gadmm_10(1:399),'LineWidth',3);
hold on
plot(cumulative_com_qgadmm_final(1,:), qgadmm_10(1:399),'LineWidth',3);

plot(cumulative_com_gadmm_final(2,:), gadmm_20(1:399),'LineWidth',3);
hold on
plot(cumulative_com_qgadmm_final(2,:), qgadmm_20(1:399),'LineWidth',3);

plot(cumulative_com_gadmm_final(3,:), gadmm_30(1:399),'LineWidth',3);
hold on
plot(cumulative_com_qgadmm_final(3,:), qgadmm_30(1:399),'LineWidth',3);

xlabel({'Total number of transmitted bits'},'fontsize',16,'fontname','Times New Roman')
ylabel('Test accuracy','fontsize',16,'fontname','Times New Roman')
legend('GADMM-10 workers', 'QSGADMM-8bits-10 workers','GADMM-20 workers', 'QSGADMM-8bits-20 workers','GADMM-30 workers', 'QSGADMM-8bits-30 workers');%,'ADIANA'
%ylim([9E-2 10^0])
xlim([0 5E9])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
gadmm=readNPY('gadmm_dnn_max_test_rho20.npy');
qgadmm=readNPY('qgadmm_dnn_max_test_rho20.npy');
compTime_gadmm=readNPY('gadmm_dnn_compTime.npy');
compTime_qgadmm=readNPY('qgadmm_dnn_compTime.npy');

figure(5);
plot(compTime_gadmm, gadmm,'k-','LineWidth',3);
hold on
plot(compTime_qgadmm, qgadmm,'r--','LineWidth',3);
xlabel({'Local computation time';(b)},'fontsize',16,'fontname','Times New Roman')
ylabel('Test accuracy','fontsize',16,'fontname','Times New Roman')
legend('GADMM', 'QSGADMM-8bits');
%ylim([9E-2 10^0])
%xlim([0 5E9])



