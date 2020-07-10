function [obj_GADMM, loss_GADMM, Iter, cumulative_com_GADMM,feas, comp_time]=quantized_group_ADMM...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, acc,alpha,bitsToSend)

Iter= num_iter; 

%bitsToSend_w=

s1=num_feature;
s2=noSamples;


lambda_1 = zeros(s1,no_workers);
lambda_2 = zeros(s1,no_workers);

out=zeros(s1,no_workers);
prev_out=zeros(s1,no_workers);
q_out=zeros(s1,no_workers);
%q_out_r=zeros(s1,no_workers);
%delta=zeros(s1,no_workers);

max_iter = num_iter;
cumulative_com_GADMM(1)=0;
comp_time(1)=0;
 for i = 1:max_iter
    % tic

     if(i > 1)
        cumulative_com_GADMM(i)=cumulative_com_GADMM(i-1);
        comp_time(i)=comp_time(i-1);
     end

    
     for ii =1:2:no_workers
        
%          if(ismember(ii,oddrnd(i,:)))
%             cumulative_com_GADMM(i)=cumulative_com_GADMM(i)+2;
%          end
         
         
      %if(ismember(ii,oddrnd(i,:)))   
         %cvx_begin quiet 
     %cvx_precision high
     %cvx_solver  mosek%SDPT3 %mosek
     %variable x(s1)
     
         if ii==1
             C1 = zeros(s1,1);
             term_1 =0;
         else
             C1= lambda_1(:,ii);
             term_1=rho*q_out(:,ii-1);
         end
         
         if ii == no_workers
             C2 = zeros(s1,1);
             term_2 = 0;
         else
             C2= lambda_2(:,ii);
             term_2 = rho*q_out(:,ii+1);
             
         end
         first = (ii-1)*s2+1;
         last = first+s2-1;
        
         H=XX(first:last,1:s1);
         Y=YY(first:last);
         %temp1=(H'*H+2*rho*eye(s1,s1));
         %temp2=(H'*Y+C1-C2+term_1+term_2);
         %x=(H'*H+2*rho*eye(s1,s1))\(H'*Y+C1-C2+term_1+term_2);
         if(ii==1||ii==no_workers)
             if(ii==1)
                 tic
             end
            x=(H'*H+rho*eye(s1,s1))\(H'*Y+C1-C2+term_1+term_2);   
        else
            x=(H'*H+2*rho*eye(s1,s1))\(H'*Y+C1-C2+term_1+term_2);
        end
        %cvx_solver mosek
        
        %minimize(0.5*sum_square(XX(first:last,1:s1)*x - YY(first:last))- C1'*x+C2'*x + term_1...
            %+term_2)
         %cvx_end
        out(:,ii) =x;
        nq_out(:,ii)=x;
        %if(i~=1)
                                                                
            [q_out(:,ii),number_of_bits_toSend]=gadmm_stochasticQ_v2(q_out(:,ii),out(:,ii),prev_out(:,ii),bitsToSend);
            
            if(ii==1)
                 comp_time(i)=comp_time(i)+toc;
             end
            
            cumulative_com_GADMM(i)=cumulative_com_GADMM(i)+number_of_bits_toSend;
            %aa1=q_out_l(:,ii)

        %end
       out(:,ii) = q_out(:,ii);
       prev_out(:,ii) = q_out(:,ii); 
       
        
     end
    
    
        
     for ii =2:2:no_workers
       %if(ismember(ii+1,oddrnd(i,:)) && ismember(ii-1,oddrnd(i,:)))  
          %cvx_begin quiet 
        %cvx_precision high
        %cvx_solver mosek%SDPT3 %mosek
        %variable x(s1)
         C1= lambda_1(:,ii);
         if ii == no_workers
             C2 = zeros(s1,1);
             term_2 = 0;
         else
             C2= lambda_2(:,ii);
             term_2 = rho*q_out(:,ii+1);
         end
         term_1=rho*q_out(:,ii-1);
         first = (ii-1)*s2+1;
         last = first+s2-1;
         
        H=XX(first:last,1:s1);
        Y=YY(first:last);
        %x=(H'*H+2*rho*eye(s1,s1))\(H'*Y+C1-C2+term_1+term_2);
        if(ii==no_workers)
            x=(H'*H+rho*eye(s1,s1))\(H'*Y+C1-C2+term_1+term_2);   
        else
            x=(H'*H+2*rho*eye(s1,s1))\(H'*Y+C1-C2+term_1+term_2);
        end
        out(:,ii) =x;
        nq_out(:,ii)=x;
        %if(i~=1)
            
            [q_out(:,ii),number_of_bits_toSend]=gadmm_stochasticQ_v2(q_out(:,ii),out(:,ii),prev_out(:,ii),bitsToSend);
            cumulative_com_GADMM(i)=cumulative_com_GADMM(i)+number_of_bits_toSend;
        

        %end
        out(:,ii) = q_out(:,ii);
       prev_out(:,ii) = out(:,ii); 
     %end
     end
     
        for ii=1:no_workers
            %if(ismember(ii,oddrnd(i,:)) || ismember(ii+1,oddrnd(i,:)))
                if(ii==1)
                    lambda_1(:,ii) = zeros;
                else
                    lambda_1(:,ii) = lambda_1(:,ii) + rho*alpha*(q_out(:,ii-1)-out(:,ii));
                end
                if(ii==no_workers)
                    
                    lambda_2(:,ii) = zeros;
                else                        
                lambda_2(:,ii) = lambda_2(:,ii) + rho*alpha*(out(:,ii)-q_out(:,ii+1));
                end
%                 lambda_1(:,ii) = lambda_1(:,ii) + rho*alpha*(out(:,ii)-q_out_r(:,ii+1));
%                 lambda_2(:,ii) = lambda_2(:,ii) + rho*alpha*(q_out_l(:,ii)-out(:,ii+1));
            %end
        end
        
        
        final_obj = 0;
        for ii =1:no_workers
            first = (ii-1)*s2+1;
            last = first+s2-1;
            final_obj = final_obj + 0.5*sum_square(XX(first:last,1:s1)*out(:,ii) - YY(first:last));
        end
        obj_GADMM(i)=final_obj;
        loss_GADMM(i)=abs(final_obj-obj0);%*100;
        toPrintOut_qgadmm=loss_GADMM(i)
        
           temp=0;
         for ii =1:no_workers-1
            temp=temp+sum(abs(out(:,ii)-out(:,ii+1)));
         end 
            feas(i)=temp;
            

          if(loss_GADMM(i) < acc)
            Iter = i;
            break;
        end
       
    

 end

 end
     





    
                
    

