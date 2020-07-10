function [obj_GADMM, loss_GADMM, Iter, feas, comp_time]=group_ADMM_closedForm_2...
    (XX,YY, rho, no_workers, num_feature, noSamples, num_iter, obj0, acc, alpha)
Iter= num_iter;   
           
s1=num_feature;
s2=noSamples;
lambda = zeros(s1,no_workers);
out=zeros(s1,no_workers);
comp_time(1)=0;
max_iter = num_iter;
 for i = 1:max_iter
     
     if (i>1)
         comp_time(i)=comp_time(i-1);
     end
     
     for ii =1:2:no_workers
         if ii==1
             C1 = zeros(s1,1);
             term_1 =0;
         else
             C1= lambda(:,ii-1);
             term_1=rho*out(:,ii-1);
         end
         
         if ii == no_workers
             C2 = zeros(s1,1);
             term_2 = 0;
         else
             C2= lambda(:,ii);
             term_2 = rho*out(:,ii+1);
             
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
            if(ii==1)
                 comp_time(i)=comp_time(i)+toc;
             end
        else
            x=(H'*H+2*rho*eye(s1,s1))\(H'*Y+C1-C2+term_1+term_2);
        end
        
        out(:,ii) =x;
        
     end
    
    
        
     for ii =2:2:no_workers
         C1= lambda(:,ii-1);
         if ii == no_workers
             C2 = zeros(s1,1);
             term_2 = 0;
         else
             C2= lambda(:,ii);
             term_2 = rho*out(:,ii+1);
         end
         term_1=rho*out(:,ii-1);
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
     end
        for ii=1:no_workers-1
            lambda(:,ii) = lambda(:,ii) + rho*alpha*(out(:,ii)-out(:,ii+1));
        end
        final_obj = 0;
        for ii =1:no_workers
            first = (ii-1)*s2+1;
            last = first+s2-1;
            final_obj = final_obj + 0.5*sum_square(XX(first:last,1:s1)*out(:,ii) - YY(first:last));
        end
        %i
        obj_GADMM(i)=final_obj;
        loss_GADMM(i)=abs(final_obj-obj0);%*100;
        i;
        toPrintOut_gadmm=loss_GADMM(i);
        %loss_GADMM(i)=abs(final_obj-obj0);
            temp=0;
         for ii =1:no_workers-1
            temp=temp+sum(abs(out(:,ii)-out(:,ii+1)));
         end 
            feas(i)=temp;
            
        if(i > 100 && loss_GADMM(i) < acc)
            Iter = i;
            break;
        end
       
    

end
     




