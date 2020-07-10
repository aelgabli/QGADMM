function [quantized,number_of_bits_toSend]=ADIANA_Q(current,prev,bitsToSend)

b=bitsToSend;
tau=1/(2^b-1);
if(b==1)
    
    R=1E-12;
for i=1:length(current)
    temp=abs(prev(i)-current(i));
    if(temp > R)
        R=temp;
    end
end
R=R/2;
number_of_bits_toSend =32+32*length(current);
for i=1:length(current)
    Q(i)=(current(i)-prev(i)+R)/(2*tau*R);
    Q(i)=round(Q(i));
%     p=(ceil(Q(i))-Q(i));%/(ceil(Q(i))-floor(Q(i)));
%     temp=rand;
%     if(temp <=p)
%         Q(i)=ceil(Q(i));
%     else
%         Q(i)=floor(Q(i));
%     end
    quantized(i)=2*tau*Q(i)*R-R;
    number_of_bits_toSend = number_of_bits_toSend + b;
    
end
    
else


R=1E-12;
for i=1:length(current)
    temp=abs(prev(i)-current(i));
    if(temp > R)
        R=temp;
    end
end

number_of_bits_toSend =32+32*length(current);
for i=1:length(current)
    Q(i)=(current(i)-prev(i)+R)/(2*tau*R);
    %p=(ceil(Q(i))-Q(i));%/(ceil(Q(i))-floor(Q(i)));
    Q(i)=round(Q(i));
    %temp=rand;
%     if(temp <=p)
%         Q(i)=ceil(Q(i));
%     else
%         Q(i)=floor(Q(i));
%     end
    quantized(i)=2*tau*Q(i)*R-R;
    number_of_bits_toSend = number_of_bits_toSend + b;
    
end
end    
