function [quantized,number_of_bits_toSend]=gadmm_stochasticQ(quantized,current,prev,bitsToSend)


b=bitsToSend;
tau=1/(2^b-1);
if(b<3)
R=1E-12;
for i=1:length(current)
    temp=abs(prev(i)-current(i));
    if(temp > R)
        R=temp;
    end
end

% Quantization
number_of_bits_toSend =32;% the number of bits to send the value of R.
%b=2; % for level quantization
for i=1:length(current)
    Q(i)=(current(i)-prev(i)+R)/(2*tau*R);
    Q(i)=round(Q(i));
    quantized(i)=quantized(i)+2*tau*Q(i)*R-R;
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
% Quantization
number_of_bits_toSend =32;% the number of bits to send the value of R.
%b=2; % for level quantization
for i=1:length(current)
    Q(i)=(current(i)-prev(i)+R)/(2*tau*R);
    p=(ceil(Q(i))-Q(i));%/(ceil(Q(i))-floor(Q(i)));
    temp=rand;
    if(temp <=p)
        Q(i)=ceil(Q(i));
    else
        Q(i)=floor(Q(i));
    end
    quantized(i)=quantized(i)+2*tau*Q(i)*R-R;
    number_of_bits_toSend = number_of_bits_toSend + b;
    
end
end
end    

