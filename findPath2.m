function [path, pathCost_final, d_square, P_central_uplink, P_central_downlink, center]=findPath2(N, R, B)
X=250;
Y=250;
eta=1E-6;

d_square=zeros(N,N);
%d_square_central=zeros(N,N);
for n=1:N
  x(n)=X*rand;
  y(n)=Y*rand;
end

%center=sqrt(25^2+25^2);

min_dist=1E15;
for n=1:N
    dist=sqrt((x(n)-X/2)^2+(y(n)-Y/2)^2);
    if (dist < min_dist)
        min_dist=dist;
        center=n;
    end
end

%P_central_downlink=0;
for n=1:N
    d_square_central(n)=(x(n)-x(center))^2+(y(n)-y(center))^2;
    P_central_uplink(n)= 1/2*d_square_central(n)*eta*B*(2^(2*R/B)-1);
    P_central_downlink(n)= N/2*B*d_square_central(n)*eta*(2^(2*R/(N*B))-1);
    %temp= d_square_central(n)*eta*B*2^(R/B);
%     if(temp > P_central_downlink)
%         P_central_downlink=temp;
%     end
end


        
for n=1:N
    for m=1:N
        if(n~=m)
            d_square(n,m)=(x(n)-x(m))^2+(y(n)-y(m))^2;
            P(n,m)= d_square(n,m)*eta*B*(2^(R/B)-1);
        end
    end
end

    
path=[1];
pathCost=[];
while(length(path) < N)
    n=path(length(path));
    min=1E12;
    for m=1:N
        temp=d_square(n,m);
        if(m~=n && ~ismember(m,path) && temp < min)
            min=temp;
            min_idx=m;
            cost=P(n,m);
        end
    end
    path=[path, min_idx];
    %pathCost=[pathCost, min];
    pathCost=[pathCost, cost];
end
pathCost=[pathCost, cost];
pathCost_final(1)=pathCost(1);
pathCost_final(N)=pathCost(N);
for n=2:N-1
    pathCost_final(n)=max(pathCost(n),pathCost(n-1));
end
    


%kkk=1;

        