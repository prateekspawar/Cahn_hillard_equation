clc
close all
clear all
m=50 ; % input('Enter number of nodes:');
L=1.0e-6 ; 
dx=L/m; % meter input('Enter space step:');
t=5 ; % secs   input('Enter final time:');
dt=0.01; %input('Enter time step:');

Diff=7.08e-15; % m2/s
R_cons=8.314;  
Temp=300;     % Kelvin

M=Diff/(R_cons*Temp);%input('Enter Constant M');
j2=0.001; 
j1=0.00;
ci=0.00;
%ci=%input('Enter initial value:');
%j1=%input('Enter initial flux:');
%j2=%input('Enter end flux');

n=t/dt+1;

a=-M/(2*dx*dx);
b=1/dt + M/(dx*dx);

for i=1:m
    center_diag(i)=b;
end

for i=1:m-1
    upper_diag(i)=a;
end

K=diag(center_diag)+diag(upper_diag,-1)+diag(upper_diag,1);
ki=inv(K);



c=zeros(m+2,n);
c(2:m+1,1)=ci;

dc=zeros(m,n);
dc(:,:)=0;
R=zeros(m,n);


for j=2:1:n
     c(1,j-1)=c(3,j-1)+2*dx/M*j1;
     c(m+2,j-1)=c(m,j-1)-2*dx*j2/M;
    fprintf('Time step %d\n',j);
    while(1)  
%         fprintf('dc');
%         dc(:,j)
        c(2:m+1,j)=c(2:m+1,j-1)+dc(:,j);     
        c(1,j)=c(3,j)+2*dx/M*j1;
        c(m+2,j)=c(m,j)-2*dx*j2/M;
%         fprintf('c');
%         c(:,j)
        for i=1:1:m
            R(i,j)=r(M,dx,dt,c(i,j),c(i+1,j),c(i+2,j),c(i,j-1),c(i,j-1),c(i,j-1));
        end
%         fprintf('R');
%         R(:,j)
        dc(:,j)=(-1)*K\R(:,j);
        
        rr1=0;rr2=0;    
        for i=1:1:m
            rr1=rr1+dc(i,j)*dc(i,j);
            rr2=rr2+c(i+1,j)*c(i+1,j);
        end
        rr=rr1/rr2;
        rr=sqrt(rr);
        fprintf('rr=%f \n',rr);
        if(rr<0.00001)
            break;
        end
       
     end
    
end
c(:,n)

function re=r(M,dx,dt,cd0,cd1,cd2,c0,c1,c2)
re=(cd1-c1)/dt-M/(2*dx*dx)*(cd2-2*cd1+cd0+c2-2*c1+c0);
    end

