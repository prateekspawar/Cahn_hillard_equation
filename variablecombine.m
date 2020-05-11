clc
close all
clear all

N=50 ;
L=1 ; 
h=L/N; % meter 
chi=0.0;   
delta_t=7.08e-5; 
A=1; % 1 for cylinder,2 for sphere
Diff=7.08e-15; % m2/s
R_cons=8.314;  
Temp=300;     % Kelvin

M=1;% Diff/(R_cons*Temp);%input('Enter Constant M');
ka=0.01; %kappa
j1_1=0.0;
j3_1=0.0;
j1_n=-0.6168;
j3_n=0.0;

cn(3:N+2)=0.1;
ck(3:N+2)=0.1;

numloadSteps=100;
r1=h:h:1;
r(3:N+2)=r1(1:N);
R=zeros(N+2);
for j=1:numloadSteps

  %  numstep=numstep+j
  
    error=1000;
    tol=0.1;
        
while error>tol


J1_n=0;
Jn_n=-0.6168;    
J1_k=0;
Jn_k=-0.6168;    

cn(1)=(2*h*h*h)/(ka*cn(3)*(1-cn(3)))*(j1_1*2*chi*cn(3)*(1-cn(3))-j1_1-J1_n)...
    +cn(5)-4*(j1_1)*h ;
cn(2)=cn(4)+2*h/(2*chi*cn(3)*(1-cn(3))-1)*(j3_1*ka*cn(3)*(1-cn(3))-J1_n);
cn(N+4)=(-2*h*h*h)/(ka*cn(N+2)*(1-cn(N+2)))*(2*j1_n*chi*cn(N+2)*(1-cn(N+2))-Jn_n)...
    +cn(N)+4*h*j1_n;
cn(N+3)=cn(N+1)-(2*h)/(2*chi*cn(N+2)*(1-cn(N+2))-1)*(ka*j3_n*cn(N+2)*(1-cn(N+2)))-Jn_n;

ck(1)=(2*h*h*h)/(ka*ck(3)*(1-ck(3)))*(j1_1*2*chi*ck(3)*(1-ck(3))-j1_1-J1_k)...
    +ck(5)-4*(j1_1)*h ;
ck(2)=ck(4)+2*h/(2*chi*ck(3)*(1-ck(3))-1)*(j3_1*ka*ck(3)*(1-ck(3))-J1_k);
ck(N+4)=(-2*h*h*h)/(ka*ck(N+2)*(1-ck(N+2)))*(2*j1_n*chi*ck(N+2)*(1-ck(N+2))-Jn_k)...
    +ck(N)+4*h*j1_n;
ck(N+3)=ck(N+1)-(2*h)/(2*chi*ck(N+2)*(1-ck(N+2))-1)*(ka*j3_n*ck(N+2)*(1-ck(N+2)))-Jn_k;

for i=3:N+2
R(i)=(ck(i)-cn(i))/delta_t+ ... 
    1/2*(-A/r(i)*(1-2*chi*ck(i)*(1-ck(i)))*(ck(i+1)-ck(i-1))/2/h...
    +2*chi*(1-2*ck(i))*((ck(i+1)-ck(i-1))/(2*h))^2 ...
    -(1-2*chi*ck(i)*(1-ck(i)))*((ck(i+1)-2*ck(i)+ck(i-1))/(h*h)) ...
    +A/r(i)*ka*ck(i)*(1-ck(i))*((ck(i+2)-2*ck(i+1)+2*ck(i-1)-ck(i-2))/(2*h*h*h)) ...
    +ka*(1-2*ck(i))*(ck(i+1)-ck(i-1))/2/h*((ck(i+2)-2*ck(i+1)+2*ck(i-1)-ck(i-2))/(2*h*h*h)) ...
    +ka*ck(i)*(1-ck(i))*(ck(i+2)-4*ck(i+1)+6*ck(i)-4*ck(i-1)+ck(i-2))/(h^4) ...    
    -A/r(i)*(1-2*chi*cn(i)*(1-cn(i)))*(cn(i+1)-cn(i-1))/2/h...
    +2*chi*(1-2*cn(i))*((cn(i+1)-cn(i-1))/(2*h))^2 ...
    -(1-2*chi*cn(i)*(1-cn(i)))*((cn(i+1)-2*cn(i)+cn(i-1))/(h*h)) ...
    +A/r(i)*ka*cn(i)*(1-cn(i))*((cn(i+2)-2*cn(i+1)+2*cn(i-1)-cn(i-2))/(2*h*h*h)) ...
    +ka*(1-2*cn(i))*(cn(i+1)-cn(i-1))/2/h*((cn(i+2)-2*cn(i+1)+2*cn(i-1)-cn(i-2))/(2*h*h*h)) ...
    +ka*cn(i)*(1-cn(i))*(cn(i+2)-4*cn(i+1)+6*cn(i)-4*cn(i-1)+cn(i-2))/(h^4));
end
Rc(1:N)=R(3:N+2);

for i=1:N
    center(i)=1/delta_t+1/2*(-A/r(i+2)*(ck(i+1+2)-ck(i-1+2))/2/h*(2*chi*(1-ck(i+2))) ...
        +2*chi*(-2)*(ck(i+1+2)-ck(i-1+2))/2/h+(2*chi*(1-2*ck(i+2)))/h/h*(ck(i+1+2)-2*ck(i+2)+ck(i-1+2)) ...
        +(2/h/h)*(1-2*chi*ck(i+2)*(1-ck(i+2)))+A/r(i+2)*ka*(1-2*ck(i+2))/(2*h^3)*(ck(i+2+2)-2*ck(i+1+2)+2*ck(i-1+2)-ck(i-2+2)) ...
        -2*ka*(ck(i+1+2)-ck(i-1+2))/(2*h)*(ck(i+2+2)-2*ck(i+1+2)+2*ck(i-1+2)-ck(i-2+2))/(2*h^3) ...
        +ka*(1-2*ck(i+2))/(h^4)*(ck(i+2+2)-4*ck(i+1+2)+6*ck(i+2)-4*ck(i-1+2)+ck(i-2+2)) ...
        +ka*ck(i+2)*(1-ck(i+2))*6/(h^4));
end
for i=1:N-1
    up_1(i)=1/2*(-A/r(i+2)/(2*h)*(1-2*chi*ck(i+2)*(1-ck(i+2)))+2*chi*2*(1-2*ck(i+2))/(2*h)*(ck(i+1+2)-ck(i-1+2))/(2*h) ...
       -(1-2*chi*ck(i+2)*(1-ck(i+2)))/(h*h)+A/r(i+2)*ka*ck(i+2)*(1-ck(i+2))*(-2)/(2*h^3) ...
       +ka*(1-2*ck(i+2))*(ck(i+1+2)-ck(i-1+2))*(-2)/(4*h^4)+ka*(1-2*ck(i+2))/(2*h^3)/(2*h) ...
       *(ck(i+2+2)-2*ck(i+1+2)+2*ck(i-1+2)-ck(i-2+2))+ka*ck(i+2)*(1-ck(i+2))*(-4)/(h^4));
end
for i=2:N
    down_11(i)=1/2*(-A/r(i+2)*(1-2*chi*ck(i+2)*(1-ck(i+2)))*(-1/2/h)+2*chi*(1-2*ck(i+2))/(h)*(-1/2/h) ...
        *(ck(i+1+2)-ck(i-1+2))-(1-2*chi*ck(i+2)*(1-ck(i+2)))/(h*h)+A/r(i+2)*ka*ck(i+2)*(1-ck(i+2))/(h^3) ...
        +ka*(1-2*ck(i+2))/(2*h)/(h^3)*(ck(i+1+2)-ck(i-1+2))+ka*(1-2*ck(i+2))/(2*h^3)*(-1/2/h) ...
        *(ck(i+2+2)-2*ck(i+1+2)+2*ck(i-1+2)-ck(i-2+2))+ka*ck(i+2)*(1-ck(i+2))*(-4/(h^4)));
end
for i=1:N-2
    up_2(i)=1/2*(A/r(i+2)*ka*ck(i+2)*(1-ck(i+2))/(2*h^3)+ka*(1-2*ck(i+2))*(ck(i+1+2)-ck(i-1+2))/(4*h^4) ...
        +ka*ck(i+2)*(1-ck(i+2))/(h^4));
end
for i=3:N
    down_22(i)=1/2*(-A/r(i+2)*ka*ck(i+2)*(1-ck(i+2))/(2*h^3)-ka*(1-2*ck(i+2))*(ck(i+1+2)-ck(i-1+2))/(4*h^4) ...
        +ka*ck(i+2)*(1-ck(i+2))/(h^4));
end
down_1(1:N-1)=down_11(2:N);
down_2(1:N-2)=down_22(3:N);

K=diag(center)+diag(down_1,-1)+diag(up_1,1)+diag(up_2,2)+diag(down_2,-2);
dck=-inv(K)*Rc';
CK(1:N)=ck(3:N+2);
CK=CK+dck';
ck(3:N+2)=CK(1:N);

sumdC=0.0;
sumC=0.0;

for i=1:N
sumdC=sumdC+dck(i)*dck(i);
sumC=sumC+CK(i)*CK(i);
end
error=sqrt(sumdC/sumC);
fprintf('Error=%f \n',error);

end

fprintf('Numloadsteps=%d\n',j)

cn=ck;
CN=CK;
end


%  x1=0:1.0e-6/(N-1):1e-6;
%  t=0;
%  q=-1e-4;
%  delta_t=0.01;
%  time=delta_t:delta_t:numloadSteps;
%  for i=numloadSteps;
%       t=time(i);
%  X=x1./(2*sqrt(Diff*t));
%  X1=2*sqrt(t/pi)*exp(-x1.^2/(4*Diff*t));
%  X2=(x1./sqrt(Diff)).*erfc(X);
%  Y = (-q/sqrt(Diff))*(X1-X2);
% 
%  %plot(0.01:0.01:1.0,C(i,1:length(cn)),'-b','linewidth',2); hold on;
%  %hold on
%  end
%  
%  Y=sort(Y,'ascend');
%  
%  plot(x1/1.0e-6,Y/22900,'or'); hold on;
 plot(CN,'-k');
 
 function u=mu(chi,c)
 u=chi*(1-2*c)+log(c/(1-c));
 end
 
