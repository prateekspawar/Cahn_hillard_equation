clc
close all
clear all

N=100 ; % input('Enter number of nodes:');
L=1 ; 
h=L/N; % meter input('Enter space step:');
t=5 ; % secs   input('Enter final time:');
delta_t=7.08e-5; %input('Enter time step:');

Diff=7.08e-15; % m2/s
R_cons=8.314;  
Temp=300;     % Kelvin

M=1% Diff/(R_cons*Temp);%input('Enter Constant M');
Jp=-0.6168; 

%ci=%input('Enter initial value:');
%j1=%input('Enter initial flux:');
%j2=%input('Enter end flux');

%numloadsteps=t/delta_t+1;

a=-M/(2*h*h);
b= 1/delta_t + M/(h*h);

for i=1:N
    center_diag(i)=b;
end

for i=1:N-1
    upper_diag(i)=a;
end

K=diag(center_diag)+diag(upper_diag,-1)+diag(upper_diag,1);

cn(1:N)=0.0;
ck(1:N)=0.0;

numloadSteps=3000;

for j=1:numloadSteps;

  %  numstep=numstep+j
  
    error=1000;
    tol=0.001;
        
while error>tol
    
cL_n=cn(2);
cL_k=ck(2);

cR_n=cn(N-1)-2*h*Jp/M;
cR_k=ck(N-1)-2*h*Jp/M;

R(1)=(ck(1)-cn(1))/delta_t- ... 
    M/(2*h*h)*(ck(3)-2*ck(2)+ck(1)+cn(3)-2*cn(2)+ck(1));
R(N)=(ck(N)-cn(N))/delta_t- ...
    M/(2*h*h)*(ck(N-1)-2*ck(N)+cR_k+cn(N-1)-2*cn(N)+cR_n);

for i=2:N-1
R(i)=(ck(i)-cn(i))/delta_t- ...
    M/(2*h*h)*(ck(i-1)-2*ck(i)+ck(i+1)+cn(i-1)-2*cn(i)+cn(i+1));
end

dck=-inv(K)*R';
ck=ck+dck';

sumdC=0.0;
sumC=0.0;
for i=1:N;
sumdC=sumdC+dck(i)*dck(i);
sumC=sumC+ck(i)*ck(i);
end
error=sqrt(sumdC/sumC);


end

fprintf('Numloadsteps=%d\n',j)

cn=ck;

end


 x1=0:1.0e-6/(N-1):1e-6;
 t=0;
 q=-1e-4;
 delta_t=0.01;
 time=delta_t:delta_t:numloadSteps;
 for i=numloadSteps;
      t=time(i);
 X=x1./(2*sqrt(Diff*t));
 X1=2*sqrt(t/pi)*exp(-x1.^2/(4*Diff*t));
 X2=(x1./sqrt(Diff)).*erfc(X);
 Y = (-q/sqrt(Diff))*(X1-X2);

 %plot(0.01:0.01:1.0,C(i,1:length(cn)),'-b','linewidth',2); hold on;
 %hold on
 end
 
 Y=sort(Y,'ascend');
 
 plot(x1/1.0e-6,Y/22900,'or'); hold on;
 plot(x1/1.0e-6,cn,'-k');
 
 
