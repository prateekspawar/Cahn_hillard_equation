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
Jp=-0.6168; 

cn(1:N)=0.1;
ck(1:N)=0.1;

numloadSteps=1000;
r=h:h:1;
for j=1:numloadSteps

  %  numstep=numstep+j
  
    error=1000;
    tol=0.001;
        
while error>tol
    
 
cL_n=cn(2);
cL_k=ck(2);

cR_n=cn(N-1)-2*h*Jp/(1-2*chi*cn(N)*(1-cn(N)));
cR_k=ck(N-1)-2*h*Jp/(1-2*chi*ck(N)*(1-ck(N)));

ul_n=mu(chi,cL_n);
ul_k=mu(chi,cL_k);

ur_n=mu(chi,cR_n);
ur_k=mu(chi,cR_k);

Rc(1)=(ck(1)-cn(1))/delta_t- ... 
    1/2*((mu(chi,ck(2))-ul_k)/(2*h)*(1+A*ck(1)*...
    (1-ck(1))/r(1))-2*ck(1)*(ck(2)-cL_k)*(mu(chi,ck(2))-ul_k)/(4*h*h)...
    +(ck(1)-ck(1)*ck(1))/(h*h)*(mu(chi,ck(2))-2*mu(chi,ck(1))+ul_k)...
    +(mu(chi,cn(2))-ul_n)/(2*h)*(1+A*cn(1)*...
    (1-cn(1))/r(1))-2*cn(1)*(cn(2)-cL_n)*(mu(chi,cn(2))-ul_n)/(4*h*h)...
    +(cn(1)-cn(1)*cn(1))/(h*h)*(mu(chi,cn(2))-2*mu(chi,cn(1))+ul_n));
Rc(N)=(ck(N)-cn(N))/delta_t- ... 
    1/2*((ur_k-mu(chi,ck(N-1)))/(2*h)*(1+A*ck(N)*...
    (1-ck(N))/r(N))-2*ck(N)*(cR_k-ck(N-1))*(ur_k-mu(chi,ck(N-1)))/(4*h*h)...
    +(ck(N)-ck(N)*ck(N))/(h*h)*(mu(chi,ck(N-1))-2*mu(chi,ck(N))+ur_k)...
    +(ur_n-mu(chi,cn(N-1)))/(2*h)*(1+A*cn(N)*...
    (1-cn(N))/r(N))-2*cn(N)*(cR_n-cn(N-1))*(ur_n-mu(chi,cn(N-1)))/(4*h*h)...
    +(cn(N)-cn(N)*cn(N))/(h*h)*(mu(chi,cn(N-1))-2*mu(chi,cn(N))+ur_n));
for i=2:N-1
Rc(i)=(ck(i)-cn(i))/delta_t- ... 
    1/2*((mu(chi,ck(i+1))-mu(chi,ck(i-1)))/(2*h)*(1+A*ck(i)*...
    (1-ck(i))/r(i))-2*ck(i)*(ck(i+1)-ck(i-1))*(mu(chi,ck(i+1))-mu(chi,ck(i-1)))/(4*h*h)...
    +(ck(i)-ck(i)*ck(i))/(h*h)*(mu(chi,ck(i+1))-2*mu(chi,ck(i))+mu(chi,ck(i-1)))...
    +(mu(chi,cn(i+1))-mu(chi,cn(i-1)))/(2*h)*(1+A*cn(i)*...
    (1-cn(i))/r(i))-2*cn(i)*(cn(i+1)-cn(i-1))*(mu(chi,cn(i+1))-mu(chi,cn(i-1)))/(4*h*h)...
    +(cn(i)-cn(i)*cn(i))/(h*h)*(mu(chi,cn(i+1))-2*mu(chi,cn(i))+mu(chi,cn(i-1))));
end


for i=1:N
    if i==1
    center_diag(i)=1/delta_t-1/2*((mu(chi,ck(2))-ul_k)/(2*h*r(1))*A*(1-2*ck(1))...
            -2*(ck(2)-cL_k)*(mu(chi,ck(2))-ul_k)/(4*h*h)+(1-2*ck(1))/(h*h)*...
            (mu(chi,ck(2))-2*mu(chi,ck(1))+ul_k));
    
    else if i==N
        center_diag(i)=1/delta_t-1/2*((ur_k-mu(chi,ck(N-1)))/(2*h*r(N))*A*(1-2*ck(N))...
            -2*(cR_k-ck(N-1))*(ur_k-mu(chi,ck(N-1)))/(4*h*h)+(1-2*ck(N))/(h*h)*...
            (mu(chi,ck(N-1))-2*mu(chi,ck(N))+ur_k));
        else
            center_diag(i)=1/delta_t-1/2*((mu(chi,ck(i+1))-mu(chi,ck(i-1)))/(2*h*r(i))*A*(1-2*ck(i))...
            -2*(ck(i+1)-ck(i-1))*(mu(chi,ck(i+1))-mu(chi,ck(i-1)))/(4*h*h)+(1-2*ck(i))/(h*h)*...
            (mu(chi,ck(i+1))-2*mu(chi,ck(i))+mu(chi,ck(i-1))));
        end
    end
                
end

for i=1:N-1
    if i==1
    upper_diag(i)=1/(4*h*h)*(ck(1)*(mu(chi,ck(2))-ul_k));
    lower_diag(i)=(-1)/(4*h*h)*ck(i+1)*(mu(chi,ck(i+2))-mu(chi,ck(i)));
    else if i==(N-1)
        lower_diag(i)=(-1)/(4*h*h)*ck(N)*(ul_k-mu(chi,ck(N-1)));
        upper_diag(i)=(1)/(4*h*h)*ck(i)*(mu(chi,ck(i+1))-mu(chi,ck(i-1)));
        else
            lower_diag(i)=(-1)/(4*h*h)*ck(i+1)*(mu(chi,ck(i+2))-mu(chi,ck(i)));
             upper_diag(i)=(1)/(4*h*h)*ck(i)*(mu(chi,ck(i+1))-mu(chi,ck(i-1)));
        end
    end
end

K=diag(center_diag)+diag(lower_diag,-1)+diag(upper_diag,1);
dck=-inv(K)*Rc';
ck=ck+dck';

sumdC=0.0;
sumC=0.0;

for i=1:N
sumdC=sumdC+dck(i)*dck(i);
sumC=sumC+ck(i)*ck(i);
end
error=sqrt(sumdC/sumC);


end

fprintf('Numloadsteps=%d\n',j)

cn=ck;

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
 plot(cn,'-k');
 
 function u=mu(chi,c)
 u=chi*(1-2*c)+log(c/(1-c));
 end
 
