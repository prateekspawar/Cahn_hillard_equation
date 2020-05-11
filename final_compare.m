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
j1_1=0.0;
j1_n=-0.6168;
cn=zeros(N,1);
ck=zeros(N,1);
cn(1:N)=0.1;
ck(1:N)=0.1;
numloadSteps=500;
r=h:h:1;
R=zeros(N,1);
for j=1:numloadSteps

  %  numstep=numstep+j
  
    error=1000;
    tol=0.001;
        
while abs(error)>tol
cn_l=cn(2);
cn_r=cn(N-1)+(2*h)*j1_n/(2*chi*cn(N)*(1-cn(N))-1);
ck_l=ck(2);
ck_r=ck(N-1)+(2*h)*j1_n/(2*chi*ck(N)*(1-ck(N))-1);

R(1)=(ck(1)-cn(1))/delta_t+ ... 
    1/2*(-A/r(1)*(1-2*chi*ck(1)*(1-ck(1)))*(ck(2)-ck_l)/2/h...
    +2*chi*(1-2*ck(1))*((ck(2)-ck_l)/(2*h))^2 ...
    -(1-2*chi*ck(1)*(1-ck(1)))*((ck(2)-2*ck(1)+ck_l)/(h*h)) ...
    -A/r(1)*(1-2*chi*cn(1)*(1-cn(1)))*(cn(2)-cn_l)/2/h...
    +2*chi*(1-2*cn(1))*((cn(2)-cn_l)/(2*h))^2 ...
    -(1-2*chi*cn(1)*(1-cn(1)))*((cn(2)-2*cn(1)+cn_l)/(h*h)));

R(N)=(ck(N)-cn(N))/delta_t+ ... 
    1/2*(-A/r(N)*(1-2*chi*ck(N)*(1-ck(N)))*(ck_r-ck(N-1))/2/h...
    +2*chi*(1-2*ck(N))*((ck_r-ck(N-1))/(2*h))^2 ...
    -(1-2*chi*ck(N)*(1-ck(N)))*((ck_r-2*ck(N)+ck(N-1))/(h*h)) ...
    -A/r(N)*(1-2*chi*cn(N)*(1-cn(N)))*(cn_r-cn(N-1))/2/h...
    +2*chi*(1-2*cn(N))*((cn_r-cn(N-1))/(2*h))^2 ...
    -(1-2*chi*cn(N)*(1-cn(N)))*((cn_r-2*cn(N)+cn(N-1))/(h*h)));

for i=2:N-1
R(i)=(ck(i)-cn(i))/delta_t+ ... 
    1/2*(-A/r(i)*(1-2*chi*ck(i)*(1-ck(i)))*(ck(i+1)-ck(i-1))/2/h...
    +2*chi*(1-2*ck(i))*((ck(i+1)-ck(i-1))/(2*h))^2 ...
    -(1-2*chi*ck(i)*(1-ck(i)))*((ck(i+1)-2*ck(i)+ck(i-1))/(h*h)) ...
    -A/r(i)*(1-2*chi*cn(i)*(1-cn(i)))*(cn(i+1)-cn(i-1))/2/h...
    +2*chi*(1-2*cn(i))*((cn(i+1)-cn(i-1))/(2*h))^2 ...
    -(1-2*chi*cn(i)*(1-cn(i)))*((cn(i+1)-2*cn(i)+cn(i-1))/(h*h)));
end


for i=1:N
    if(i==1)
        center(i)=1/delta_t+1/2*(-A/r(1)*(ck(2)-ck_l)/2/h*(2*chi*(1-ck(1))) ...
        +2*chi*(-2)*(ck(2)-ck_l)/2/h+(2*chi*(1-2*ck(1)))/h/h*(ck(2)-2*ck(1)+ck_l) ...
        +(2/h/h)*(1-2*chi*ck(1)*(1-ck(1))));
    else
        if(i==N)
            center(i)=1/delta_t+1/2*(-A/r(N)*(ck_r-ck(N-1))/2/h*(2*chi*(1-ck(N))) ...
             +2*chi*(-2)*(ck_r-ck(N-1))/2/h+(2*chi*(1-2*ck(N)))/h/h*(ck_r-2*ck(N)+ck(N-1)) ...
            +(2/h/h)*(1-2*chi*ck(N)*(1-ck(N))));
        else
            center(i)=1/delta_t+1/2*(-A/r(i)*(ck(i+1)-ck(i-1))/2/h*(2*chi*(1-ck(i))) ...
             +2*chi*(-2)*(ck(i+1)-ck(i-1))/2/h+(2*chi*(1-2*ck(i)))/h/h*(ck(i+1)-2*ck(i)+ck(i-1)) ...
             +(2/h/h)*(1-2*chi*ck(i)*(1-ck(i))));
        end
    end
    
end
for i=1:N-1
    if(i==1)
        up(i)=1/2*(-A/r(1)/(2*h)*(1-2*chi*ck(1)*(1-ck(1)))+2*chi*2*(1-2*ck(1))/(2*h)*(ck(2)-ck_l)/(2*h) ...
       -(1-2*chi*ck(1)*(1-ck(1)))/(h*h));
    else
        up(i)=1/2*(-A/r(i)/(2*h)*(1-2*chi*ck(i)*(1-ck(i)))+2*chi*2*(1-2*ck(i))/(2*h)*(ck(i+1)-ck(i-1))/(2*h) ...
       -(1-2*chi*ck(i)*(1-ck(i)))/(h*h));
    end
end
for i=2:N
    if(i==N)
       down(i)=1/2*(-A/r(i)*(1-2*chi*ck(i)*(1-ck(i)))*(-1/2/h)+2*chi*(1-2*ck(i))/(h)*(-1/2/h) ...
        *(ck_r-ck(i-1))-(1-2*chi*ck(i)*(1-ck(i)))/(h*h));
    else
        down(i)=1/2*(-A/r(i)*(1-2*chi*ck(i)*(1-ck(i)))*(-1/2/h)+2*chi*(1-2*ck(i))/(h)*(-1/2/h) ...
        *(ck(i+1)-ck(i-1))-(1-2*chi*ck(i)*(1-ck(i)))/(h*h));
    end
end
down1(1:N-1)=down(2:N);
K=diag(center)+diag(down1,-1)+diag(up,1);
dck=-inv(K)*R;
ck=ck+dck;
sumdC=0.0;
sumC=0.0;

for i=1:N
sumdC=sumdC+dck(i)*dck(i);
sumC=sumC+ck(i)*ck(i);
end
error=sqrt(sumdC/sumC);
fprintf('error=%f\n',error);

end

fprintf('Numloadsteps=%d\n',j)

cn=ck;

end

 x1=0:1.0e-6/(N-1):1e-6;
 t=0;
 q=-1e-4;
 delta_t=0.01;
 time=delta_t:delta_t:numloadSteps;
 for i=numloadSteps
      t=time(i);
 X=x1./(2*sqrt(Diff*t));
 X1=2*sqrt(t/pi)*exp(-x1.^2/(4*Diff*t));
 X2=(x1./sqrt(Diff)).*erfc(X);
 Y = (-q/sqrt(Diff))*(X1-X2);

 %plot(0.01:0.01:1.0,C(i,1:length(cn)),'-b','linewidth',2); hold on;
 %hold on
 end
 
 Y=sort(Y,'ascend');
 figure('Name','compare')
% plot(x1/1.0e-6,Y/22900,'or','DisplayName','exact'); hold on;
 plot(x1/1.0e-6,cn,'-k','DisplayName','combined 2nd order');hold on;
 s=sep();

 function s=sep()
 
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
     Jdd_1=0.0;
     Jdd_n=0.0;
     cn=zeros(N+4,1);
     ck=zeros(N+4,1);
     CN=zeros(N,1);
     CK=zeros(N,1);
     cn(3:N+2)=0.1;
     ck(3:N+2)=0.1;
     CK=ck(3:N+2);
     CN=cn(3:N+2);
     ka=0.00;
     numloadSteps=500;
     r(3:N+2)=h:h:1;
     R=zeros(N,1);
     for j=1:numloadSteps

         %  numstep=numstep+j

         error=1000;
         tol=0.001;

         while error>tol


             cn(2)=cn(4);
             ck(2)=ck(4);
             cn(1)=cn(5)-2*cn(4)+2*cn(2)-Jdd_1*2*h^3;
             ck(1)=ck(5)-2*ck(4)+2*ck(2)-Jdd_1*2*h^3;

             cn(N+3)=cn(N+1)-2*h*Jp/(1-2*chi*cn(N+2)*(1-cn(N+2)));
             cn(N+4)=cn(N)-2*cn(N+1)+2*cn(N+3)+Jdd_n*2*h^3;

             ck(N+3)=ck(N+1)-2*h*Jp/(1-2*chi*ck(N+2)*(1-ck(N+2)));
             ck(N+4)=ck(N)-2*ck(N+1)+2*ck(N+3)+Jdd_n*2*h^3;


             for i=3:N+2
                 Rc(i)=(ck(i)-cn(i))/delta_t- ...
                     1/2*((mu(chi,ck(i+1),ck(i+2),ck(i),ka,h)-mu(chi,ck(i-1),ck(i),ck(i-2),ka,h))/(2*h)*(1+A*ck(i)*...
                     (1-ck(i))/r(i))-2*ck(i)*(ck(i+1)-ck(i-1))*(mu(chi,ck(i+1),ck(i+2),ck(i),ka,h)-mu(chi,ck(i-1),ck(i),ck(i-2),ka,h))/(4*h*h)...
                     +(ck(i)-ck(i)*ck(i))/(h*h)*(mu(chi,ck(i+1),ck(i+2),ck(i),ka,h)-2*mu(chi,ck(i),ck(i+1),ck(i-1),ka,h)+mu(chi,ck(i-1),ck(i),ck(i-2),ka,h))...
                     +(mu(chi,cn(i+1),cn(i+2),cn(i),ka,h)-mu(chi,cn(i-1),cn(i),cn(i-2),ka,h))/(2*h)*(1+A*cn(i)*...
                     (1-cn(i))/r(i))-2*cn(i)*(cn(i+1)-cn(i-1))*(mu(chi,cn(i+1),cn(i+2),cn(i),ka,h)-mu(chi,cn(i-1),cn(i),cn(i-2),ka,h))/(4*h*h)...
                     +(cn(i)-cn(i)*cn(i))/(h*h)*(mu(chi,cn(i+1),cn(i+2),cn(i),ka,h)-2*mu(chi,cn(i),cn(i+1),cn(i-1),ka,h)+mu(chi,cn(i-1),cn(i),cn(i-2),ka,h)));
             end
             R(1:N)=Rc(3:N+2);
             for i=3:N+2
                 center_diag(i)=1/delta_t-1/2*((mu(chi,ck(i+1),ck(i+2),ck(i),ka,h)-mu(chi,ck(i-1),ck(i),ck(i-2),ka,h))/(2*h*r(i))*A*(1-2*ck(i))...
                     -2*(ck(i+1)-ck(i-1))*(mu(chi,ck(i+1),ck(i+2),ck(i),ka,h)-mu(chi,ck(i-1),ck(i),ck(i-2),ka,h))/(4*h*h)+(1-2*ck(i))/(h*h)*...
                     (mu(chi,ck(i+1),ck(i+2),ck(i),ka,h)-2*mu(chi,ck(i),ck(i+1),ck(i-1),ka,h)+mu(chi,ck(i-1),ck(i),ck(i-2),ka,h)));

             end
             center(1:N)=center_diag(3:N+2);
             for i=3:N-1+2
                 upper_diag(i)=1/(4*h*h)*(ck(i)*(mu(chi,ck(i+1),ck(i+2),ck(i),ka,h)-mu(chi,ck(i-1),ck(i),ck(i-2),ka,h)));
             end
             upper(1:N-1)=upper_diag(3:N+1);
             for i=4:N+2
                 lower_diag(i)=(-1)/(4*h*h)*(ck(i)*(mu(chi,ck(i+1),ck(i+2),ck(i),ka,h)-mu(chi,ck(i-1),ck(i),ck(i-2),ka,h)));
             end
             lower(1:N-1)=lower_diag(4:N+2);
             K=diag(center)+diag(lower,-1)+diag(upper,1);
             dck=-inv(K)*R;
             CK=CK+dck;
             ck(3:N+2)=CK;

             sumdC=0.0;
             sumC=0.0;

             for i=1:N
                 sumdC=sumdC+dck(i)*dck(i);
                 sumC=sumC+CK(i)*CK(i);
             end
             error=sqrt(sumdC/sumC);


         end

         fprintf('Numloadsteps=%d\n',j)

         CN=CK;
         cn(3:N+2)=CN;
     end

    x1=0:1.0e-6/(N-1):1e-6;
     plot(x1/1.0e-6,CN,'-r');
     legend('show','DisplayName','separate variables')
     s=1;
 end
 function u=mu(chi,c,a,b,ka,h)
 u=chi*(1-2*c)+log(c/(1-c))-ka*((a+b-2*c)/(2*h));
 end
 