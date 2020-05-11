clc
close all
clear all
l=input('Enter length of rod:');
le=input('Enter Element Lenght:');
n=(l/le)+1;
tend=input('Enter end time:');
ts=input('Enter  time steps:');
tn=(tend/ts)+1;
q1=input('Enter input flux:');
q2=input('Enter Output flux:');
t0=input('Enter initial Temp :');
k=input('Enter conductivity:');
a=k*ts/(le*le)
q1s=q1*le/k;
q2s=q2*le/k;
for s=1:1:n-2
    d1z(s)=2+2*a;
end
for s=1:1:n-3
    d0z(s)=-a;
    d2z(s)=-a;
end
dz=diag(d1z)+diag(d0z,-1)+diag(d2z,1);
for s=1:1:n-2
    d1y(s)=2-2*a;
end
for s=1:1:n-3
    d0y(s)=-a;
    d2y(s)=-a;
end
dy=diag(d1y)+diag(d0y,-1)+diag(d2y,1);
X=inv(dz);
u=zeros(n,tn);
u(2:n-1,1)=t0;

for j=2:1:tn
    u(1,j-1)=u(2,j-1)+q1s;
    u(n,j-1)=u(n-1,j-1)-q2s;
    for i=2:1:n-2
        U(i,1)=u(i,j);
    end
     for i=2:1:n-2
      Ui(i,1)=u(i,j-1);           
     end
    ini=zeros(n-2,1);
    ini(1)=u(1,j-1);
    ini(n-2)=u(n,j-1);
    U=X*dy*Ui+2.*X*a*ini;
    for i=2:1:n-2
        u(i,j)=U(i,1);
    end
    if(j==tn)
        u(1,j)=u(2,j)+q1s;
        u(n,j)=u(n-1,j)-q2s;
    end
    
    N=1:1:n;
    plot(N,u(:,j));
    hold on;
    
end
for i=1:1:n
    fprintf('Step %d',i);
        u(:,j)
    end


        
        
        
    