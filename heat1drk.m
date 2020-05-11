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
a=k*ts/(le*le)/2
q1s=q1*le/k
q2s=q2*le/k
u=zeros(n,tn);
ud=zeros(n,tn-1);
for i=2:1:n-1
    u(i,1)=t0;
end
for j=2:1:tn
    
        u(1,j-1)=u(2,j-1)+q1s;
        u(n,j-1)=u(n-1,j-1)-q2s;
        for i=2:1:n-1
        ud(i,j-1)=a*u(i+1,j-1)+(1-2*a)*u(i,j-1)+a*u(i-1,j-1);
        end
        ud(1,j-1)=ud(2,j-1)+q1s;
        ud(n,j-1)=ud(n-1,j-1)-q2s;
        for i=2:1:n-1
       u(i,j)=a*ud(i+1,j-1)+(1-2*a)*ud(i,j-1)+ a*ud(i-1,j-1);
        end
          
end
u(1,tn)=u(2,tn)+q1s;
u(n,tn)=u(n-1,tn)-q2s;


for j=1:1:tn
    fprintf('Step %d',j);
        u(:,j)
end
for j=1:1:tn-1
    fprintf('Intertep %d',j);
        ud(:,j)
end
for j=1:1:tn
    N=1:1:n;
    plot(N,u(:,j));
    hold on;
end


        
        
        
    