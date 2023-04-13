clc
clear
close all

% interpolation for f in x
% you can change f , x and custmize code

x = -2*pi:0.4:2*pi;
f = sin(x) + exp(.2*x) .* cos(2*x);

min_x=min(x);
max_x=max(x);
n=length(x);

% Implementation Natural Cubic Spline Method

M(1)=0; M(n)=0;
h=x(2:end)-x(1:end-1);
ff=(f(2:end)-f(1:end-1))./h;
F=6*(ff(2:end)-ff(1:end-1));
for i=1:n-2
   A(i,i)=2*(h(i)+h(i+1));
   if i<n-2
        A(i,i+1)=h(i+1);
        A(i+1,i)=h(i+1);
   end
end
M(2:n-1)=A/F;
for j=1:n-1
   t1=[1, -x(j)];
   t2=[-1, x(j+1)];
   T1=conv(conv(t1,t1),t1);
   T2=conv(conv(t2,t2),t2);
   s(j,:)=1/h(j)*(T2*M(j)/6+T1*M(j+1)/6+...
       (f(j)-h(j)^2*M(j)/6)*[zeros(1,length(T2)-length(t2)),t2]+ ...
       (f(j+1)-h(j)^2*M(j+1)/6)*[zeros(1,length(T1)-length(t1)),t1]);
end
xxplot=[]; Splot=[];
for j=1:n-1
    xplot(j,:)=linspace(x(j),x(j+1),n);
    splot(j,:)=polyval(s(j,:),xplot(j,:));
    xxplot=[xxplot xplot(j,:)];
    Splot=[Splot splot(j,:)];
end

% plot results
t=linspace(min_x,max_x,n-1);
for i=1:n-1
   fx(i)=f(i)+(f(i+1)-f(i))/(x(i+1)-x(i))*((x(i+1)+x(i))/2-x(i));
end
d=Splot(fix(n/2)+1:n:end);
e=fx-d; % find error
m=min_x+.5:max_x-.5;
figure; % plot f(x) and S(x)
plot(x,f,'-ro',xxplot,Splot,'Linewidth',1.25)
grid on
title('f(x) & S(x)');
legend('support points','Natural Cubic Spline interpolant S(x)',2)
% figure; % plot error
% plot(t,e,'Linewidth',1.25) 
% title('Error');
% grid on
% also you can type Splot in command window and see S(x) values

