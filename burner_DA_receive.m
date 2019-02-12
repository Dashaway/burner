%180916
%Dashaway
%ȼ�ղ�������

%dp���㹫ʽ����
%���ֲ����޺�����Դ

clear;
close all;
%��ʼֵ�趨
%����
%��������
long = 5000;        %���鳤��
n = 1:1:long;       %��ͼ������
t = 0;      %ʱ��
dt = 0.01;      %����ֵ
%��֪����
p0 = 1.02e5;    %��ʼѹǿ
D = 0.132;       %�⾶
d = 0.06;       %��ʼ�ھ�
ep = (D - d)/2;     %�����
Lp = 0.22;      %װҩ����
gamma = 1.2;        %���ȱ�
rho_p = 1730;       %�ܶ�
np = 0.302;      %ѹǿָ��
c = 1600;       %�����ٶ�
alpha = 1.41e-3;
phi = 1; 
%����ó���
Gamma = ( (2/(gamma + 1) )^( (gamma + 1)/(2*(gamma - 1) ) ) )*sqrt(gamma);      %���ȱȺ���
%δȷ������
phi_m = 1;
At = 0.002;

%����
%�����ñ������Ѹ���ֵ��
e = 0*ones(1,long);     %����ȥ���
p = p0*ones(1,long);        %ʵ��ѹǿ
Ab = pi*d*ones(1,long);     %ȼ�����
Ap = (pi/4)*(d^2)*ones(1,long);     %ͨ�����
Vg = Ap*Lp;     %�������
%ѭ������
i = 1;


%��ֵ����
%ѭ��ǰ����
a = rho_p*1.41e-3*Gamma^2*c^2*phi;
b = 1*Gamma^2*c*At;

%ѭ������
while e <= ep
    dp = 5e4;
    %dp = -1e4;
    %dp = a*Ab(i)*p(i)^np/Vg(i)-b*p(i)/Vg(i);
    i = i + 1;
    p(i) = p(i - 1) + Runge_Kutta_4th(dp,dt);
    e(i) = e(i - 1) + alpha*p(i)^np*dt;
    %p(i) = p(i - 1) + dp;
    %e(i) = e(i - 1) + 1e-3;
    %p(i) = p(i - 1) + Runge_Kutta_4th(dp,dt);
    %e(i) = e(i - 1) + alpha*p(i)^np*dt;
    
    %��������
    Ab(i) = pi*(d+2*e(i) )*Lp;
    Ap(i) = (pi/4)*(d+2*e(i) )^2;
    Vg(i) = Ap(i)*Lp;
    if i >= long
        break;
    end
end

%ѭ������
i = i + 1;
n(i:1:long) = [];
e(i:1:long) = [];
p(i:1:long) = [];
Ab(i:1:long) = [];
Ap(i:1:long) = [];
Vg(i:1:long) = [];


%�������

figure;
hold on;
subplot(2,1,1);
plot(n,p);
grid on,axis tight ;
title('ȼ����ѹ��')

subplot(2,1,2);
plot(n,e);
grid on,axis tight ;
title('��ȥ���')

figure;
hold on;
subplot(2,1,1);
plot(n,Ab);
grid on,axis tight ;
title('ȼ�����')

subplot(2,1,2);
plot(n,Ap);
grid on,axis tight ;
title('ͨ�����')



