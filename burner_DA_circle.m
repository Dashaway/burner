%181016
%Dashaway
%ȼ�ղ�������
%Բ��װҩ������ȼ��

%190117
%����ע��
%��ΪԲ��ҩ����ҩ����������

clear;
close all;
%��ʼֵ�趨
%����
%��������
long = 100000;        %���鳤��
n = 1:1:long;       %��ͼ������
dt = 5e-5;      %����ֵ
%��֪����
p0 = 1.02e5;    %��ʼѹǿ(Pa)
Dr = 0.132;       %ȼ�����⾶(m)
D = 0.0265;         %ҩ����ʼ�⾶(m)
d = 0.0125;       %ҩ����ʼ�ھ�(m)
ep = (D - d)/2;     %�����(m)
n_s = 17;    %ҩ������
Lp = 0.22;      %װҩ����(m)
gamma = 1.2;        %���ȱ�
rho_p = 1730;       %�ܶ�(kg/m^3)
n_p = 0.302;      %ѹǿָ��
c = 1600;       %�����ٶ�(m/s)
At = 1.884785e-3;     %��ܺ����(m^2)
r0 = 5e-3;      %?��ʼȼ��(m/s)
alpha_r = 1.7e-4;        %?ȼ��ϵ��
phi_alpha = 1;       %?��ʴ����
phi_m = 1;      %?

%����ó���

Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���
%�ܳ�����
s0 = n_s*pi*d;      %ȼ�����ʼ�߳�
    %s =  n_s*pi*(d + e);
%ͨ���������
Ap0 = pi*Dr^2 / 4 - n_s*pi*D^2 / 4 + n_s*pi*d^2 / 4;
    %Ap = pi*Dr^2 / 4 - n_s*pi*D^2 / 4 + n_s*pi*(d + e)^2 / 4;


%����
%�����ñ������Ѹ���ֵ��
p = p0*ones(1,long);        %ʵ��ѹǿ(Pa)
r = r0*ones(1,long);     %ȼ��(m/s)
e = 0*ones(1,long);     %����ȥ���(m)
s = s0*ones(1,long);     %ȼ����ʵ�ʱ߳�(m)
m_b = rho_p*s0*Lp*r0*ones(1,long);     %ȼ��������(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %��������(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %����(?)
Ab = s0*Lp*ones(1,long);     %ȼ�����(m^2)
Ap = Ap0*ones(1,long);     %ͨ�����(m^2)
Vg = Ap*Lp;     %�������(m^3)
%ѭ������
i = 1;


%��ֵ����
%ѭ��ǰ����������
%ѹǿ���в���
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;


%ѭ������
while e <= ep
    %��������𲽼��㣬������µ�ѹǿֵ
    %dp = p_a*(Ab(i) / Vg(i))*p(i)^n_p  - p_b*p(i) / Vg(i);
    k1 = p_a*(Ab(i) / Vg(i))*( p(i)^n_p ) - p_b*p(i) / Vg(i);
    k2 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k1/2)^n_p ) - p_b*(p(i) + dt*k1/2) / Vg(i);
    k3 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k2/2)^n_p ) - p_b*(p(i) + dt*k2/2) / Vg(i);
    k4 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k3)^n_p )  - p_b*(p(i) + dt*k3) / Vg(i);
    dp = (k1 + 2*k2 + 2*k3 + k4)/6;
    i = i + 1;
    p(i) = p(i - 1) + dp*dt;

    %��������
    %�������������������ֵ
    r(i) = alpha_r*p(i)^n_p;
    e(i) = e(i - 1) + r(i)*dt;
    
    s(i) = n_s*pi*(d + e(i));
    Ap(i) = pi*Dr^2 / 4 - n_s*pi*D^2 / 4 + n_s*pi*(d + e(i))^2 / 4;
    Ab(i) = s(i)*Lp;
    Vg(i) = Ap(i)*Lp;
    m_b(i) = rho_p*Ab(i)*r(i);
    m_p(i) = phi_m*p(i)*At / c;
    F(i) = 2000*m_p(i) + 4.4*At*(p(i)-p0);
    if i >= long
        break;
    end
end

%ѭ������
i = i + 1;
n(i:1:long) = [];
p(i:1:long) = [];
r(i:1:long) = [];
e(i:1:long) = [];
s(i:1:long) = [];
m_b(i:1:long) = [];
m_p(i:1:long) = [];
F(i:1:long) = [];
Ab(i:1:long) = [];
Ap(i:1:long) = [];
Vg(i:1:long) = [];
t = dt*n;

%�������

figure;
hold on;
subplot(2,1,1);
plot(t,p);
grid on,axis tight ;
title('ȼ����ѹ��')

subplot(2,1,2);
plot(t,e);
grid on,axis tight ;
title('��ȥ���')

figure;
hold on;
subplot(2,1,1);
plot(t,Ab);
grid on,axis tight ;
title('ȼ�����')

subplot(2,1,2);
plot(t,Ap);
grid on,axis tight ;
title('ͨ�����')

figure;
hold on;
subplot(2,1,1);
plot(t,m_b);
grid on,axis tight ;
title('ȼ��������')

subplot(2,1,2);
plot(t,m_p);
grid on,axis tight ;
title('��������')

figure;
hold on;
plot(t,F);
grid on,axis tight ;
title('����')

figure;
hold on;
plot(t,p);
grid on,axis tight ;
title('ȼ����ѹ��')

