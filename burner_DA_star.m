%181016
%Dashaway
%ȼ�ղ�������
%����װҩ������ȼ��

%181114
%����ȼ�������ʺ�������������
%190117
%����ע��

clear;
close all;
%��ʼֵ�趨
%����
%��������
long = 50000;        %���鳤��
n = 1:1:long;       %��ͼ������
dt = 5e-5;      %����ֵ
%��֪����
p0 = 1.02e5;    %��ʼѹǿ(Pa)
D = 0.132;       %�⾶(m)
ep = 0.04;     %�����(m)
r_s = 0.003;        %�Ǽ�Բ���뾶(m)
r1_s = 0.003;        %�Ǹ�����Բ���뾶(m)
l_s = 0.023;        %��������(m)
n_s = 7;        %�ǽ���
beta_s = pi / n_s;        %�ȷֽ�(rad)
epsilon_s = 1;      %�Ƿ���
theta_s = 0.620465;      %�Ǹ����(rad)
Lp = 0.22;      %װҩ����(m)
gamma = 1.2;        %���ȱ�
rho_p = 1730;       %�ܶ�(kg/m^3)
n_p = 0.302;      %ѹǿָ��
c = 1600;       %�����ٶ�(m/s)
At = 5e-4;     %��ܺ����(m^2)
r0 = 5e-3;      %?��ʼȼ��(m/s)
alpha_r = 1.7e-4;        %?ȼ��ϵ��
phi_alpha = 1;       %?��ʴ����
phi_m = 1;      %?

%����ó���
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���
%�ܳ�����
s_a = (sin(epsilon_s*beta_s)) / sin(theta_s) + (1 - epsilon_s)*beta_s;
s_b = (pi / 2 + beta_s - theta_s - cot(theta_s));
s_c = (1 - epsilon_s)*beta_s;
s_d =  (l_s*sin(epsilon_s*beta_s));
s0 = 2*n_s*(l_s*s_a + (r_s + r1_s)*s_b  - beta_s*r1_s);      %ȼ�����ʼ�߳�
%si1 = l_s*s_a + (r_s + r1_s)*s_b - beta_s*(r1_s - e);
%si2 = l_s*s_a + (r_s + e)*s_b;
%si3 = l_s*s_c + (r_s + e)*( beta_s + asin(s_d / (r_s + e)) );
%ͨ���������
Ap_a = ((l_s^2) / 2)*( (1 - epsilon_s)*beta_s + (sin(epsilon_s*beta_s)) ...
    *( (cos(epsilon_s*beta_s)) - (sin(epsilon_s*beta_s))*cot(theta_s) ) );
Ap_b = s_a;
Ap_c = s_b;
Ap_d = cot(theta_s) + theta_s - pi / 2;
Ap_e = s_c;
Ap_f = s_d;
Ap_g = epsilon_s*beta_s;
Api0 = Ap_a + l_s*r_s*Ap_b + ((r_s^2) / 2)*Ap_c + ((r1_s^2) / 2)*Ap_d;
Ap0 = 2*n_s*Api0;
%Api1 = Ap_a + l_s*(r_s + e)*Ap_b + (((r_s + e)^2) / 2)*Ap_c ...
%   + (((r1_s - e)^2) / 2)*Ap_d;
%Api2 = Ap_a + l_s*(r_s + e)*Ap_b + (((r_s + e)^2) / 2)*Ap_c;
%Api3 = ( ((l_s + r_s + e)^2)*Ap_e ...
%   + ((r_s + e)^2)*(Ap_g + asin(Ap_f / (r_s + e))) ...
%   + Ap_f*(sqrt((r_s + e)^2 - Ap_f^2)  + l_s*cos(Ap_g)) ) / 2;

Vp0 = (pi*(D^2) / 4 - Ap0)*Lp;       %ҩ�����(m^3)
mp = Vp0*rho_p;     %ҩ������(kg)


%����
%�����ñ������Ѹ���ֵ��
p = p0*ones(1,long);        %ʵ��ѹǿ(Pa)
r = r0*ones(1,long);     %ȼ��(m/s)
e = 0*ones(1,long);     %����ȥ���(m)
s = s0*ones(1,long);     %ȼ����߳�(m)
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

e_a = r1_s;
e_b = l_s*(sin(epsilon_s*beta_s)) / cos(theta_s) - r_s;

%ѭ������
while (e(i) <= ep)
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
    if(e < e_a)
        s(i) = 2*n_s*(l_s*s_a + (r_s + r1_s)*s_b - beta_s*(r1_s - e(i)));
        Api1 = Ap_a + l_s*(r_s + e(i))*Ap_b + (((r_s + e(i))^2) / 2)*Ap_c ...
            + (((r1_s - e(i))^2) / 2)*Ap_d;
        Ap(i) = 2*n_s*Api1;
    elseif(e < e_b)
        s(i) = 2*n_s*(l_s*s_a + (r_s + e(i))*s_b);
        Api2 = Ap_a + l_s*(r_s + e(i))*Ap_b + (((r_s + e(i))^2) / 2)*Ap_c;
        Ap(i) = 2*n_s*Api2;
    else
        s(i) = 2*n_s*(l_s*s_c + (r_s + e(i))*( beta_s + asin(s_d / (r_s + e(i))) ) );
        Api3 = ( ((l_s + r_s + e(i))^2)*Ap_e ...
            + ((r_s + e(i))^2)*(Ap_g + asin(Ap_f / (r_s + e(i)))) ...
            + Ap_f*(sqrt((r_s + e(i))^2 - Ap_f^2)  + l_s*cos(Ap_g)) ) / 2;
        Ap(i) = 2*n_s*Api3;
    end
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

