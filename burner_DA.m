%181016
%Dashaway
%ȼ�ղ�������


%190117
%����ע��
%��ΪԲ��ҩ����ҩ����������

%190212
%��Ϊ�ܵ����ڿ�
%���Ĳ���

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
Dr = 0.150;       %ȼ�����⾶(m)
r = 0.00625;        %΢Բ��(m)
l = 0.04375;        %���ľ�(m)
R = 0.05;       %���Ŀ�(m)

n_s = 6;    %ҩ������
Lp = 0.22;      %װҩ����(m)
gamma = 1.2;        %���ȱ�
rho_p = 1730;       %�ܶ�(kg/m^3)
n_p = 0.302;      %ѹǿָ��
c = 1600;       %�����ٶ�(m/s)
At = 1.884785e-3;     %��ܺ����(m^2)
rb_0 = 5e-3;      %?��ʼȼ��(m/s)
alpha_r = 1.7e-4;        %?ȼ��ϵ��
phi_alpha = 1;       %?��ʴ����
phi_m = 1;      %?

%����ó���

Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���
%�б�����
e_a = l*sqrt((1 - cos(pi / n_s)) / 2) - r;       %��һ�׶�
e_b = Dr/2 - r - l;
%�ܳ�����
s_a0 = 2*pi*(R + l + n_s*r) + 2*pi*n_s*(r);      %��һ�׶�ȼ�����ʼ�߳�
    %s_a =  2*pi*(R + l + e) + 2*pi*n_s*(r + e);    %��һ�׶�ȼ���ܳ���ʽ
    %s_b = 4*n_s*(r + e)*asin( (l*sin(pi / (2*n_s))) / (r + e) ) ...   %�ڶ��׶�ȼ���ܳ���ʽ
        %+ 2*pi*(l + R + e)
%ͨ���������
Ap_a0 = 2*pi*l*r + n_s*pi*r^2 + pi*R^2;     %��һ�׶γ�ʼͨ�����
    %Ap_a = 2*pi*l*(r + e) + n_s*pi*(r + e)^2 + pi*(R + e)^2;   %��һ�׶�ͨ�����
    %Ap_b = 2*n_s*((r + e)^2)*asin( (l*sin(pi / (2*n_s))) / (r + e) ) ...  %��һ�׶�ͨ�����
        %+ 2*pi*l*(r + e) + pi*(R + e)^2;
        

%����
%�����ñ������Ѹ���ֵ��
p = p0*ones(1,long);        %ʵ��ѹǿ(Pa)
rb = rb_0*ones(1,long);     %ȼ��(m/s)
e = 0*ones(1,long);     %����ȥ���(m)
s = s_a0*ones(1,long);     %ȼ����ʵ�ʱ߳�(m)
m_b = rho_p*s_a0*Lp*rb_0*ones(1,long);     %ȼ��������(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %��������(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %����(?)
Ab = s_a0*Lp*ones(1,long);     %ȼ�����(m^2)
Ap = Ap_a0*ones(1,long);     %ͨ�����(m^2)
Vg = Ap*Lp;     %�������(m^3)
%ѭ������
i = 1;


%��ֵ����
%ѭ��ǰ����������
%ѹǿ���в���
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;




%ѭ������
while e <= e_b
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
    rb(i) = alpha_r*p(i)^n_p;
    e(i) = e(i - 1) + rb(i)*dt;
    %�ֽ׶μ���ȼ���ܳ���ͨ�����
    if(e < e_a)
        s(i) = 2*pi*(R + l + e(i)) + 2*pi*n_s*(r + e(i));
        Ap(i) = 2*pi*l*(r + e(i)) + n_s*pi*(r + e(i))^2 + pi*(R + e(i))^2;
    elseif(e < e_b)
        tr =  asin( (l*sin(pi / (2*n_s))) / (r + e(i)) ) ;
        s(i) = 4*n_s*(r + e(i))*asin( (l*sin(pi / (2*n_s))) / (r + e(i)) ) ...
                + 2*pi*(l + R + e(i));
        Ap(i) = 2*n_s*((r + e(i))^2)*asin( (l*sin(pi / (2*n_s))) / (r + e(i)) ) ...
                + 2*pi*l*(r + e(i)) + pi*(R + e(i))^2;
    else
        s(i) = 2*pi*(R + l + e(i)) + 2*pi*n_s*(r + e(i));
        Ap(i) = 2*pi*l*(r + e(i)) + n_s*pi*(r + e(i))^2 + pi*(R + e(i))^2;
    end
    
    Ab(i) = s(i)*Lp;
    Vg(i) = Ap(i)*Lp;
    m_b(i) = rho_p*Ab(i)*rb(i);
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
rb(i:1:long) = [];
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

