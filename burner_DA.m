%181016
%Dashaway
%ȼ�ղ�������


%190219
%����ʽ
%����ע��

%190321
%�ֶ�ȼ��


clear;
close all;
%��ʼֵ�趨


%����

%���㸨������
long = 100000;        %���鳤��
n = 1:1:long;       %��ͼ������
dt = 2e-5;      %����ֵ
%ȼ���Ҳ���
Dr = 0.132;       %ȼ�����⾶(m)
% At = 1.884785e-3;     %��ܺ����(m^2)
At = 5e-4;     %��ܺ����(m^2)
%װҩ����
%�ǿ�װҩ����
Ds = 0.132;       %�⾶(m)
ep_s = 0.04;     %���(m)
r_s = 0.003;        %�Ǽ�Բ���뾶(m)
r1_s = 0.003;        %�Ǹ�����Բ���뾶(m)
l_s = 0.023;        %��������(m)
n_s = 7;        %�ǽ���
epsilon_s = 1;      %�Ƿ���
theta_s = 0.620465;      %�Ǹ����(rad)
%Բ��װҩ����
Dc = 0.132;     %�⾶(m)
rc = 0.05;      %�ھ�(m)
ep_c = (Dc - rc) / 2;       %���(m)
%�ֶβ���
Lp = 0.22;      %ҩ������(m)
pers = 0.5;     %�ǿ׶�ռ��
Lp_s = Lp*pers;     %�ǿ׶γ���
Lp_c = Lp - Lp_s;       %Բ���γ���
%��֪����
p0 = 1.02e5;    %��ʼѹǿ(Pa)
gamma = 1.2;        %���ȱ�
rho_p = 1730;       %�ܶ�(kg/m^3)
n_p = 0.302;      %ѹǿָ��
c = 1600;       %�����ٶ�(m/s)
rb_0 = 5e-3;      %?��ʼȼ��(m/s)
alpha_r = 1.7e-4;        %?ȼ��ϵ��
phi_alpha = 1;       %?��ʴ����
phi_m = 1;      %?

%����ó���
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���
beta_s = pi / n_s;        %�ȷֽ�(rad)

%�����м�ֵ
%�ǿ�װҩ
%ȼ���ܳ����в���
s_a = (sin(epsilon_s*beta_s)) / sin(theta_s) + (1 - epsilon_s)*beta_s;
s_b = (pi / 2 + beta_s - theta_s - cot(theta_s));
s_c = (1 - epsilon_s)*beta_s;
s_d =  (l_s*sin(epsilon_s*beta_s));
%ͨ��������в���
Ap_a = ((l_s^2) / 2)*( (1 - epsilon_s)*beta_s + (sin(epsilon_s*beta_s)) ...
    *( (cos(epsilon_s*beta_s)) - (sin(epsilon_s*beta_s))*cot(theta_s) ) );
Ap_b = s_a;
Ap_c = s_b;
Ap_d = cot(theta_s) + theta_s - pi / 2;
Ap_e = s_c;
Ap_f = s_d;
Ap_g = epsilon_s*beta_s;
%ѹǿ���в���
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;

%��ʼ����
%�ǿ�װҩ
%�ܳ�����
s_s_a0 = 2*n_s*(l_s*s_a + (r_s + r1_s)*s_b  - beta_s*r1_s);      
%ͨ���������
Ap_si0 = Ap_a + l_s*r_s*Ap_b + ((r_s^2) / 2)*Ap_c + ((r1_s^2) / 2)*Ap_d;
Ap_s_a0 = 2*n_s*Ap_si0;
%Բ��װҩ
s_c_0 = 2*pi*rc;
Ap_c_0 = pi*(rc^2);


%Լ������
%����Լ������

%�ֽ׶��б�����
e_a = r1_s;     %�Ǹ������ʧ
e_b = l_s*(sin(epsilon_s*beta_s)) / cos(theta_s) - r_s;     %�Ǳ���ʧ

%���׶ζ�Ӧ����
%     �ǿ׵�һ�׶Σ�Բ����һ�׶�
%     e11 = ((e <= e_a) & (e <= e_b) & (e <= ep_s) & (e <= ep_c));
%     sw = 1;
%     �ǿ׵ڶ��׶Σ�Բ����һ�׶�
%     e21 = ((e > e_a) & (e < e_b) & (e <= ep_s) & (e <= ep_c));
%     sw = 3;
%     �ǿ׵����׶Σ�Բ����һ�׶�
%     e31 = ((e > e_a) & (e > e_b) & (e <= ep_s) & (e <= ep_c));
%     sw = 7;
%     �ǿ׽�����Բ����һ�׶�
%     e41 = ((e > e_a) & (e > e_b) & (e > ep_s) & (e <= ep_c));
%     sw = 15;
%     �ǿ׵�һ�׶Σ�Բ������
%     e12 = ((e <= e_a) & (e <= e_b) & (e <= ep_s) & (e > ep_c));
%     sw = 17;
%     �ǿ׵ڶ��׶Σ�Բ������
%     e22 = ((e > e_a) & (e < e_b) & (e <= ep_s) & (e > ep_c));
%     sw = 19;
%     �ǿ׵����׶Σ�Բ������
%     e32 = ((e > e_a) & (e > e_b) & (e <= ep_s) & (e > ep_c));
%     sw = 23;
%�����������
ep = max(ep_s,ep_c);
%     e >= ep;


%����
%�����ñ������飨�Ѹ���ֵ��
p = p0*ones(1,long);        %ʵ��ѹǿ(Pa)
rb = rb_0*ones(1,long);     %ȼ��(m/s)
e = 0*ones(1,long);     %����ȥ���(m)
s_s = s_s_a0*ones(1,long);     %�ǿ�ȼ����ʵ�ʱ߳�(m)
s_c = s_c_0*ones(1,long);     %Բ��ȼ����ʵ�ʱ߳�(m)
m_b = rho_p*s_a0*Lp*rb_0*ones(1,long);     %ȼ��������(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %��������(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %����(N)
Ab = s_a0*Lp*ones(1,long);     %ȼ�����(m^2)
Ap_s = Ap_s_a0*ones(1,long);     %�ǿ�ͨ�����(m^2)
Ap_c = Ap_c_0*ones(1,long);     %Բ��ͨ�����(m^2)
Vg = Ap*Lp;     %�������(m^3)
%���ֵ����
p_max = p0;        %���ѹǿ(Pa)
e_max = 0;     %�������ȥ���(m)
s_max = s_a0;     %ȼ��������ܳ�(m)
m_b_max = rho_p*s_a0*Lp*rb_0;     %���ȼ��������(kg/s)
m_p_max = (phi_m*p0*At / c);     %�����������(kg/s)
F_max = 2000*(p0*At / c);     %�������(?)
Ab_max = s_a0*Lp;     %���ȼ�����(m^2)
Ap_max = Ap_a0;     %���ͨ�����(m^2)
Vg_max = Ap_max*Lp;     %����������(m^3)
%��Сֵ����
p_min = p0;        %��Сѹǿ(Pa)
e_min = 0;     %��С����ȥ���(m)
s_min = s_a0;     %ȼ������С�ܳ�(m)
m_b_min = rho_p*s_a0*Lp*rb_0;     %��Сȼ��������(kg/s)
m_p_min = (phi_m*p0*At / c);     %��С��������(kg/s)
F_min = 2000*(p0*At / c);     %��С����(?)
Ab_min = s_a0*Lp;     %��Сȼ�����(m^2)
Ap_min = Ap_a0;     %��Сͨ�����(m^2)
Vg_min = Ap_min*Lp;     %��С�������(m^3)
%ѭ������
i = 1;
j = 1;
sw = 0*ones(1,long);        %�׶�ѡ��
swc = 0*ones(1,100);        %�׶α仯��

%�����м�ֵ

%���׶�ȼ���ܳ�
%     si1 = l_s*s_a + (r_s + r1_s)*s_b - beta_s*(r1_s - e);
%     si2 = l_s*s_a + (r_s + e)*s_b;
%     si3 = l_s*s_c + (r_s + e)*( beta_s + asin(s_d / (r_s + e)) );
%     s = 2*n_s*si;
%���׶�ͨ�����
%     Api1 = Ap_a + l_s*(r_s + e)*Ap_b + (((r_s + e)^2) / 2)*Ap_c ...
%        + (((r1_s - e)^2) / 2)*Ap_d;
%     Api2 = Ap_a + l_s*(r_s + e)*Ap_b + (((r_s + e)^2) / 2)*Ap_c;
%     Api3 = ( ((l_s + r_s + e)^2)*Ap_e ...
%        + ((r_s + e)^2)*(Ap_g + asin(Ap_f / (r_s + e))) ...
%        + Ap_f*(sqrt((r_s + e)^2 - Ap_f^2)  + l_s*cos(Ap_g)) ) / 2;
%     Ap = 2*n_s*Api;

%����������㹫ʽ
%     dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;

%��ֵ����
%ѭ��ǰ����������



%ѭ������
while (e(i) <= ep)
    %��������𲽼��㣬������µ�ѹǿֵ
    %dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;
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

    %�жϴ�ʱ�����ĸ��׶�
    sw(i) = 1;
    if(e(i) <= e_a)
        sw(i) = sw(i) + 0*2;
    else
        sw(i) = sw(i) + 1*2;
    end
    if(e(i) <= e_b)
        sw(i) = sw(i) + 0*4;
    else
        sw(i) = sw(i) + 1*4;
    end
    if(e(i) <= ep_s)
        sw(i) = sw(i) + 0*8;
    else
        sw(i) = sw(i) + 1*8;
    end
    if(e(i) <= ep_c)
        sw(i) = sw(i) + 0*16;
    else
        sw(i) = sw(i) + 1*16;
    end
    %��¼�׶θı�ʱ��ѭ����ֵ
    if(sw(i) ~= sw(i - 1))
        swc(j) = i;
        j = j + 1;
    end
    %�ֽ׶μ���ȼ���ܳ���ͨ�����
    switch sw(i)
        case{1}     %�ǿ׵�һ�׶Σ�Բ����һ�׶�
            s(i) = 2*n_s*(l_s*s_a + (r_s + r1_s)*s_b - beta_s*(r1_s - e(i)));
            Api1 = Ap_a + l_s*(r_s + e(i))*Ap_b + (((r_s + e(i))^2) / 2)*Ap_c ...
                + (((r1_s - e(i))^2) / 2)*Ap_d;
            Ap(i) = 2*n_s*Api1;
        case{3}     %�ǿ׵ڶ��׶Σ�Բ����һ�׶�
            s(i) = 2*n_s*(l_s*s_a + (r_s + e(i))*s_b);
            Api2 = Ap_a + l_s*(r_s + e(i))*Ap_b + (((r_s + e(i))^2) / 2)*Ap_c;
            Ap(i) = 2*n_s*Api2;
        case{7}     %�ǿ׵����׶Σ�Բ����һ�׶�
            s(i) = 2*n_s*(l_s*s_c + (r_s + e(i))*( beta_s + asin(s_d / (r_s + e(i))) ) );
            Api3 = ( ((l_s + r_s + e(i))^2)*Ap_e ...
                + ((r_s + e(i))^2)*(Ap_g + asin(Ap_f / (r_s + e(i)))) ...
                + Ap_f*(sqrt((r_s + e(i))^2 - Ap_f^2)  + l_s*cos(Ap_g)) ) / 2;
            Ap(i) = 2*n_s*Api3;
        case{15}     %�ǿ׽�����Բ����һ�׶�
            s(i) = 2*n_s*(l_s*s_c + (r_s + e(i))*( beta_s + asin(s_d / (r_s + e(i))) ) );
            Api3 = ( ((l_s + r_s + e(i))^2)*Ap_e ...
                + ((r_s + e(i))^2)*(Ap_g + asin(Ap_f / (r_s + e(i)))) ...
                + Ap_f*(sqrt((r_s + e(i))^2 - Ap_f^2)  + l_s*cos(Ap_g)) ) / 2;
            Ap(i) = 2*n_s*Api3;
        case{17}     %�ǿ׵�һ�׶Σ�Բ������
            s(i) = 2*n_s*(l_s*s_a + (r_s + r1_s)*s_b - beta_s*(r1_s - e(i)));
            Api1 = Ap_a + l_s*(r_s + e(i))*Ap_b + (((r_s + e(i))^2) / 2)*Ap_c ...
                + (((r1_s - e(i))^2) / 2)*Ap_d;
            Ap(i) = 2*n_s*Api1;
        case{19}     %�ǿ׵ڶ��׶Σ�Բ������
            s(i) = 2*n_s*(l_s*s_a + (r_s + e(i))*s_b);
            Api2 = Ap_a + l_s*(r_s + e(i))*Ap_b + (((r_s + e(i))^2) / 2)*Ap_c;
            Ap(i) = 2*n_s*Api2;
        case{23}     %�ǿ׵����׶Σ�Բ������
            s(i) = 2*n_s*(l_s*s_c + (r_s + e(i))*( beta_s + asin(s_d / (r_s + e(i))) ) );
            Api3 = ( ((l_s + r_s + e(i))^2)*Ap_e ...
                + ((r_s + e(i))^2)*(Ap_g + asin(Ap_f / (r_s + e(i)))) ...
                + Ap_f*(sqrt((r_s + e(i))^2 - Ap_f^2)  + l_s*cos(Ap_g)) ) / 2;
            Ap(i) = 2*n_s*Api3;
        otherwise     %�����׶�
            s(i) = s(i - 1);
            Ap(i) = Ap(i - 1);
    end
    %������������
    Ab(i) = s(i)*Lp;
    Vg(i) = Ap(i)*Lp;
    m_b(i) = rho_p*Ab(i)*rb(i);
    m_p(i) = phi_m*p(i)*At / c;
    F(i) = 2000*m_p(i) + 4.4*At*(p(i)-p0);
    
    %��¼���ֵ
    if(p(i) > p_max)
        p_max = p(i);
    end
    if(e(i) > e_max)
        e_max = e(i);
    end
    if(s(i) > s_max)
        s_max = s(i);
    end
    if(m_b(i) > m_b_max)
        m_b_max = m_b(i);
    end
    if(m_p(i) > m_p_max)
        m_p_max = m_p(i);
    end
    if(F(i) > F_max)
        F_max = F(i);
    end
    if(Ab(i) > Ab_max)
        Ab_max = Ab(i);
    end
    if(Ap(i) > Ap_max)
        Ap_max = Ap(i);
    end
    if(Vg(i) > Vg_max)
        Vg_max = Vg(i);
    end
    
    %��¼��Сֵ
    if(p(i) < p_min)
        p_min = p(i);
    end
    if(e(i) < e_min)
        e_min = e(i);
    end
    if(s(i) < s_min)
        s_min = s(i);
    end
    if(m_b(i) < m_b_min)
        m_b_min = m_b(i);
    end
    if(m_p(i) < m_p_min)
        m_p_min = m_p(i);
    end
    if(F(i) < F_min)
        F_min = F(i);
    end
    if(Ab(i) < Ab_min)
        Ab_min = Ab(i);
    end
    if(Ap(i) < Ap_min)
        Ap_min = Ap(i);
    end
    if(Vg(i) < Vg_min)
        Vg_min = Vg(i);
    end
    
    %���鳤��Լ��
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
sw(i:1:long) = [];
t = dt*n;
t_max = dt*(i - 1);     %ȼ����ʱ��
prx = 1.05;
pri = -0.05;

str = [' ȼ����ʱ�䣺 ',num2str(t_max),'s'];
disp(str);

%�������
%д���ĵ�


%����ͼ��
figure;
hold on;
subplot(2,1,1);
plot(t,p);
axis ([pri*t_max,prx*t_max,(pri*(p_max - p_min) + p_min), ...
    (prx*(p_max - p_min) + p_min)]);
title('ȼ����ѹ��');
xlabel('ʱ��(s)');
ylabel('ѹǿ(Pa)');
% legend('����2');

subplot(2,1,2);
plot(t,e);
axis ([pri*t_max,prx*t_max,(pri*(e_max - e_min) + e_min), ...
    (prx*(e_max - e_min) + e_min)]);
title('��ȥ���');
xlabel('ʱ��(s)');
ylabel('����(m)');
% legend('����2');

figure;
hold on;
subplot(2,1,1);
plot(t,Ab);
axis ([pri*t_max,prx*t_max,(pri*(Ab_max - Ab_min) + Ab_min), ...
    (prx*(Ab_max - Ab_min) + Ab_min)]);
title('ȼ�����');
xlabel('ʱ��(s)');
ylabel('���(m^2)');
% legend('����2');


subplot(2,1,2);
plot(t,Ap);
axis ([pri*t_max,prx*t_max,(pri*(Ap_max - Ap_min) + Ap_min), ...
    (prx*(Ap_max - Ap_min) + Ap_min)]);
title('ͨ�����');
xlabel('ʱ��(s)');
ylabel('���(m^2)');
% legend('����2');

figure;
hold on;
subplot(2,1,1);
plot(t,m_b);
axis ([pri*t_max,prx*t_max,(pri*(m_b_max - m_b_min) + m_b_min), ...
    (prx*(m_b_max - m_b_min) + m_b_min)]);
title('ȼ��������');
xlabel('ʱ��(s)');
ylabel('��������(kg/s)');
% legend('����2');

subplot(2,1,2);
plot(t,m_p);
axis ([pri*t_max,prx*t_max,(pri*(m_p_max - m_p_min) + m_p_min), ...
    (prx*(m_p_max - m_p_min) + m_p_min)]);
title('��������');
xlabel('ʱ��(s)');
ylabel('��������(kg/s)');
% legend('����2');

figure;
hold on;
plot(t,F);
axis ([pri*t_max,prx*t_max,(pri*(F_max - F_min) + F_min), ...
    (prx*(F_max - F_min) + F_min)]);
title('����');
xlabel('ʱ��(s)');
ylabel('��(N)');
legend('����2');

figure;
hold on;
plot(t,p);
axis ([pri*t_max,prx*t_max,(pri*(p_max - p_min) + p_min), ...
    (prx*(p_max - p_min) + p_min)]);
title('ȼ����ѹ��');
xlabel('ʱ��(s)');
ylabel('ѹǿ(Pa)');
legend('����2');


%����
