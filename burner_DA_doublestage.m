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
long = 150000;        %���鳤��
n = 1:1:long;       %��ͼ������
dt = 2e-5;      %����ֵ
%ȼ���Ҳ���
Dr = 0.082;       %ȼ�����⾶(m)
% At = 1.884785e-3;     %��ܺ����(m^2)
At = 1.103e-3;     %��ܺ����(m^2)
%װҩ����
%�ǿ�װҩ����
Ds = 0.082;       %�⾶(m)
ep_s = 0.018;     %���(m)
r_s = 0.0012;        %�Ǽ�Բ���뾶(m)
r1_s = 0.001;        %�Ǹ�����Բ���뾶(m)
l_s = 0.022;        %��������(m)
n_s = 6;        %�ǽ���
epsilon_s = 0.7;      %�Ƿ���
theta_s = pi*25/180;      %�Ǹ����(rad)
%Բ��װҩ����
Dc = 0.082;     %�⾶(m)
dc = 0.01;      %�ھ�(m)
ep_c = (Dc - dc) / 2;       %���(m)
%�ֶβ���
Lp = 0.49;      %ҩ������(m)
pers = 0.7;     %�ǿ׶�ռ��
Lp_s = Lp*pers;     %�ǿ׶γ���
Lp_c = Lp - Lp_s;       %Բ���γ���
%��֪����
p0 = 1.02e5;    %��ʼѹǿ(Pa)
gamma = 1.2;        %���ȱ�
rho_p = 1700;       %�ܶ�(kg/m^3)
n_p = 0.4;      %ѹǿָ��
c = 1584;       %�����ٶ�(m/s)
rb_0 = 9e-3;      %?��ʼȼ��(m/s)
alpha_r = 9e-5;        %?ȼ��ϵ��
phi_alpha = 1;       %?��ʴ����
phi_m = 1;      %?

%����ó���
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���
beta_s = pi / n_s;        %�ȷֽ�(rad)

%�����м�ֵ
%�ǿ�װҩ
%ȼ���ܳ����в���
s_sa = (sin(epsilon_s*beta_s)) / sin(theta_s) + (1 - epsilon_s)*beta_s;
s_sb = (pi / 2 + beta_s - theta_s - cot(theta_s));
s_sc = (1 - epsilon_s)*beta_s;
s_sd =  (l_s*sin(epsilon_s*beta_s));
%ͨ��������в���
Ap_sa = ((l_s^2) / 2)*( (1 - epsilon_s)*beta_s + (sin(epsilon_s*beta_s)) ...
    *( (cos(epsilon_s*beta_s)) - (sin(epsilon_s*beta_s))*cot(theta_s) ) );
Ap_sb = s_sa;
Ap_sc = s_sb;
Ap_sd = cot(theta_s) + theta_s - pi / 2;
Ap_se = s_sc;
Ap_sf = s_sd;
Ap_sg = epsilon_s*beta_s;
%ѹǿ���в���
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;

%��ʼ����
%�ǿ�װҩ
%�ܳ�����
s_s_0 = 2*n_s*(l_s*s_sa + (r_s + r1_s)*s_sb  - beta_s*r1_s);      
%ͨ���������
Ap_si0 = Ap_sa + l_s*r_s*Ap_sb + ((r_s^2) / 2)*Ap_sc + ((r1_s^2) / 2)*Ap_sd;
Ap_s_0 = 2*n_s*Ap_si0;
%Բ��װҩ
s_c_0 = pi*dc;
Ap_c_0 = pi*(dc^2) / 4;

Vp_s0 = (pi*(Dr^2) / 4 - Ap_s_0)*Lp_s;       %�ǿ�ҩ�����(m^3)
Vp_c0 = (pi*(Dr^2) / 4 - Ap_c_0)*Lp_c;       %Բ��ҩ�����(m^3)
Vp0 = Vp_s0 + Vp_c0;      %��ҩ�����(m^3)
mp_s0 = Vp_s0*rho_p;     %�ǿ�ҩ������(kg)
mp_c0 = Vp_c0*rho_p;     %Բ��ҩ������(kg)
mp = mp_s0 + mp_c0;     %��ҩ������(kg)

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
s_s = s_s_0*ones(1,long);     %�ǿ�ȼ����ʵ�ʱ߳�(m)
s_c = s_c_0*ones(1,long);     %Բ��ȼ����ʵ�ʱ߳�(m)
m_b = rho_p*(s_s_0*Lp_s + s_c_0*Lp_c)*rb_0*ones(1,long);     %ȼ��������(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %��������(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %����(N)
Ab = (s_s_0*Lp_s + s_c_0*Lp_c)*ones(1,long);     %ȼ�����(m^2)
Ap_s = Ap_s_0*ones(1,long);     %�ǿ�ͨ�����(m^2)
Ap_c = Ap_c_0*ones(1,long);     %Բ��ͨ�����(m^2)
Vg = (Ap_s*Lp_s + Ap_c*Lp_c);     %�������(m^3)
%���ֵ����
p_max = p0;        %���ѹǿ(Pa)
e_max = 0;     %�������ȥ���(m)
s_s_max = s_s_0;     %�ǿ�ȼ��������ܳ�(m)
s_c_max = s_c_0;     %Բ��ȼ��������ܳ�(m)
m_b_max = rho_p*(s_s_0*Lp_s + s_c_0*Lp_c)*rb_0;     %���ȼ��������(kg/s)
m_p_max = (phi_m*p0*At / c);     %�����������(kg/s)
F_max = 2000*(p0*At / c);     %�������(N)
Ab_max = (s_s_0*Lp_s + s_c_0*Lp_c);     %���ȼ�����(m^2)
Ap_s_max = Ap_s_0;     %�ǿ����ͨ�����(m^2)
Ap_c_max = Ap_c_0;     %Բ�����ͨ�����(m^2)
Vg_max = (Ap_s*Lp_s + Ap_c*Lp_c);     %����������(m^3)
%��Сֵ����
p_min = p0;        %��Сѹǿ(Pa)
e_min = 0;     %��С����ȥ���(m)
s_s_min = s_s_0;     %ȼ������С�ܳ�(m)
s_c_min = s_c_0;     %ȼ������С�ܳ�(m)
m_b_min = rho_p*(s_s_0*Lp_s + s_c_0*Lp_c)*rb_0;     %��Сȼ��������(kg/s)
m_p_min = (phi_m*p0*At / c);     %��С��������(kg/s)
F_min = 2000*(p0*At / c);     %��С����(N)
Ab_min = (s_s_0*Lp_s + s_c_0*Lp_c);     %��Сȼ�����(m^2)
Ap_s_min = Ap_s_0;     %�ǿ���Сͨ�����(m^2)
Ap_c_min = Ap_c_0;     %Բ����Сͨ�����(m^2)
Vg_min = (Ap_s*Lp_s + Ap_c*Lp_c);     %��С�������(m^3)
%ѭ������
i = 1;
j = 1;
sw = zeros(1,long);        %�׶�ѡ��
swc = zeros(2,100);        %�׶α仯��

%�����м�ֵ

%���׶�ȼ���ܳ�
%     ssi1 = l_s*s_sa + (r_s + r1_s)*s_sb - beta_s*(r1_s - e);
%     ssi2 = l_s*s_sa + (r_s + e)*s_sb;
%     ssi3 = l_s*s_sc + (r_s + e)*( beta_s + asin(s_sd / (r_s + e)) );
%     s_s = 2*n_s*ssi;
%���׶�ͨ�����
%     Apsi1 = Ap_sa + l_s*(r_s + e)*Ap_sb + (((r_s + e)^2) / 2)*Ap_sc ...
%        + (((r1_s - e)^2) / 2)*Ap_sd;
%     Apsi2 = Ap_sa + l_s*(r_s + e)*Ap_sb + (((r_s + e)^2) / 2)*Ap_sc;
%     Apsi3 = ( ((l_s + r_s + e)^2)*Ap_se ...
%        + ((r_s + e)^2)*(Ap_sg + asin(Ap_sf / (r_s + e))) ...
%        + Ap_sf*(sqrt((r_s + e)^2 - Ap_sf^2)  + l_s*cos(Ap_sg)) ) / 2;
%     Ap_s = 2*n_s*Apsi;

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
        swc(1,j) = i;
        swc(2,j) = sw(i);
        j = j + 1;
    end
    %�ֽ׶μ���ȼ���ܳ���ͨ�����
    switch sw(i)
        case{1}     %�ǿ׵�һ�׶Σ�Բ����һ�׶�
            s_s(i) = 2*n_s*(l_s*s_sa + (r_s + r1_s)*s_sb - beta_s*(r1_s - e(i)));
            Apsi1 = Ap_sa + l_s*(r_s + e(i))*Ap_sb + (((r_s + e(i))^2) / 2)*Ap_sc ...
                + (((r1_s - e(i))^2) / 2)*Ap_sd;
            Ap_s(i) = 2*n_s*Apsi1;
            s_c(i) = pi*(dc + 2*e(i));
            Ap_c(i) = pi*(dc + 2*e(i))^2 / 4;
        case{3}     %�ǿ׵ڶ��׶Σ�Բ����һ�׶�
            s_s(i) = 2*n_s*(l_s*s_sa + (r_s + e(i))*s_sb);
            Apsi2 = Ap_sa + l_s*(r_s + e(i))*Ap_sb + (((r_s + e(i))^2) / 2)*Ap_sc;
            Ap_s(i) = 2*n_s*Apsi2;
            s_c(i) = pi*(dc + 2*e(i));
            Ap_c(i) = pi*(dc + 2*e(i))^2 / 4;
        case{7}     %�ǿ׵����׶Σ�Բ����һ�׶�
            s_s(i) = 2*n_s*(l_s*s_sc + (r_s + e(i))*( beta_s + asin(s_sd / (r_s + e(i))) ) );
            Apsi3 = ( ((l_s + r_s + e(i))^2)*Ap_se ...
                + ((r_s + e(i))^2)*(Ap_sg + asin(Ap_sf / (r_s + e(i)))) ...
                + Ap_sf*(sqrt((r_s + e(i))^2 - Ap_sf^2)  + l_s*cos(Ap_sg)) ) / 2;
            Ap_s(i) = 2*n_s*Apsi3;
            s_c(i) = pi*(dc + 2*e(i));
            Ap_c(i) = pi*(dc + 2*e(i))^2 / 4;
        case{15}     %�ǿ׽�����Բ����һ�׶�
            s_s(i) = 0;
            Ap_s(i) = pi*Dr^2 / 4;
            s_c(i) = pi*(dc + 2*e(i));
            Ap_c(i) = pi*(dc + 2*e(i))^2 / 4;
        case{17}     %�ǿ׵�һ�׶Σ�Բ������
            s_s(i) = 2*n_s*(l_s*s_sa + (r_s + r1_s)*s_sb - beta_s*(r1_s - e(i)));
            Apsi1 = Ap_sa + l_s*(r_s + e(i))*Ap_sb + (((r_s + e(i))^2) / 2)*Ap_sc ...
                + (((r1_s - e(i))^2) / 2)*Ap_sd;
            Ap_s(i) = 2*n_s*Apsi1;
            s_c(i) = 0;
            Ap_c(i) = pi*Dr^2 / 4;
        case{19}     %�ǿ׵ڶ��׶Σ�Բ������
            s_s(i) = 2*n_s*(l_s*s_sa + (r_s + e(i))*s_sb);
            Apsi2 = Ap_sa + l_s*(r_s + e(i))*Ap_sb + (((r_s + e(i))^2) / 2)*Ap_sc;
            Ap_s(i) = 2*n_s*Apsi2;
            s_c(i) = 0;
            Ap_c(i) = pi*Dr^2 / 4;
        case{23}     %�ǿ׵����׶Σ�Բ������
            s_s(i) = 2*n_s*(l_s*s_sc + (r_s + e(i))*( beta_s + asin(s_sd / (r_s + e(i))) ) );
            Apsi3 = ( ((l_s + r_s + e(i))^2)*Ap_se ...
                + ((r_s + e(i))^2)*(Ap_sg + asin(Ap_sf / (r_s + e(i)))) ...
                + Ap_sf*(sqrt((r_s + e(i))^2 - Ap_sf^2)  + l_s*cos(Ap_sg)) ) / 2;
            Ap_s(i) = 2*n_s*Apsi3;
            s_c(i) = 0;
            Ap_c(i) = pi*Dr^2 / 4;
        otherwise     %�����׶�
            s_s(i) = s_s(i - 1);
            Ap_s(i) = Ap_s(i - 1);
            s_c(i) = s_c(i - 1);
            Ap_c(i) = Ap_c(i - 1);
    end
    %������������
    Ab(i) = s_s(i)*Lp_s + s_c(i)*Lp_c;
    Vg(i) = Ap_s(i)*Lp_s + Ap_c(i)*Lp_c;
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
    if(s_s(i) > s_s_max)
        s_s_max = s_s(i);
    end
    if(s_c(i) > s_c_max)
        s_c_max = s_c(i);
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
    if(Ap_s(i) > Ap_s_max)
        Ap_s_max = Ap_s(i);
    end
    if(Ap_c(i) > Ap_c_max)
        Ap_c_max = Ap_c(i);
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
    if(s_s(i) < s_s_min)
        s_s_min = s_s(i);
    end
    if(s_c(i) < s_c_min)
        s_c_min = s_c(i);
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
    if(Ap_s(i) < Ap_s_min)
        Ap_s_min = Ap_s(i);
    end
    if(Ap_c(i) < Ap_c_min)
        Ap_c_min = Ap_c(i);
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
s_s(i:1:long) = [];
s_c(i:1:long) = [];
m_b(i:1:long) = [];
m_p(i:1:long) = [];
F(i:1:long) = [];
Ab(i:1:long) = [];
Ap_s(i:1:long) = [];
Ap_c(i:1:long) = [];
Vg(i:1:long) = [];
sw(i:1:long) = [];
swc(:,j:100) = [];
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


subplot(2,2,3);
plot(t,Ap_s);
axis ([pri*t_max,prx*t_max,(pri*(Ap_s_max - Ap_s_min) + Ap_s_min), ...
    (prx*(Ap_s_max - Ap_s_min) + Ap_s_min)]);
title('�ǿ�ͨ�����');
xlabel('ʱ��(s)');
ylabel('���(m^2)');
subplot(2,2,4);
plot(t,Ap_c);
axis ([pri*t_max,prx*t_max,(pri*(Ap_c_max - Ap_c_min) + Ap_c_min), ...
    (prx*(Ap_c_max - Ap_c_min) + Ap_c_min)]);
title('Բ��ͨ�����');
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
% legend('����2');

figure;
hold on;
plot(t,p);
axis ([pri*t_max,prx*t_max,(pri*(p_max - p_min) + p_min), ...
    (prx*(p_max - p_min) + p_min)]);
title('ȼ����ѹ��');
xlabel('ʱ��(s)');
ylabel('ѹǿ(Pa)');
% legend('����2');


%˫Y��ͼ��
%��ͼ����������y������
left_color = [0 0 0];
right_color = [0 0 0];
set(figure,'defaultAxesColorOrder',[left_color;right_color]);
hold on;
%�������
yyaxis left;
plot(t,p,'-','LineWidth',1,'color','k');   
ylabel('ѹǿ(Pa)')
%���ÿ̶�
per = 0.6;      %������
cou = p_max;        %ѡ����
mag = floor(log10(cou / per));        %��������
fir = floor((cou / per) / (10^mag));      %ȡ��һλ
sec = floor((cou / per - fir*10^mag) / (10^(mag - 1)));       %ȡ�ڶ�λ
ran = fir*10^mag + sec*10^(mag - 1);        %������
cal = fir*10^(mag - 1);     %���̶�
sta = -0.5*fir*10^(mag - 1);        %�����
axis([pri*t_max,prx*t_max,sta,ran]);
set(gca,'YTick',0:cal:ran);
%�����Ҳ�
yyaxis right;
plot(t,F,'-.','LineWidth',1,'color','k');
ylabel('��(N)');
%���ÿ̶�
per = 0.8;      %������
cou = F_max;        %ѡ����
mag = floor(log10(cou / per));        %��������
fir = floor((cou / per) / (10^mag));      %ȡ��һλ
sec = floor((cou / per - fir*10^mag) / (10^(mag - 1)));       %ȡ�ڶ�λ
ran = fir*10^mag + sec*10^(mag - 1);        %������
cal = fir*10^(mag - 1);     %���̶�
sta = -0.5*fir*10^(mag - 1);        %�����
axis([pri*t_max,prx*t_max,sta,ran]);
set(gca,'YTick',0:cal:ran);
%����X��ͱ���
xlabel('ʱ��(s)');
title('ѹǿ������');
legend('ѹǿ', '����');

%����
