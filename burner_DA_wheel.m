%181016
%Dashaway
%ȼ�ղ�������


%190327
%��ʽװҩ


clear;
close all;
%��ʼֵ�趨


%����

%���㸨������
long = 150000;        %���鳤��
n = 1:1:long;       %��ͼ������
dt = 5e-5;      %����ֵ
%ȼ���Ҳ���
Dr = 0.240;       %ȼ�����⾶(m)
% At = 1.884785e-3;     %��ܺ����(m^2)
At = 7.853982e-3;     %��ܺ����(m^2)
%װҩ����
%��ʽװҩ����
Dc = 0.238;     %�⾶(m)
dc = 0.0868;      %�ھ�(m)
e_s1 = 0.0155;       %�����(m)
e_s2 = 0.0033;       %�����(m)
n_s = 12;       %�ֿ���
l_s1 = 0.0995;      %��Բ����������(m)
l_s2 = 0.053;       %��Բ����������(m)
r_s1 = 0.004;       %�Ϲ���Բ���뾶(m)
r_s2 = 0.003;       %�¹���Բ���뾶(m)
Lp = 0.345;      %ҩ������(m)
%�ƽ�������
n_p = 0.302;      %ѹǿָ��
rho_p = 1730;       %�ܶ�(kg/m^3)
c = 1600;       %�����ٶ�(m/s)
rb_0 = 5e-3;      %?��ʼȼ��(m/s)
alpha_r = 1.7e-4;        %?ȼ��ϵ��(m/s*MPa)

%��֪����
p0 = 1.02e5;    %��ʼѹǿ(Pa)
gamma = 1.2;        %���ȱ�
phi_alpha = 1;       %?��ʴ����
phi_m = 1;      %?

%����ó���
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���
beta_s = 2*pi / n_s;    %�ֿ׽�(rad)

%�����м�ֵ
as1 = asin((e_s1 + r_s1) / l_s1);       %���ַ��ϲ��Ƿ���
as2 = asin((e_s1 + r_s2) / l_s2);       %���ַ��²��Ƿ���
as3 = asin((e_s2 + r_s1) / l_s1);       %���ַ��ϲ��Ƿ���
as4 = asin((e_s2 + r_s2) / l_s2);       %���ַ��²��Ƿ���
h_s1 = l_s1*cos(as1) - l_s2*cos(as2);       %���ַ��߶�
h_s2 = l_s1*cos(as3) - l_s2*cos(as4);       %���ַ��߶�
%ѹǿ���в���
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;
%���׶�ȼ���ܳ�
% si1 = (l_s1 + r_s1 + e)*(beta_s - as1 - as3) + (r_s1 + e)*(pi + as1 + as2) ...
%     + h_s1 + h_s2 + (r_s2 + e)*(pi - as2 - as4) ...
%     + (l_s2 - r_s2 - e)*(beta_s - as2 - as4) + (e + dc/2)*beta_s;
% si2 = (l_s1 + r_s1 + e)*(beta_s - as1 - as3) + (r_s1 + e)*(pi/2 + as1) ...
%     + l_s1*cos(as1) - (e + dc/2) ...
%     + (e + dc/2)*asin((l_s2*sin(as2) - r_s2 - e) / (e + dc/2)) ...
%     + (r_s1 + e)*(as3 + asin( l_s1*sin(as3) / (r_s1 + e) ));
% s = n_s*si*Lp;
%���׶�ͨ�����
% Api1 = (l_s1^2)*((sin(2*as1) + sin(2*as2)) / 4) ...
%     - (l_s2^2)*((sin(2*as2) + sin(2*as4)) / 4) ...
%     + ((l_s1 + r_s1 + e)^2)*((beta_s - as1 - as3) / 2) ...
%     - ((l_s2 - r_s2 - e)^2)*((beta_s - as2 - as4) / 2) ...
%     + ((r_s1 + e)^2)*((pi + as1 + as3) / 2) ...
%     + ((r_s2 + e)^2)*((pi - as2 - as4) / 2) ...
%     - h_s1*(e_s1 - e) - h_s2*(e_s2 - e) + ((e + dc/2)^2)*beta_s/2;
% Api2 = ((l_s1 + r_s1 + e)^2)*((beta_s - as1 - as3) / 2) ...
%     + (((r_s1 + e)^2) / 2)*(as3 + asin( l_s1*sin(as3) / (r_s1 + e) )) ...
%     + (l_s1*sin(as3) / 2)*( sqrt( (r_s1 + e)^2 - (l_s1*sin(as3))^2 ) + l_s1*cos(as3) ) ...
%     + (((r_s1 + e)^2) / 2)*(pi/2 + as1) + ((l_s1^2) / 2)*sin(as1)*cos(as1) ...
%     - (l_s1*sin(as1) - r_s1 - e)*(l_s1*cos(as1) - e - dc/2);
% Ap = n_s*Api;


%��ʼ����
%Բ��װҩ
%�ܳ�����
si0 = (l_s1 + r_s1)*(beta_s - as1 - as3) + r_s1*(pi + as1 + as2) ...
    + h_s1 + h_s2 + r_s2*(pi - as2 - as4) ...
    + (l_s2 - r_s2)*(beta_s - as2 - as4) + (dc/2)*beta_s;
s_0 = n_s*si0*Lp;
%ͨ���������
Api0 = (l_s1^2)*((sin(2*as1) + sin(2*as2)) / 4) ...
    - (l_s2^2)*((sin(2*as2) + sin(2*as4)) / 4) ...
    + ((l_s1 + r_s1)^2)*((beta_s - as1 - as3) / 2) ...
    - ((l_s2 - r_s2)^2)*((beta_s - as2 - as4) / 2) ...
    + (r_s1^2)*((pi + as1 + as3) / 2) ...
    + (r_s2^2)*((pi - as2 - as4) / 2) ...
    - h_s1*e_s1 - h_s2*e_s2 + ((dc/2)^2)*beta_s/2;
Ap_0 = n_s*Api0;

Vp0 = (pi*(Dr^2) / 4 - Ap_0)*Lp;       %ҩ�����(m^3)
mp = Vp0*rho_p;     %ҩ������(kg)
%ѹǿ��


%Լ������
%����Լ������

%�ֽ׶��б�����


%���׶ζ�Ӧ����
%     �������׶�
%     e1 = ((e <= e_s1) & (e <= e_s2));
%     sw = 1;
%     С�����׶�
%     e2 = ((e <= e_s1) & (e > e_s2));
%     sw = 3;
%�����������
ep = e_s1;
%     e >= ep;


%����
%�����ñ������飨�Ѹ���ֵ��
p = p0*ones(1,long);        %ʵ��ѹǿ(Pa)
rb = rb_0*ones(1,long);     %ȼ��(m/s)
e = 0*ones(1,long);     %����ȥ���(m)
s = s_0*ones(1,long);     %ȼ����ʵ���ܳ�(m)
m_b = rho_p*s_0*Lp*rb_0*ones(1,long);     %ȼ��������(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %��������(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %����(N)
Ab = s_0*Lp*ones(1,long);     %ȼ�����(m^2)
Ap = Ap_0*ones(1,long);     %ͨ�����(m^2)
Vg = Ap*Lp;     %�������(m^3)
%���ֵ����
p_max = p0;        %���ѹǿ(Pa)
e_max = 0;     %�������ȥ���(m)
s_max = s_0;     %ȼ��������ܳ�(m)
m_b_max = rho_p*s_0*Lp*rb_0;     %���ȼ��������(kg/s)
m_p_max = (phi_m*p0*At / c);     %�����������(kg/s)
F_max = 2000*(p0*At / c);     %�������(N)
Ab_max = s_0*Lp;     %���ȼ�����(m^2)
Ap_max = Ap_0;     %���ͨ�����(m^2)
Vg_max = Ap_0*Lp;     %����������(m^3)
%��Сֵ����
p_min = p0;        %��Сѹǿ(Pa)
e_min = 0;     %��С����ȥ���(m)
s_min = s_0;     %ȼ������С�ܳ�(m)
m_b_min = rho_p*s_0*Lp*rb_0;     %��Сȼ��������(kg/s)
m_p_min = (phi_m*p0*At / c);     %��С��������(kg/s)
F_min = 2000*(p0*At / c);     %��С����(N)
Ab_min = s_0*Lp;     %��Сȼ�����(m^2)
Ap_min = Ap_0;     %Բ����Сͨ�����(m^2)
Vg_min = Ap_0*Lp;     %��С�������(m^3)
%ѭ������
i = 1;
j = 1;
sw = zeros(1,long);        %�׶�ѡ��
swc = zeros(2,100);        %�׶α仯��

%�����м�ֵ

%���׶�ȼ���ܳ�
%���׶�ͨ�����


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
    if(e(i) <= e_s2)
        sw(i) = sw(i) + 0*2;
    else
        sw(i) = sw(i) + 1*2;
    end
    %��¼�׶θı�ʱ��ѭ����ֵ
    if(sw(i) ~= sw(i - 1))
        swc(1,j) = i;
        swc(2,j) = sw(i);
        j = j + 1;
    end
    %�ֽ׶μ���ȼ���ܳ���ͨ�����
    switch sw(i)
        case{1}     %�������׶�
            si1 = (l_s1 + r_s1 + e(i))*(beta_s - as1 - as3) + (r_s1 + e(i))*(pi + as1 + as2) ...
                + h_s1 + h_s2 + (r_s2 + e(i))*(pi - as2 - as4) ...
                + (l_s2 - r_s2 - e(i))*(beta_s - as2 - as4) + (e(i) + dc/2)*beta_s;
            s(i) = n_s*si1*Lp;
            Api1 = (l_s1^2)*((sin(2*as1) + sin(2*as2)) / 4) ...
                - (l_s2^2)*((sin(2*as2) + sin(2*as4)) / 4) ...
                + ((l_s1 + r_s1 + e(i))^2)*((beta_s - as1 - as3) / 2) ...
                - ((l_s2 - r_s2 - e(i))^2)*((beta_s - as2 - as4) / 2) ...
                + ((r_s1 + e(i))^2)*((pi + as1 + as3) / 2) ...
                + ((r_s2 + e(i))^2)*((pi - as2 - as4) / 2) ...
                - h_s1*(e_s1 - e(i)) - h_s2*(e_s2 - e(i)) + ((e(i) + dc/2)^2)*beta_s/2;
            Ap(i) = n_s*Api1;
        case{3}     %С�����׶�
            si2 = (l_s1 + r_s1 + e(i))*(beta_s - as1 - as3) ...
                + (r_s1 + e(i))*(pi/2 + as1) + l_s1*cos(as1) - (e(i) + dc/2) ...
                + (e(i) + dc/2)*asin((l_s2*sin(as2) - r_s2 - e(i)) / (e(i) + dc/2)) ...
                + (r_s1 + e(i))*(as3 + asin( l_s1*sin(as3) / (r_s1 + e(i)) ));
            s(i) = n_s*si2*Lp;
            Api2 = ((l_s1 + r_s1 + e(i))^2)*((beta_s - as1 - as3) / 2) ...
                + (((r_s1 + e(i))^2) / 2)*(as3 + asin( l_s1*sin(as3) / (r_s1 + e(i)) )) ...
                + (l_s1*sin(as3) / 2)*( sqrt( (r_s1 + e(i))^2 - (l_s1*sin(as3))^2 ) + l_s1*cos(as3) ) ...
                + (((r_s1 + e(i))^2) / 2)*(pi/2 + as1) + ((l_s1^2) / 2)*sin(as1)*cos(as1) ...
                - (l_s1*sin(as1) - r_s1 - e(i))*(l_s1*cos(as1) - e(i) - dc/2);
            Ap(i) = n_s*Api2;
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
