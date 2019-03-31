%181016
%Dashaway
%ȼ�ղ�������


%190321
%�ֶ�ȼ��

%190326
%����㲻ͬװҩ


clear;
close all;
%��ʼֵ�趨


%����

%���㸨������
long = 150000;        %���鳤��
n = 1:1:long;       %��ͼ������
dt = 1e-4;      %����ֵ
%ȼ���Ҳ���
Dr = 0.300;       %ȼ�����⾶(m)
At = 5.026548e-3;     %��ܺ����(m^2)
%װҩ����
%����Բ��װҩ����
Dc = 0.300;     %�⾶(m)
dm = 0.100;     %�о�(m)
dc = 0.050;      %�ھ�(m)
ep_m = (dm - dc) / 2;       %�����(m)
ep_c = (Dc - dc) / 2;       %�����(m)
Lp = 0.500;      %ҩ������(m)
%���ƽ�������
n_p1 = 0.4;      %ѹǿָ��
rho_p1 = 1700;       %�ܶ�(kg/m^3)
c1 = 1584;       %�����ٶ�(m/s)
rb1_0 = 8.77e-5;      %?��ʼȼ��(m/s)
alpha_r1 = 3e-5;        %?ȼ��ϵ��
%���ƽ�������
n_p2 = 0.302;      %ѹǿָ��
rho_p2 = 1730;       %�ܶ�(kg/m^3)
c2 = 1600;       %�����ٶ�(m/s)
rb2_0 = 5e-3;      %?��ʼȼ��(m/s)
alpha_r2 = 1.71e-4;        %?ȼ��ϵ��
%��֪����
p0 = 1.02e5;    %��ʼѹǿ(Pa)
gamma = 1.2;        %���ȱ�
phi_alpha = 1;       %?��ʴ����
phi_m = 1;      %?

%����ó���
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���

%�����м�ֵ

%ѹǿ���в���
p_a1 = rho_p1*alpha_r1*phi_alpha*Gamma^2*c1^2;
p_b1 = phi_m*Gamma^2*c1*At;
p_a2 = rho_p2*alpha_r2*phi_alpha*Gamma^2*c2^2;
p_b2 = phi_m*Gamma^2*c2*At;

%��ʼ����
%Բ��װҩ
%�ܳ�����
s_0 = pi*dc;
%ͨ���������
Ap_0 = pi*(dc^2) / 4;


Vp10 = (pi*(dm^2) / 4 - pi*(dc^2) / 4)*Lp;       %��ҩ�����(m^3)
Vp20 = (pi*(Dr^2) / 4 - pi*(dm^2) / 4)*Lp;       %��ҩ�����(m^3)
Vp0 = Vp10 + Vp20;      %��ҩ�����(m^3)
mp1 = Vp10*rho_p1;     %��ҩ������(kg)
mp2 = Vp20*rho_p2;     %��ҩ������(kg)
mp = mp1 + mp2;     %��ҩ������(kg)

%ѹǿ��
p_a = p_a1;
p_b = p_b1;
n_p = n_p1;
alpha_r = alpha_r1;

%Լ������
%����Լ������

%�ֽ׶��б�����


%���׶ζ�Ӧ����
%     �ǿ׵�һ�׶Σ�Բ����һ�׶�
%     e1 = ((e <= ep_m) & (e <= ep_c));
%     sw = 1;
%     �ǿ׵ڶ��׶Σ�Բ����һ�׶�
%     e2 = ((e > ep_m) & (e <= ep_c));
%     sw = 3;
%�����������
ep = ep_c;
%     e >= ep;


%����
%�����ñ������飨�Ѹ���ֵ��
p = p0*ones(1,long);        %ʵ��ѹǿ(Pa)
rb = rb1_0*ones(1,long);     %ȼ��(m/s)
e = 0*ones(1,long);     %����ȥ���(m)
s = s_0*ones(1,long);     %ȼ����ʵ���ܳ�(m)
m_b = rho_p1*s_0*Lp*rb1_0*ones(1,long);     %ȼ��������(kg/s)
m_p = (phi_m*p0*At / c1)*ones(1,long);     %��������(kg/s)
F = 2000*(p0*At / c1)*ones(1,long);     %����(N)
Ab = s_0*Lp*ones(1,long);     %ȼ�����(m^2)
Ap = Ap_0*ones(1,long);     %ͨ�����(m^2)
Vg = Ap*Lp;     %�������(m^3)
%���ֵ����
p_max = p0;        %���ѹǿ(Pa)
e_max = 0;     %�������ȥ���(m)
s_max = s_0;     %ȼ��������ܳ�(m)
m_b_max = rho_p1*s_0*Lp*rb1_0;     %���ȼ��������(kg/s)
m_p_max = (phi_m*p0*At / c1);     %�����������(kg/s)
F_max = 2000*(p0*At / c1);     %�������(N)
Ab_max = s_0*Lp;     %���ȼ�����(m^2)
Ap_max = Ap_0;     %���ͨ�����(m^2)
Vg_max = Ap_0*Lp;     %����������(m^3)
%��Сֵ����
p_min = p0;        %��Сѹǿ(Pa)
e_min = 0;     %��С����ȥ���(m)
s_min = s_0;     %ȼ������С�ܳ�(m)
m_b_min = rho_p1*s_0*Lp*rb1_0;     %��Сȼ��������(kg/s)
m_p_min = (phi_m*p0*At / c1);     %��С��������(kg/s)
F_min = 2000*(p0*At / c1);     %��С����(N)
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
    s(i) = pi*(dc + 2*e(i));
    Ap(i) = pi*(dc + 2*e(i))^2 / 4;
    %�жϴ�ʱ�����ĸ��׶�
    sw(i) = 1;
    if(e(i) <= ep_m)
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
        case{1}     %��װҩ�׶�
            p_a = p_a1;
            p_b = p_b1;
            n_p = n_p1;
            alpha_r = alpha_r1;
            m_b(i) = rho_p1*Ab(i)*rb(i);
            m_p(i) = phi_m*p(i)*At / c1;
        case{3}     %��װҩ�׶�
            p_a = p_a2;
            p_b = p_b2;
            n_p = n_p2;
            alpha_r = alpha_r2;
            m_b(i) = rho_p2*Ab(i)*rb(i);
            m_p(i) = phi_m*p(i)*At / c2;
        otherwise     %�����׶�
            m_b(i) = m_b(i - 1);
            m_p(i) = m_p(i - 1);
    end
    %������������
    Ab(i) = s(i)*Lp;
    Vg(i) = Ap(i)*Lp;
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


%˫y��2014��
figure;
hold on;
%AX(1)�� AX(2)�ֱ������� axes �ľ���������� set()��������
[AX,H1,H2] = plotyy(t,p,t,F); 
set(AX,'Xlim',[pri*t_max,prx*t_max]);
xlabel('ʱ��(s)');
title('ѹǿ������');
set(AX(:),'Ycolor','k');

set(get(AX(1),'Ylabel'),'string','ѹǿ(Pa)','color','k','linewidth',1.0); 
set(get(AX(2),'Ylabel'),'string','��(N)','color','k','linewidth',1.0); 

per = 0.6;      %������
cou = p_max;        %ѡ����
mag = floor(log10(cou / per));        %��������
fir = floor((cou / per) / (10^mag));      %ȡ��һλ
sec = floor((cou / per - fir*10^mag) / (10^(mag - 1)));       %ȡ�ڶ�λ
ran = fir*10^mag + sec*10^(mag - 1);        %������
cal = fir*10^(mag - 1);     %���̶�
sta = -0.5*fir*10^(mag - 1);        %�����
set(AX(1),'Ylim',[sta,ran],'yTick',0:cal:ran);

per = 0.8;      %������
cou = F_max;        %ѡ����
mag = floor(log10(cou / per));        %��������
fir = floor((cou / per) / (10^mag));      %ȡ��һλ
sec = floor((cou / per - fir*10^mag) / (10^(mag - 1)));       %ȡ�ڶ�λ
ran = fir*10^mag + sec*10^(mag - 1);        %������
cal = fir*10^(mag - 1);     %���̶�
sta = -0.5*fir*10^(mag - 1);        %�����
set(AX(2),'Ylim',[sta,ran],'yTick',0:cal:ran);

box off;
set(H1,'LineStyle','-','color','k','linewidth',1.0);
set(H2,'LineStyle','-.','color','k','linewidth',1.0);
legend('ѹǿ','����');

%����
