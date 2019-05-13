%181016
%Dashaway
%ȼ�ղ�������


%190311
%����˫Y��ͼ��

%190329
%˫Y�ᰴ������ʾ


clear;
close all;
%��ʼֵ�趨
FileName = '����ϵ��.xlsx';
SheetName = 'sheet1';
read = xlsread(FileName,SheetName,'C4:D163');
v_f = read(:,1)';
c_f = read(:,2)';

%����

%���㸨������
long = 150000;        %���鳤��
n = 1:1:long;       %��ͼ������
dt = 5e-5;      %����ֵ
%�������
m_c = 120;      %�̶�����(kg)
s_f = 0.1018;        %���������(m^2)
V0 = 0;     %��ʼ�ٶ�(m/s)
%ȼ���Ҳ���
Dr = 0.300;       %ȼ�����⾶(m)
At = 5.026548e-3;     %��ܺ����(m^2)
%װҩ����
D = 0.300;         %ҩ����ֱ��(m)
r = 0.010;        %΢Բ���뾶(m)
l = 0.130;        %���ľ�(m)
R = 0.01;       %���Ŀ�(m)
m_s = 8;        %����ҩ��
n_s = 3;        %������
Lp = 0.500;      %ҩ������(m)
%�ƽ�������
n_p = 0.302;      %ѹǿָ��
rho_p = 1730;       %�ܶ�(kg/m^3)
c = 1600;       %�����ٶ�(m/s)
rb_0 = 5e-3;      %?��ʼȼ��(m/s)
alpha_r = 1.71e-4;        %?ȼ��ϵ��
%��֪����
p0 = 1.02e5;    %��ʼѹǿ(Pa)
gamma = 1.2;        %���ȱ�
phi_alpha = 1;       %?��ʴ����
phi_m = 1;      %?
g = 9.8;        %�������ٶ�(m/s^2)
Ma = 330;       %����(m/s)
rho_air = 1.1691;       %�����ܶ�(kg/m^3)



%����ó���
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���
trd = (2*pi) / (n_s*(m_s + 1));     %ҩԲ�Ľ�
trr = (2*pi*m_s) / (n_s*(m_s + 1));     %��Բ�Ľ�

s_a0 = 2*n_s*trr*l + 2*pi*(R + n_s*r);      %ȼ�����ʼ�ܳ�
Ap_a0 = 2*n_s*trr*l*r + n_s*pi*r^2 + pi*R^2;     %��ʼͨ�����

Vc = (pi*(Dr^2) / 4)*Lp;        %ǻ�����(m^3)
Vp0 = (pi*(Dr^2) / 4 - Ap_a0)*Lp;       %ҩ����ʼ���(m^3)
mp = Vp0*rho_p;     %ҩ������(kg)
v_f = v_f*Ma;
c_f = c_f*1e6;
C_x = c_f(1);
%Լ������
%����Լ������

%�ֽ׶��б�����
e_a = D/2 - r - l;              %�������
e_b = l*sin(trd/2) - r;       %������Բ�ཻ
e_c = sqrt( (D/2 - l*cos(trd/2))^2 + (l*sin(trd/2))^2 ) - r;    %����Բ�����
e_d = (l - r - R)/2;            %������Բ
e_e = ( l^2 + R^2 - r^2 - 2*l*R*cos(trd/2) ) / ( 2*(l*cos(trd/2) + r - R) );    %����Բ����Բ
%���׶ζ�Ӧ����
%     ��һ�׶�
%     e1 = ((e <= e_a) & (e <= e_b));
%     sw = 1;
%     �ڶ�A�׶�
%     e2a = ((e > e_a) & (e < e_b));
%     sw = 3;
%     �ڶ�B�׶�
%     e2b = ((e <= e_a) & (e > e_b));
%     sw = 5;
%     �����׶�
%     e3 = ((e > e_b) & (e <= e_c) & (e <= e_d));
%     sw = 7;
%     ���Ľ׶�
%     e4 = ((e > e_b) & (e > e_c) & (e <= e_d));
%     sw = 15;
%     ����׶�
%     e5 = ((e > e_c) & (e > e_d) & (e <= e_e));
%     sw = 31;
%�����������
ep = max(max( max(e_a,e_b),max(e_c,e_d) ),e_e);


%����
%�����ñ������飨�Ѹ���ֵ��
p = p0*ones(1,long);        %ʵ��ѹǿ(Pa)
rb = rb_0*ones(1,long);     %ȼ��(m/s)
e = 0*ones(1,long);     %����ȥ���(m)
s = s_a0*ones(1,long);     %ȼ����ʵ�ʱ߳�(m)
m_b = rho_p*s_a0*Lp*rb_0*ones(1,long);     %ȼ��������(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %��������(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %����(N)
Ab = s_a0*Lp*ones(1,long);     %ȼ�����(m^2)
Ap = Ap_a0*ones(1,long);     %ͨ�����(m^2)
Vg = Ap*Lp;     %�������(m^3)
m_z = (m_c + mp)*ones(1,long);      %������(kg)
a = ( 2000*(p0*At / c) / (m_c + mp))*ones(1,long);     %���ٶ�(N/m^2)
V = V0*ones(1,long);      %�ٶ�(m/s)
f = zeros(1,long);      %����(N)
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
m_z_max = (m_c + mp);
a_max = ( 2000*(p0*At / c) / (m_c + mp));
V_max = V0;
f_max = 0;
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
m_z_min = (m_c + mp);
a_min = ( 2000*(p0*At / c) / (m_c + mp));
V_min = V0;
f_min = 0;
%ѭ������
i = 1;
j = 1;
k = 1;
sw = zeros(1,long);        %�׶�ѡ��
swc = zeros(2,100);        %�׶α仯��

%�����м�ֵ
%     as1 =  asin( (l*sin(trd/2)) / (r + e) );
%     ac1 = acos( (4*l^2 + D^2 - 4*(r + e)^2) / (4*l*D) ); 
%     ac2 = acos( (4*(r + e)^2 + 4*l^2 - D^2) / (8*l*(r + e)) ); 
%     ac3 = acos( (l^2 + (r + e)^2 - (R + e)^2) / (2*l*(r + e)) );
%     ac4 = acos( (l^2 + (R + e)^2 - (r + e)^2) / (2*l*(R + e)) );
%���׶�ȼ���ܳ�
%     s1 =  2*n_s*trr*l + 2*pi*(n_s*(r + e) + (R + e));
%     s2a = 2*n_s*ac2*(r + e) + 2*n_s*(trr/2)*(l - r - e) + 2*pi*(R + e);
%     s2b =  2*n_s*trr*l + 2*pi*(R + e) + 4*n_s*as1*(r + e);
%     s3 = 2*n_s*(r + e)*(ac2 + 2*as1 - pi) + 2*pi*(R + e) + n_s*trr*(l - r - e);
%     s4 = 2*n_s*(r + e)*(as1 - trd/2) + 2*pi*(R + e) + n_s*trr*(l - r - e);
%     s5 = 2*n_s*(r + e)*(as1 - trd/2 - ac3) + 2*n_s*(R + e)*(trd/2 - ac4);  
%���׶�ͨ�����
%     Ap1 = 2*n_s*trr*l*(r + e) + n_s*pi*(r + e)^2 + pi*(R + e)^2;
%     Ap2a = pi*(R + e)^2 + n_s*ac2*(r + e)^2 ...  
%         + n_s*(((D/2)^2)*ac1 - n_s*l*(D/2)*sin(ac1) ...
%         + (n_s/2)*trr*( (D/2)^2 - (l - r - e)^2 );
%     Ap2b = 2*n_s*trr*l*(r + e) + pi*(R + e)^2 ...
%         + n_s*((r + e)^2)*(2*as1 + sin(2*as1));
%     Ap3 = pi*(R + e)^2 + n_s*ac2*(r + e)^2 ...  
%         + n_s*((D/2)^2)*ac1 - n_s*l*(D/2)*sin(ac1) ...
%         + (n_s/2)*trr*( (D/2)^2 - (l - r - e)^2 ) ....
%         - n_s*((r + e)^2)*(pi - 2*as1 - sin(2*as1));
%     Ap4 = pi*((D/2)^2) + pi*(R + e)^2 - (trr/2)*n_s*((l - r - e)^2);
%         - n_s*l*(r + e)*sin(as1 - trd/2) + n_s*((r + e)^2)*(as1 - trd/2);
%     Ap5 = pi*((D/2)^2);
%����������㹫ʽ
%     dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;

%��ֵ����
%ѭ��ǰ����������
%ѹǿ���в���
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;


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
    if(e(i) <= e_c)
        sw(i) = sw(i) + 0*8;
    else
        sw(i) = sw(i) + 1*8;
    end
    if(e(i) <= e_d)
        sw(i) = sw(i) + 0*16;
    else
        sw(i) = sw(i) + 1*16;
    end
    if(e(i) <= e_e)
        sw(i) = sw(i) + 0*32;
    else
        sw(i) = sw(i) + 1*32;
    end
    %��¼�׶θı�ʱ��ѭ����ֵ
    if(sw(i) ~= sw(i - 1))
        swc(1,j) = i;
        swc(2,j) = sw(i);
        j = j + 1;
    end
    %�ֽ׶μ���ȼ���ܳ���ͨ�����
    switch sw(i)
        case{1}     %��һ�׶�
            s(i) =  2*n_s*trr*l + 2*pi*(n_s*(r + e(i)) + (R + e(i)));
            Ap(i) = 2*n_s*trr*l*(r + e(i)) + n_s*pi*(r + e(i))^2 + pi*(R + e(i))^2;
        case{3}     %�ڶ�A�׶�
            ac1 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / (4*l*D) );
            ac2 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) );
            s(i) = 2*n_s*ac2*(r + e(i)) + 2*n_s*(trr/2)*(l - r - e(i)) + 2*pi*(R + e(i));
            Ap(i) = pi*(R + e(i))^2 + n_s*ac2*(r + e(i))^2 ...
                + n_s*((D/2)^2)*ac1 - n_s*l*(D/2)*sin(ac1) ...
                + (n_s/2)*trr*( (D/2)^2 - (l - r - e(i))^2 );
        case{5}     %�ڶ�B�׶�
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            s(i) =  2*n_s*trr*l + 2*pi*(R + e(i)) + 4*n_s*as1*(r + e(i));
            Ap(i) = 2*n_s*trr*l*(r + e(i))  + pi*(R + e(i))^2 ...
                + n_s*((r + e(i))^2)*(2*as1 + sin(2*as1));
        case{7}     %�����׶�
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            ac1 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / (4*l*D) );
            ac2 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) );
            s(i) = 2*n_s*(r + e(i))*(ac2 + 2*as1 - pi) + 2*pi*(R + e(i)) ...
                + n_s*trr*(l - r - e(i));
            Ap(i) = pi*(R + e(i))^2 + n_s*ac2*(r + e(i))^2 ...
                + n_s*((D/2)^2)*ac1 - n_s*l*(D/2)*sin(ac1) ...
                + (n_s/2)*trr*( (D/2)^2 - (l - r - e(i))^2 ) ....
                - n_s*((r + e(i))^2)*(pi - 2*as1 - sin(2*as1));
        case{15}     %���Ľ׶�
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            s(i) = 2*n_s*(r + e(i))*(as1 - trd/2) + 2*pi*(R + e(i)) + n_s*trr*(l - r - e(i));
            Ap(i) = pi*((D/2)^2) + pi*((R + e(i))^2) - (trr/2)*n_s*((l - r - e(i))^2) ...
                - n_s*l*(r + e(i))*sin(as1 - trd/2) + n_s*((r + e(i))^2)*(as1 - trd/2);
        case{31}     %����׶�
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            ac3 = acos( (l^2 + (r + e(i))^2 - (R + e(i))^2) / (2*l*(r + e(i))) );
            ac4 = acos( (l^2 + (R + e(i))^2 - (r + e(i))^2) / (2*l*(R + e(i))) );
            s(i) = 2*n_s*(r + e(i))*(as1 - trd/2 - ac3) + 2*n_s*(R + e(i))*(trd/2 - ac4);
            Ap(i) = pi*((D/2)^2);
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
    m_z(i) = m_c + (Vc - Vg(i))*rho_p;
    a(i) = (F(i) - f(i - 1)) / m_z(i) ;
    V(i) = V(i - 1) + a(i)*dt;
    
    %��������ϵ��
    if(V(i) > v_f(k))
        if(V(i) < v_f(160))
            k = k + 1;
            C_x = c_f(k);
        end
    end
    
    f(i) = (s_f*rho_air / 2)*C_x*((V(i) / Ma)^2);
    
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
    if(m_z(i) > m_z_max)
        m_z_max = m_z(i);
    end
    if(a(i) > a_max)
        a_max = a(i);
    end
    if(V(i) > V_max)
        V_max = V(i);
    end
    if(f(i) > f_max)
        f_max = f(i);
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
    if(m_z(i) < m_z_min)
        m_z_min = m_z(i);
    end
    if(a(i) < a_min)
        a_min = a(i);
    end
    if(V(i) < V_min)
        V_min = V(i);
    end
    if(f(i) < f_min)
        f_min = f(i);
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
f(i:1:long) = [];
m_z(i:1:long) = [];
a(i:1:long) = [];
V(i:1:long) = [];
sw(i:1:long) = [];
swc(:,j:100) = [];
t = dt*n;
t_max = dt*(i - 1);     %ȼ����ʱ��
prx = 1.05;
pri = -0.05;

str = [' t_max = ',num2str(t_max)];
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
% text(swc(2)*dt,p(swc(2)),['(',num2str(swc(2)*dt),',',num2str(p(swc(2))),')'],'color','g');
% text(swc(3)*dt,p(swc(3)),['(',num2str(swc(3)*dt),',',num2str(p(swc(3))),')'],'color','b');
% text(swc(4)*dt,p(swc(4)),['(',num2str(swc(4)*dt),',',num2str(p(swc(4))),')'],'color','r');


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


figure;
hold on;
plot(t,m_z);
axis ([pri*t_max,prx*t_max,(pri*(m_z_max - m_z_min) + m_z_min), ...
    (prx*(m_z_max - m_z_min) + m_z_min)]);
title('������');
xlabel('ʱ��(s)');
ylabel('����(kg)');
% legend('����2');

figure;
hold on;
plot(t,a);
axis ([pri*t_max,prx*t_max,(pri*(a_max - a_min) + a_min), ...
    (prx*(a_max - a_min) + a_min)]);
title('���ٶ�');
xlabel('ʱ��(s)');
ylabel('���ٶ�(m/s^2)');
% legend('����2');

figure;
hold on;
plot(t,V);
axis ([pri*t_max,prx*t_max,(pri*(V_max - V_min) + V_min), ...
    (prx*(V_max - V_min) + V_min)]);
title('�ٶ�');
xlabel('ʱ��(s)');
ylabel('�ٶ�(m/s)');
% legend('����2');

figure;
hold on;
plot(t,f);
axis ([pri*t_max,prx*t_max,(pri*(f_max - f_min) + f_min), ...
    (prx*(f_max - f_min) + f_min)]);
title('����');
xlabel('ʱ��(s)');
ylabel('��(N)');
% legend('����2');


%����
