%181016
%Dashaway
%ȼ�ղ�������


%190512
%���������������ٶ�


%190825
%���Ļ�ͼ�ã�����װҩ��ʽ�������


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
long = 250000;        %���鳤��
n = 1:1:long;       %��ͼ������
dt = 5e-5;      %����ֵ
%�������
m_c = 120;      %�̶�����(kg)
s_f = 0.1018;        %���������(m^2)
V0 = 0;     %��ʼ�ٶ�(m/s)
%ȼ���Ҳ���
Dr = 0.302;       %ȼ�����ھ�(m)
At = 5.026548e-3;     %��ܺ����(m^2)
%���װҩ����
D = 0.300;         %ҩ����ֱ��(m)
r = 0.010;        %΢Բ���뾶(m)
l = 0.130;        %���ľ�(m)
R = 0.01;       %���Ŀ�(m)
m_s = 8;        %����ҩ��
n_s_mh = 3;        %������
Lp = 1.000;      %ҩ������(m)
%˫�׶�װҩ����
%�ǿ�װҩ����
Ds = 0.300;       %�⾶(m)
ep_s = 0.042;     %���(m)
r_s = 0.005;        %�Ǽ�Բ���뾶(m)
r1_s = 0.001;        %�Ǹ�����Բ���뾶(m)
l_s = 0.103;        %��������(m)
n_s = 12;        %�ǽ���
epsilon_s = 0.7;      %�Ƿ���
theta_s = 0.43;      %�Ǹ����(rad)
%Բ��װҩ����
Dc = 0.300;     %�⾶(m)
dc = 0.115;      %�ھ�(m)
ep_c = (Dc - dc) / 2;       %���(m)
%�ֶβ���
% Lp = 1.000;      %ҩ������(m)
pers = 0.5;     %�ǿ׶�ռ��
Lp_s = Lp*pers;     %�ǿ׶γ���
Lp_c = Lp - Lp_s;       %Բ���γ���
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


%%%%%%%%%%%%%%%%%%%%%%
%���װҩ���ݼ���


%����ó���
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���
trd = (2*pi) / (n_s_mh*(m_s + 1));     %ҩԲ�Ľ�
trr = (2*pi*m_s) / (n_s_mh*(m_s + 1));     %��Բ�Ľ�

s_a0 = 2*n_s_mh*trr*l + 2*pi*(R + n_s_mh*r);      %ȼ�����ʼ�ܳ�
Ap_a0 = 2*n_s_mh*trr*l*r + n_s_mh*pi*r^2 + pi*R^2;     %��ʼͨ�����

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
%     s1 =  2*n_s_mh*trr*l + 2*pi*(n_s_mh*(r + e) + (R + e));
%     s2a = 2*n_s_mh*ac2*(r + e) + 2*n_s_mh*(trr/2)*(l - r - e) + 2*pi*(R + e);
%     s2b =  2*n_s_mh*trr*l + 2*pi*(R + e) + 4*n_s_mh*as1*(r + e);
%     s3 = 2*n_s_mh*(r + e)*(ac2 + 2*as1 - pi) + 2*pi*(R + e) + n_s_mh*trr*(l - r - e);
%     s4 = 2*n_s_mh*(r + e)*(as1 - trd/2) + 2*pi*(R + e) + n_s_mh*trr*(l - r - e);
%     s5 = 2*n_s_mh*(r + e)*(as1 - trd/2 - ac3) + 2*n_s_mh*(R + e)*(trd/2 - ac4);  
%���׶�ͨ�����
%     Ap1 = 2*n_s_mh*trr*l*(r + e) + n_s_mh*pi*(r + e)^2 + pi*(R + e)^2;
%     Ap2a = pi*(R + e)^2 + n_s_mh*ac2*(r + e)^2 ...  
%         + n_s_mh*(((D/2)^2)*ac1 - n_s_mh*l*(D/2)*sin(ac1) ...
%         + (n_s_mh/2)*trr*( (D/2)^2 - (l - r - e)^2 );
%     Ap2b = 2*n_s_mh*trr*l*(r + e) + pi*(R + e)^2 ...
%         + n_s_mh*((r + e)^2)*(2*as1 + sin(2*as1));
%     Ap3 = pi*(R + e)^2 + n_s_mh*ac2*(r + e)^2 ...  
%         + n_s_mh*((D/2)^2)*ac1 - n_s_mh*l*(D/2)*sin(ac1) ...
%         + (n_s_mh/2)*trr*( (D/2)^2 - (l - r - e)^2 ) ....
%         - n_s_mh*((r + e)^2)*(pi - 2*as1 - sin(2*as1));
%     Ap4 = pi*((D/2)^2) + pi*(R + e)^2 - (trr/2)*n_s_mh*((l - r - e)^2);
%         - n_s_mh*l*(r + e)*sin(as1 - trd/2) + n_s_mh*((r + e)^2)*(as1 - trd/2);
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
            s(i) =  2*n_s_mh*trr*l + 2*pi*(n_s_mh*(r + e(i)) + (R + e(i)));
            Ap(i) = 2*n_s_mh*trr*l*(r + e(i)) + n_s_mh*pi*(r + e(i))^2 + pi*(R + e(i))^2;
        case{3}     %�ڶ�A�׶�
            ac1 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / (4*l*D) );
            ac2 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) );
            s(i) = 2*n_s_mh*ac2*(r + e(i)) + 2*n_s_mh*(trr/2)*(l - r - e(i)) + 2*pi*(R + e(i));
            Ap(i) = pi*(R + e(i))^2 + n_s_mh*ac2*(r + e(i))^2 ...
                + n_s_mh*((D/2)^2)*ac1 - n_s_mh*l*(D/2)*sin(ac1) ...
                + (n_s_mh/2)*trr*( (D/2)^2 - (l - r - e(i))^2 );
        case{5}     %�ڶ�B�׶�
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            s(i) =  2*n_s_mh*trr*l + 2*pi*(R + e(i)) + 4*n_s_mh*as1*(r + e(i));
            Ap(i) = 2*n_s_mh*trr*l*(r + e(i))  + pi*(R + e(i))^2 ...
                + n_s_mh*((r + e(i))^2)*(2*as1 + sin(2*as1));
        case{7}     %�����׶�
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            ac1 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / (4*l*D) );
            ac2 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) );
            s(i) = 2*n_s_mh*(r + e(i))*(ac2 + 2*as1 - pi) + 2*pi*(R + e(i)) ...
                + n_s_mh*trr*(l - r - e(i));
            Ap(i) = pi*(R + e(i))^2 + n_s_mh*ac2*(r + e(i))^2 ...
                + n_s_mh*((D/2)^2)*ac1 - n_s_mh*l*(D/2)*sin(ac1) ...
                + (n_s_mh/2)*trr*( (D/2)^2 - (l - r - e(i))^2 ) ....
                - n_s_mh*((r + e(i))^2)*(pi - 2*as1 - sin(2*as1));
        case{15}     %���Ľ׶�
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            s(i) = 2*n_s_mh*(r + e(i))*(as1 - trd/2) + 2*pi*(R + e(i)) + n_s_mh*trr*(l - r - e(i));
            Ap(i) = pi*((D/2)^2) + pi*((R + e(i))^2) - (trr/2)*n_s_mh*((l - r - e(i))^2) ...
                - n_s_mh*l*(r + e(i))*sin(as1 - trd/2) + n_s_mh*((r + e(i))^2)*(as1 - trd/2);
        case{31}     %����׶�
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            ac3 = acos( (l^2 + (r + e(i))^2 - (R + e(i))^2) / (2*l*(r + e(i))) );
            ac4 = acos( (l^2 + (R + e(i))^2 - (r + e(i))^2) / (2*l*(R + e(i))) );
            s(i) = 2*n_s_mh*(r + e(i))*(as1 - trd/2 - ac3) + 2*n_s_mh*(R + e(i))*(trd/2 - ac4);
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
    if(k > 1)
        if(V(i) < v_f(k - 1))
            k = k - 1;
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

F(i:1:long) = [];

V(i:1:long) = [];

t = dt*n;
t_max = dt*(i - 1);     %ȼ����ʱ��



str = [' ���ȼ����ʱ�䣺 ',num2str(t_max),'s'];
disp(str);

%�����ݴ�

i_mh = i;
t_mh = t;
t_mh_max = t_max;
F_mh = F;
F_max_mh = F_max;
F_min_mh = F_min;
V_mh = V;
V_max_mh = V_max;
V_min_mh = V_min;



%%%%%%%%%%%%%%%%%%%%%%


%˫�׶�װҩ���ݼ���
%���㸨������

n = 1:1:long;       %��ͼ������


%����ó���
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %���ȱȺ���
beta_s = pi / n_s;        %�ȷֽ�(rad)
% v_f = v_f*Ma;
% c_f = c_f*1e6;
% C_x = c_f(1);

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

Vc = (pi*(Dr^2) / 4)*Lp;        %ǻ�����(m^3)
Vps0 = Vc - Vp0;       %ҩ����ʼ���(m^3)

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
m_z = (m_c + mp)*ones(1,long);      %������(kg)
a = ( 2000*(p0*At / c) / (m_c + mp))*ones(1,long);     %���ٶ�(N/m^2)
V = V0*ones(1,long);      %�ٶ�(m/s)
f = zeros(1,long);      %����(N)
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
m_z_max = (m_c + mp);
a_max = ( 2000*(p0*At / c) / (m_c + mp));
V_max = V0;
f_max = 0;
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
            if(Ap_s(i - 1 ) < (pi*Dr^2 / 4) )
                Ap_s(i) = Ap_s(i - 1) +0.0008*Apsi3;
            else
                Ap_s(i) = pi*Dr^2 / 4;
            end
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
    if(k > 1)
        if(V(i) < v_f(k - 1))
            k = k - 1;
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

F(i:1:long) = [];
V(i:1:long) = [];

t = dt*n;
t_max = dt*(i - 1);     %ȼ����ʱ��
prx = 1.05;
pri = -0.05;



str = [' ȼ����ʱ�䣺 ',num2str(t_max),'s'];
disp(str);

%�����ݴ�

i_ds = i;
t_ds = t;
t_ds_max = t_max;
F_ds = F;
F_max_ds = F_max;
F_min_ds = F_min;
V_ds = V;
V_max_ds = V_max;
V_min_ds = V_min;


%%%%%%%%%%%%%%%%%%%%%%

%���ݴ���
t_max = max(t_mh_max,t_ds_max);
i_max = max(i_mh,i_ds);

F_max = max(F_max_mh,F_max_ds);
F_min = min(F_min_mh,F_min_ds);
V_max = max(V_max_mh,V_max_ds);
V_min = min(V_min_mh,V_min_ds);

% n(i_max:1:long) = [];
% t = dt*n;
% F_mh(i_max:1:long) = [];
% V_mh(i_max:1:long) = [];
% F_ds(i_max:1:long) = [];
% V_ds(i_max:1:long) = [];


%�������


figure;
hold on;
plot(t_mh,F_mh,'color','k','LineStyle','-','LineWidth',0.8);
plot(t_ds,F_ds,'color','k','LineStyle','-.','LineWidth',1);
box on;
axis ([pri*t_max,prx*t_max,(pri*(F_max - F_min) + F_min), ...
    (prx*(F_max - F_min) + F_min)]);
% title('F/N');
xlabel('t/s'); 
ylabel('F/N');
legend('case7','case8');

figure;
hold on;
plot(t_mh,V_mh,'color','k','LineStyle','-','LineWidth',0.8);
plot(t_ds,V_ds,'color','k','LineStyle','-.','LineWidth',1);
box on;
axis ([pri*t_max,prx*t_max,(pri*(V_max - V_min) + V_min), ...
    (prx*(V_max - V_min) + V_min)]);
% title('�ٶ�');
xlabel('t/s');
ylabel('V/(m/s)');
legend('case7','case8');


disp('done');
