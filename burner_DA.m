%181016
%Dashaway
%ȼ�ղ�������


%190212
%��Ϊ�ܵ����ڿ�
%���Ĳ���

%190213
%����Բ���Ƕ�


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
D = Dr;         %ҩ���⾶(m)
r = 0.003;        %΢Բ��(m)
l = 0.065;        %���ľ�(m)
R = 0.01;       %���Ŀ�(m)
m_s = 8;        %��ҩ��
n_s = 4;    %ҩ������
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
trd = (2*pi) / (n_s*(m_s + 1));     %ҩԲ�Ľ�
trr = (2*pi*m_s) / (n_s*(m_s + 1));     %��Բ�Ľ�
%�б�����
e_a = D/2 - r - l;              %1��2�׶Ρ������
e_b = l*sin(trd/2) - r;       %2��3�׶Ρ�����Բ�ཻ
e_c = (l - r - R)/2;            %4�׶Ρ�����Բ
e_d = sqrt( (D/2 - l*cos(trd/2))^2 + (l*sin(trd/2))^2 ) - r;    %3��4�׶Ρ���Բ����
%     e1 = ((e < e_b) & (e < e_a));
%     e2a = ((e < e_b) & (e > e_a));
%     e2b = ((e > e_b) & (e < e_a));
%     e3 = ((e > e_b) & (e < e_c) & (e < e_d));
%     e4 = ((e > e_b) & (e < e_c) & (e > e_d));

%�����м�ֵ
%     as1 =  asin((D/2 - l) / (r + e));
%     as2 =  asin( (l*sin(trd/2)) / (r + e) );
%     ac1 = acos( (4*l^2 + D^2 - 4*(r + e)^2) / (4*l*D) ); 
%     ac2 = acos( (4*(r + e)^2 + 4*l^2 - D^2) / (8*l*(r + e)) ); 
% 
%     s_a =  (4*pi*l*m_s) / (m_s + 1) + 2*pi*((R + e) + n_s*(r + e));    %��һ�׶�ȼ���ܳ���ʽ
%     s_b = 2*n_s*(as1 + pi/2)*(r + e) + 2*n_s*(trr/2)*(l + R - r);  %�ڶ��׶�ȼ���ܳ���ʽ
%     s_c = 2*n_s*(r + e)*(ac2 + 2*as2 - pi) + 2*pi*R + 2*n_s*trr*(l - r);  %�����׶�ȼ���ܳ���ʽ
%     s_d = 2*n_s*(r + e)*(as2 - trd/2) + 2*pi*(R + e) + 2*n_s*trr*(l - r - e);  %���Ľ׶�ȼ���ܳ���ʽ
% 
%     Ap_a = (4*pi*l*m_s*(r + e)) / (m_s + 1) ... %��һ�׶�ͨ�����
%         + n_s*pi*(r + e)^2 + pi*(R + e)^2;   
%     Ap_b = n_s*(trr/2)*(R + e)^2 + n_s*(as1 + pi/2)*(r + e)^2 ...  %�ڶ��׶�ͨ�����
%         + n_s*D*ac1 - n_s*l*(D/2)*sin(ac1) + 2*n_s*l*(r + e)*trr;
%     %Ap_c = pi*(Dr^2)/4 - pi*(R + e)^2 - n_s*Dr*(pi/2 - ac1) ...%�����׶�ͨ�����
%         %*(  Dr/2 - (r + e)*sin(as + pi / (2*n_s)) / sin(pi / (2*n_s))  );
%     Ap_c = pi*(D^2)/4 - pi*(R + e)^2;  %�����׶�ͨ�����
%     Ap_d = pi*(D^2)/4 - pi*(R + e)^2;  %���Ľ׶�ͨ�����
    
%�ܳ�����
s_a0 = (4*pi*l*m_s) / (m_s + 1) + 2*pi*(R + n_s*r);      %��һ�׶�ȼ�����ʼ�߳�

%ͨ���������
Ap_a0 = (4*pi*l*m_s*r) / (m_s + 1) + n_s*pi*r^2 + pi*R^2;     %��һ�׶γ�ʼͨ�����

    
        

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
j = 1;
sw = 0*ones(1,long);        %�׶�ѡ��
swc = 0*ones(1,100);
bit1 = 2;
bit2 = 4;
bit3 = 8;
bit4 = 16;

%��ֵ����
%ѭ��ǰ����������
%ѹǿ���в���
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;
ep = max(max(e_a,e_b),max(e_c,e_d));

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
    rb(i) = alpha_r*p(i)^n_p;
    e(i) = e(i - 1) + rb(i)*dt;

    %�жϴ�ʱ�����ĸ��׶�
    sw(i) = 1;
    if(e(i) < e_a)
        sw(i) = sw(i) + 0*bit1;
    else
        sw(i) = sw(i) + 1*bit1;
    end
    if(e(i) < e_b)
        sw(i) = sw(i) + 0*bit2;
    else
        sw(i) = sw(i) + 1*bit2;
    end
    if(e(i) < e_c)
        sw(i) = sw(i) + 0*bit3;
    else
        sw(i) = sw(i) + 1*bit3;
    end
    if(e(i) < e_d)
        sw(i) = sw(i) + 0*bit4;
    else
        sw(i) = sw(i) + 1*bit4;
    end
    %��¼�׶θı�ʱ����ֵλ��
    if(sw(i) ~= sw(i - 1))
        swc(j) = i;
        j = j + 1;
    end
    
    switch sw(i)
        case{1}     %��һ�׶�
            s(i) =  2*n_s*trr*l + 2*pi*(n_s*(r + e(i)) + (R + e(i)));       
            Ap(i) = 2*n_s*trr*l*(r + e(i)) + n_s*pi*(r + e(i))^2 + pi*(R + e(i))^2; 
        case{3}     %�ڶ�A�׶�
            ac1 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / (4*l*D) );
            ac2 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) ); 
            s(i) = 2*n_s*(ac2)*(r + e(i)) + 2*n_s*(trr/2)*(l - r - e(i)) + 2*pi*(R + e(i));
            Ap(i) = pi*(R + e(i))^2 + n_s*(ac2)*(r + e(i))^2 ...  %�ڶ��׶�ͨ�����
                + n_s*((D^2) / 4)*ac1 - n_s*l*(D/2)*sin(ac1) ...
                + (n_s/2)*trr*( (D/2)^2 - (l - r - e(i))^2 );
        case{5}     %�ڶ�B�׶�
            as2 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            s(i) =  2*n_s*trr*l + 2*pi*(n_s*(r + e(i)) + (R + e(i))) ...
                - (pi - 2*as2)*(r + e(i));
            Ap(i) = 2*n_s*trr*l*(r + e(i)) + n_s*pi*(r + e(i))^2 + pi*(R + e(i))^2 ...
                - n_s*((r  + e(i))^2)*(pi - 2*as2 - sin(2*as2)); 
        case{7}     %�����׶�
            as2 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            ac1 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / (4*l*D) ); 
            ac2 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) ); 
            s(i) = 2*n_s*(r + e(i))*(ac2 + 2*as2 - pi) + 2*pi*(R + e(i)) ... %�����׶�ȼ���ܳ���ʽ
                + n_s*trr*(l - r - e(i));
            Ap(i) = pi*(R + e(i))^2 + n_s*(ac2)*(r + e(i))^2 ...  %�����׶�ͨ�����
                + n_s*((D^2) / 4)*ac1 - n_s*l*(D/2)*sin(ac1) ...
                + (n_s/2)*trr*( (D/2)^2 - (l - r - e(i))^2 ) ....
                - n_s*((r  + e(i))^2)*(pi - 2*as2 - sin(2*as2));
            
        case{23}     %���Ľ׶�
            as2 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            s(i) = 2*n_s*(r + e(i))*(as2 - trd/2) + 2*pi*(R + e(i)) + n_s*trr*(l - r - e(i));  %���Ľ׶�ȼ���ܳ���ʽ
            Ap(i) = pi*((D/2)^2) + pi*((R + e(i))^2) - (trr/2)*n_s*((l - r - e(i))^2) ...%���Ľ׶�ͨ�����
                - n_s*(l - r - e(i))*(r + e(i))*sin(as2 - trd/2);
        otherwise
            s(i) = s(i - 1);
            Ap(i) = Ap(i - 1);
    end
    
    %�ֽ׶μ���ȼ���ܳ���ͨ�����
%     if(e < e_a)
%         s(i) = 2*pi*(R + l + e(i)) + 2*pi*n_s*(r + e(i));
%         Ap(i) = 2*pi*l*(r + e(i)) + n_s*pi*(r + e(i))^2 + pi*(R + e(i))^2;
%     elseif(e < e_b)
%         as =  asin( (l*sin(pi / (2*n_s))) / (r + e(i)) );
%         s(i) = 4*n_s*(r + e(i))*as + 2*pi*(l + R + e(i));
%         Ap(i) = 2*n_s*((r + e(i))^2)*as + 2*pi*l*(r + e(i)) + pi*(R + e(i))^2 ...
%             + n_s*((r + e(i))^2)*sin(2*as) ;
%     else
%         as =  asin( (l*sin(pi / (2*n_s))) / (r + e(i)) );
%         ac1 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) ); 
%         ac2 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / 4*l*D ); 
%         s(i) = 2*n_s*(r + e(i))*(ac1 + 2*as - pi) + 2*pi*(l + R - r);
%         Ap(i) = pi*(D^2)/4 - pi*(R + e(i))^2 - n_s*D*(pi/2 -ac2) ...
%             *( D/2 - (r + e(i))*sin(as + pi / (2*n_s)) / sin(pi / (2*n_s))  );
%     end

 
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
sw(i:1:long) = [];

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

