%181016
%Dashaway
%燃烧参数计算


%190512
%加入阻力，计算速度


%190825
%论文绘图用，两种装药方式混叠计算


clear;
close all;
%初始值设定
FileName = '阻力系数.xlsx';
SheetName = 'sheet1';
read = xlsread(FileName,SheetName,'C4:D163');
v_f = read(:,1)';
c_f = read(:,2)';
%常量

%计算辅助常量
long = 250000;        %数组长度
n = 1:1:long;       %绘图横坐标
dt = 5e-5;      %步进值
%弹体参数
m_c = 120;      %固定质量(kg)
s_f = 0.1018;        %阻力截面积(m^2)
V0 = 0;     %初始速度(m/s)
%燃烧室参数
Dr = 0.302;       %燃烧室内径(m)
At = 5.026548e-3;     %喷管喉部面积(m^2)
%多孔装药参数
D = 0.300;         %药柱外直径(m)
r = 0.010;        %微圆弧半径(m)
l = 0.130;        %弧心距(m)
R = 0.01;       %中心孔(m)
m_s = 8;        %弧孔药比
n_s_mh = 3;        %弧数量
Lp = 1.000;      %药柱长度(m)
%双阶段装药参数
%星孔装药参数
Ds = 0.300;       %外径(m)
ep_s = 0.042;     %肉厚(m)
r_s = 0.005;        %星尖圆弧半径(m)
r1_s = 0.001;        %星根过渡圆弧半径(m)
l_s = 0.103;        %特征长度(m)
n_s = 12;        %星角数
epsilon_s = 0.7;      %角分数
theta_s = 0.43;      %星根半角(rad)
%圆柱装药参数
Dc = 0.300;     %外径(m)
dc = 0.115;      %内径(m)
ep_c = (Dc - dc) / 2;       %肉厚(m)
%分段参数
% Lp = 1.000;      %药柱长度(m)
pers = 0.5;     %星孔段占比
Lp_s = Lp*pers;     %星孔段长度
Lp_c = Lp - Lp_s;       %圆柱段长度
%推进剂参数
n_p = 0.302;      %压强指数
rho_p = 1730;       %密度(kg/m^3)
c = 1600;       %特征速度(m/s)
rb_0 = 5e-3;      %?初始燃速(m/s)
alpha_r = 1.71e-4;        %?燃速系数
%已知常量
p0 = 1.02e5;    %初始压强(Pa)
gamma = 1.2;        %比热比
phi_alpha = 1;       %?侵蚀函数
phi_m = 1;      %?
g = 9.8;        %重力加速度(m/s^2)
Ma = 330;       %音速(m/s)
rho_air = 1.1691;       %空气密度(kg/m^3)


%%%%%%%%%%%%%%%%%%%%%%
%多孔装药数据计算


%计算得常量
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数
trd = (2*pi) / (n_s_mh*(m_s + 1));     %药圆心角
trr = (2*pi*m_s) / (n_s_mh*(m_s + 1));     %弧圆心角

s_a0 = 2*n_s_mh*trr*l + 2*pi*(R + n_s_mh*r);      %燃烧面初始周长
Ap_a0 = 2*n_s_mh*trr*l*r + n_s_mh*pi*r^2 + pi*R^2;     %初始通气面积

Vc = (pi*(Dr^2) / 4)*Lp;        %腔室体积(m^3)
Vp0 = (pi*(Dr^2) / 4 - Ap_a0)*Lp;       %药柱初始体积(m^3)
mp = Vp0*rho_p;     %药柱质量(kg)
v_f = v_f*Ma;
c_f = c_f*1e6;
C_x = c_f(1);
%约束条件
%参数约束条件

%分阶段判别条件
e_a = D/2 - r - l;              %弧贴外壁
e_b = l*sin(trd/2) - r;       %两弧半圆相交
e_c = sqrt( (D/2 - l*cos(trd/2))^2 + (l*sin(trd/2))^2 ) - r;    %弧半圆贴外壁
e_d = (l - r - R)/2;            %弧贴内圆
e_e = ( l^2 + R^2 - r^2 - 2*l*R*cos(trd/2) ) / ( 2*(l*cos(trd/2) + r - R) );    %弧半圆贴内圆
%各阶段对应条件
%     第一阶段
%     e1 = ((e <= e_a) & (e <= e_b));
%     sw = 1;
%     第二A阶段
%     e2a = ((e > e_a) & (e < e_b));
%     sw = 3;
%     第二B阶段
%     e2b = ((e <= e_a) & (e > e_b));
%     sw = 5;
%     第三阶段
%     e3 = ((e > e_b) & (e <= e_c) & (e <= e_d));
%     sw = 7;
%     第四阶段
%     e4 = ((e > e_b) & (e > e_c) & (e <= e_d));
%     sw = 15;
%     第五阶段
%     e5 = ((e > e_c) & (e > e_d) & (e <= e_e));
%     sw = 31;
%程序结束条件
ep = max(max( max(e_a,e_b),max(e_c,e_d) ),e_e);


%变量
%计算用变量数组（已赋初值）
p = p0*ones(1,long);        %实际压强(Pa)
rb = rb_0*ones(1,long);     %燃速(m/s)
e = 0*ones(1,long);     %已烧去肉厚(m)
s = s_a0*ones(1,long);     %燃烧面实际边长(m)
m_b = rho_p*s_a0*Lp*rb_0*ones(1,long);     %燃气生成率(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %质量流率(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %推力(N)
Ab = s_a0*Lp*ones(1,long);     %燃烧面积(m^2)
Ap = Ap_a0*ones(1,long);     %通气面积(m^2)
Vg = Ap*Lp;     %自由体积(m^3)
m_z = (m_c + mp)*ones(1,long);      %总质量(kg)
a = ( 2000*(p0*At / c) / (m_c + mp))*ones(1,long);     %加速度(N/m^2)
V = V0*ones(1,long);      %速度(m/s)
f = zeros(1,long);      %阻力(N)
%最大值变量
p_max = p0;        %最大压强(Pa)
e_max = 0;     %最大已烧去肉厚(m)
s_max = s_a0;     %燃烧面最大周长(m)
m_b_max = rho_p*s_a0*Lp*rb_0;     %最大燃气生成率(kg/s)
m_p_max = (phi_m*p0*At / c);     %最大质量流率(kg/s)
F_max = 2000*(p0*At / c);     %最大推力(?)
Ab_max = s_a0*Lp;     %最大燃烧面积(m^2)
Ap_max = Ap_a0;     %最大通气面积(m^2)
Vg_max = Ap_max*Lp;     %最大自由体积(m^3)
m_z_max = (m_c + mp);
a_max = ( 2000*(p0*At / c) / (m_c + mp));
V_max = V0;
f_max = 0;
%最小值变量
p_min = p0;        %最小压强(Pa)
e_min = 0;     %最小已烧去肉厚(m)
s_min = s_a0;     %燃烧面最小周长(m)
m_b_min = rho_p*s_a0*Lp*rb_0;     %最小燃气生成率(kg/s)
m_p_min = (phi_m*p0*At / c);     %最小质量流率(kg/s)
F_min = 2000*(p0*At / c);     %最小推力(?)
Ab_min = s_a0*Lp;     %最小燃烧面积(m^2)
Ap_min = Ap_a0;     %最小通气面积(m^2)
Vg_min = Ap_min*Lp;     %最小自由体积(m^3)
m_z_min = (m_c + mp);
a_min = ( 2000*(p0*At / c) / (m_c + mp));
V_min = V0;
f_min = 0;
%循环变量
i = 1;
j = 1;
k = 1;
sw = zeros(1,long);        %阶段选择
swc = zeros(2,100);        %阶段变化点

%计算中间值
%     as1 =  asin( (l*sin(trd/2)) / (r + e) );
%     ac1 = acos( (4*l^2 + D^2 - 4*(r + e)^2) / (4*l*D) ); 
%     ac2 = acos( (4*(r + e)^2 + 4*l^2 - D^2) / (8*l*(r + e)) ); 
%     ac3 = acos( (l^2 + (r + e)^2 - (R + e)^2) / (2*l*(r + e)) );
%     ac4 = acos( (l^2 + (R + e)^2 - (r + e)^2) / (2*l*(R + e)) );
%各阶段燃烧周长
%     s1 =  2*n_s_mh*trr*l + 2*pi*(n_s_mh*(r + e) + (R + e));
%     s2a = 2*n_s_mh*ac2*(r + e) + 2*n_s_mh*(trr/2)*(l - r - e) + 2*pi*(R + e);
%     s2b =  2*n_s_mh*trr*l + 2*pi*(R + e) + 4*n_s_mh*as1*(r + e);
%     s3 = 2*n_s_mh*(r + e)*(ac2 + 2*as1 - pi) + 2*pi*(R + e) + n_s_mh*trr*(l - r - e);
%     s4 = 2*n_s_mh*(r + e)*(as1 - trd/2) + 2*pi*(R + e) + n_s_mh*trr*(l - r - e);
%     s5 = 2*n_s_mh*(r + e)*(as1 - trd/2 - ac3) + 2*n_s_mh*(R + e)*(trd/2 - ac4);  
%各阶段通气面积
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
%龙格库塔计算公式
%     dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;

%数值计算
%循环前常参数计算
%压强项中参数
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;


%循环计算
while (e(i) <= ep)
    %龙格库塔逐步计算，计算出新的压强值
    %dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;
    k1 = p_a*(Ab(i) / Vg(i))*( p(i)^n_p ) - p_b*p(i) / Vg(i);
    k2 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k1/2)^n_p ) - p_b*(p(i) + dt*k1/2) / Vg(i);
    k3 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k2/2)^n_p ) - p_b*(p(i) + dt*k2/2) / Vg(i);
    k4 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k3)^n_p )  - p_b*(p(i) + dt*k3) / Vg(i);
    dp = (k1 + 2*k2 + 2*k3 + k4)/6;
     i = i + 1;
    p(i) = p(i - 1) + dp*dt;

    %参数修正
    %逐项计算其他参数的新值
    rb(i) = alpha_r*p(i)^n_p;
    e(i) = e(i - 1) + rb(i)*dt;

    %判断此时处于哪个阶段
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
    %记录阶段改变时的循环数值
    if(sw(i) ~= sw(i - 1))
        swc(1,j) = i;
        swc(2,j) = sw(i);
        j = j + 1;
    end
    %分阶段计算燃烧周长和通气面积
    switch sw(i)
        case{1}     %第一阶段
            s(i) =  2*n_s_mh*trr*l + 2*pi*(n_s_mh*(r + e(i)) + (R + e(i)));
            Ap(i) = 2*n_s_mh*trr*l*(r + e(i)) + n_s_mh*pi*(r + e(i))^2 + pi*(R + e(i))^2;
        case{3}     %第二A阶段
            ac1 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / (4*l*D) );
            ac2 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) );
            s(i) = 2*n_s_mh*ac2*(r + e(i)) + 2*n_s_mh*(trr/2)*(l - r - e(i)) + 2*pi*(R + e(i));
            Ap(i) = pi*(R + e(i))^2 + n_s_mh*ac2*(r + e(i))^2 ...
                + n_s_mh*((D/2)^2)*ac1 - n_s_mh*l*(D/2)*sin(ac1) ...
                + (n_s_mh/2)*trr*( (D/2)^2 - (l - r - e(i))^2 );
        case{5}     %第二B阶段
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            s(i) =  2*n_s_mh*trr*l + 2*pi*(R + e(i)) + 4*n_s_mh*as1*(r + e(i));
            Ap(i) = 2*n_s_mh*trr*l*(r + e(i))  + pi*(R + e(i))^2 ...
                + n_s_mh*((r + e(i))^2)*(2*as1 + sin(2*as1));
        case{7}     %第三阶段
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            ac1 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / (4*l*D) );
            ac2 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) );
            s(i) = 2*n_s_mh*(r + e(i))*(ac2 + 2*as1 - pi) + 2*pi*(R + e(i)) ...
                + n_s_mh*trr*(l - r - e(i));
            Ap(i) = pi*(R + e(i))^2 + n_s_mh*ac2*(r + e(i))^2 ...
                + n_s_mh*((D/2)^2)*ac1 - n_s_mh*l*(D/2)*sin(ac1) ...
                + (n_s_mh/2)*trr*( (D/2)^2 - (l - r - e(i))^2 ) ....
                - n_s_mh*((r + e(i))^2)*(pi - 2*as1 - sin(2*as1));
        case{15}     %第四阶段
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            s(i) = 2*n_s_mh*(r + e(i))*(as1 - trd/2) + 2*pi*(R + e(i)) + n_s_mh*trr*(l - r - e(i));
            Ap(i) = pi*((D/2)^2) + pi*((R + e(i))^2) - (trr/2)*n_s_mh*((l - r - e(i))^2) ...
                - n_s_mh*l*(r + e(i))*sin(as1 - trd/2) + n_s_mh*((r + e(i))^2)*(as1 - trd/2);
        case{31}     %第五阶段
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            ac3 = acos( (l^2 + (r + e(i))^2 - (R + e(i))^2) / (2*l*(r + e(i))) );
            ac4 = acos( (l^2 + (R + e(i))^2 - (r + e(i))^2) / (2*l*(R + e(i))) );
            s(i) = 2*n_s_mh*(r + e(i))*(as1 - trd/2 - ac3) + 2*n_s_mh*(R + e(i))*(trd/2 - ac4);
            Ap(i) = pi*((D/2)^2);
        otherwise     %其它阶段
            s(i) = s(i - 1);
            Ap(i) = Ap(i - 1);
    end
    %计算其他参数
    Ab(i) = s(i)*Lp;
    Vg(i) = Ap(i)*Lp;
    m_b(i) = rho_p*Ab(i)*rb(i);
    m_p(i) = phi_m*p(i)*At / c;
    F(i) = 2000*m_p(i) + 4.4*At*(p(i)-p0);
    m_z(i) = m_c + (Vc - Vg(i))*rho_p;
    a(i) = (F(i) - f(i - 1)) / m_z(i) ;
    V(i) = V(i - 1) + a(i)*dt;
    
    %计算阻力系数
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
    
    %记录最大值
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
    
    %记录最小值
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
    
    %数组长度约束
    if i >= long
        break;
    end
    
end

%循环后处理

i = i + 1;
n(i:1:long) = [];

F(i:1:long) = [];

V(i:1:long) = [];

t = dt*n;
t_max = dt*(i - 1);     %燃烧总时间



str = [' 多孔燃烧总时间： ',num2str(t_max),'s'];
disp(str);

%数据暂存

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


%双阶段装药数据计算
%计算辅助常量

n = 1:1:long;       %绘图横坐标


%计算得常量
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数
beta_s = pi / n_s;        %等分角(rad)
% v_f = v_f*Ma;
% c_f = c_f*1e6;
% C_x = c_f(1);

%计算中间值
%星孔装药
%燃烧周长项中参数
s_sa = (sin(epsilon_s*beta_s)) / sin(theta_s) + (1 - epsilon_s)*beta_s;
s_sb = (pi / 2 + beta_s - theta_s - cot(theta_s));
s_sc = (1 - epsilon_s)*beta_s;
s_sd =  (l_s*sin(epsilon_s*beta_s));
%通气面积项中参数
Ap_sa = ((l_s^2) / 2)*( (1 - epsilon_s)*beta_s + (sin(epsilon_s*beta_s)) ...
    *( (cos(epsilon_s*beta_s)) - (sin(epsilon_s*beta_s))*cot(theta_s) ) );
Ap_sb = s_sa;
Ap_sc = s_sb;
Ap_sd = cot(theta_s) + theta_s - pi / 2;
Ap_se = s_sc;
Ap_sf = s_sd;
Ap_sg = epsilon_s*beta_s;
%压强项中参数
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;

%初始参数
%星孔装药
%周长参数
s_s_0 = 2*n_s*(l_s*s_sa + (r_s + r1_s)*s_sb  - beta_s*r1_s);      
%通气面积参数
Ap_si0 = Ap_sa + l_s*r_s*Ap_sb + ((r_s^2) / 2)*Ap_sc + ((r1_s^2) / 2)*Ap_sd;
Ap_s_0 = 2*n_s*Ap_si0;
%圆柱装药
s_c_0 = pi*dc;
Ap_c_0 = pi*(dc^2) / 4;

Vp_s0 = (pi*(Dr^2) / 4 - Ap_s_0)*Lp_s;       %星孔药柱体积(m^3)
Vp_c0 = (pi*(Dr^2) / 4 - Ap_c_0)*Lp_c;       %圆柱药柱体积(m^3)
Vp0 = Vp_s0 + Vp_c0;      %总药柱体积(m^3)
mp_s0 = Vp_s0*rho_p;     %星孔药柱质量(kg)
mp_c0 = Vp_c0*rho_p;     %圆柱药柱质量(kg)
mp = mp_s0 + mp_c0;     %总药柱质量(kg)

Vc = (pi*(Dr^2) / 4)*Lp;        %腔室体积(m^3)
Vps0 = Vc - Vp0;       %药柱初始体积(m^3)

%约束条件
%参数约束条件

%分阶段判别条件
e_a = r1_s;     %星根半角消失
e_b = l_s*(sin(epsilon_s*beta_s)) / cos(theta_s) - r_s;     %星边消失

%各阶段对应条件
%     星孔第一阶段，圆柱第一阶段
%     e11 = ((e <= e_a) & (e <= e_b) & (e <= ep_s) & (e <= ep_c));
%     sw = 1;
%     星孔第二阶段，圆柱第一阶段
%     e21 = ((e > e_a) & (e < e_b) & (e <= ep_s) & (e <= ep_c));
%     sw = 3;
%     星孔第三阶段，圆柱第一阶段
%     e31 = ((e > e_a) & (e > e_b) & (e <= ep_s) & (e <= ep_c));
%     sw = 7;
%     星孔结束，圆柱第一阶段
%     e41 = ((e > e_a) & (e > e_b) & (e > ep_s) & (e <= ep_c));
%     sw = 15;
%     星孔第一阶段，圆柱结束
%     e12 = ((e <= e_a) & (e <= e_b) & (e <= ep_s) & (e > ep_c));
%     sw = 17;
%     星孔第二阶段，圆柱结束
%     e22 = ((e > e_a) & (e < e_b) & (e <= ep_s) & (e > ep_c));
%     sw = 19;
%     星孔第三阶段，圆柱结束
%     e32 = ((e > e_a) & (e > e_b) & (e <= ep_s) & (e > ep_c));
%     sw = 23;
%程序结束条件
ep = max(ep_s,ep_c);
%     e >= ep;


%变量
%计算用变量数组（已赋初值）
p = p0*ones(1,long);        %实际压强(Pa)
rb = rb_0*ones(1,long);     %燃速(m/s)
e = 0*ones(1,long);     %已烧去肉厚(m)
s_s = s_s_0*ones(1,long);     %星孔燃烧面实际边长(m)
s_c = s_c_0*ones(1,long);     %圆柱燃烧面实际边长(m)
m_b = rho_p*(s_s_0*Lp_s + s_c_0*Lp_c)*rb_0*ones(1,long);     %燃气生成率(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %质量流率(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %推力(N)
Ab = (s_s_0*Lp_s + s_c_0*Lp_c)*ones(1,long);     %燃烧面积(m^2)
Ap_s = Ap_s_0*ones(1,long);     %星孔通气面积(m^2)
Ap_c = Ap_c_0*ones(1,long);     %圆柱通气面积(m^2)
Vg = (Ap_s*Lp_s + Ap_c*Lp_c);     %自由体积(m^3)
m_z = (m_c + mp)*ones(1,long);      %总质量(kg)
a = ( 2000*(p0*At / c) / (m_c + mp))*ones(1,long);     %加速度(N/m^2)
V = V0*ones(1,long);      %速度(m/s)
f = zeros(1,long);      %阻力(N)
%最大值变量
p_max = p0;        %最大压强(Pa)
e_max = 0;     %最大已烧去肉厚(m)
s_s_max = s_s_0;     %星孔燃烧面最大周长(m)
s_c_max = s_c_0;     %圆柱燃烧面最大周长(m)
m_b_max = rho_p*(s_s_0*Lp_s + s_c_0*Lp_c)*rb_0;     %最大燃气生成率(kg/s)
m_p_max = (phi_m*p0*At / c);     %最大质量流率(kg/s)
F_max = 2000*(p0*At / c);     %最大推力(N)
Ab_max = (s_s_0*Lp_s + s_c_0*Lp_c);     %最大燃烧面积(m^2)
Ap_s_max = Ap_s_0;     %星孔最大通气面积(m^2)
Ap_c_max = Ap_c_0;     %圆柱最大通气面积(m^2)
Vg_max = (Ap_s*Lp_s + Ap_c*Lp_c);     %最大自由体积(m^3)
m_z_max = (m_c + mp);
a_max = ( 2000*(p0*At / c) / (m_c + mp));
V_max = V0;
f_max = 0;
%最小值变量
p_min = p0;        %最小压强(Pa)
e_min = 0;     %最小已烧去肉厚(m)
s_s_min = s_s_0;     %燃烧面最小周长(m)
s_c_min = s_c_0;     %燃烧面最小周长(m)
m_b_min = rho_p*(s_s_0*Lp_s + s_c_0*Lp_c)*rb_0;     %最小燃气生成率(kg/s)
m_p_min = (phi_m*p0*At / c);     %最小质量流率(kg/s)
F_min = 2000*(p0*At / c);     %最小推力(N)
Ab_min = (s_s_0*Lp_s + s_c_0*Lp_c);     %最小燃烧面积(m^2)
Ap_s_min = Ap_s_0;     %星孔最小通气面积(m^2)
Ap_c_min = Ap_c_0;     %圆柱最小通气面积(m^2)
Vg_min = (Ap_s*Lp_s + Ap_c*Lp_c);     %最小自由体积(m^3)
m_z_min = (m_c + mp);
a_min = ( 2000*(p0*At / c) / (m_c + mp));
V_min = V0;
f_min = 0;
%循环变量
i = 1;
j = 1;
k = 1;
sw = zeros(1,long);        %阶段选择
swc = zeros(2,100);        %阶段变化点

%计算中间值

%各阶段燃烧周长
%     ssi1 = l_s*s_sa + (r_s + r1_s)*s_sb - beta_s*(r1_s - e);
%     ssi2 = l_s*s_sa + (r_s + e)*s_sb;
%     ssi3 = l_s*s_sc + (r_s + e)*( beta_s + asin(s_sd / (r_s + e)) );
%     s_s = 2*n_s*ssi;
%各阶段通气面积
%     Apsi1 = Ap_sa + l_s*(r_s + e)*Ap_sb + (((r_s + e)^2) / 2)*Ap_sc ...
%        + (((r1_s - e)^2) / 2)*Ap_sd;
%     Apsi2 = Ap_sa + l_s*(r_s + e)*Ap_sb + (((r_s + e)^2) / 2)*Ap_sc;
%     Apsi3 = ( ((l_s + r_s + e)^2)*Ap_se ...
%        + ((r_s + e)^2)*(Ap_sg + asin(Ap_sf / (r_s + e))) ...
%        + Ap_sf*(sqrt((r_s + e)^2 - Ap_sf^2)  + l_s*cos(Ap_sg)) ) / 2;
%     Ap_s = 2*n_s*Apsi;

%龙格库塔计算公式
%     dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;

%数值计算
%循环前常参数计算



%循环计算
while (e(i) <= ep)
    %龙格库塔逐步计算，计算出新的压强值
    %dp = p_a*(Ab / Vg)*p^n_p  - p_b*p / Vg;
    k1 = p_a*(Ab(i) / Vg(i))*( p(i)^n_p ) - p_b*p(i) / Vg(i);
    k2 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k1/2)^n_p ) - p_b*(p(i) + dt*k1/2) / Vg(i);
    k3 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k2/2)^n_p ) - p_b*(p(i) + dt*k2/2) / Vg(i);
    k4 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k3)^n_p )  - p_b*(p(i) + dt*k3) / Vg(i);
    dp = (k1 + 2*k2 + 2*k3 + k4)/6;
    i = i + 1;
    p(i) = p(i - 1) + dp*dt;

    %参数修正
    %逐项计算其他参数的新值
    rb(i) = alpha_r*p(i)^n_p;
    e(i) = e(i - 1) + rb(i)*dt;

    %判断此时处于哪个阶段
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
    %记录阶段改变时的循环数值
    if(sw(i) ~= sw(i - 1))
        swc(1,j) = i;
        swc(2,j) = sw(i);
        j = j + 1;
    end
    %分阶段计算燃烧周长和通气面积
    switch sw(i)
        case{1}     %星孔第一阶段，圆柱第一阶段
            s_s(i) = 2*n_s*(l_s*s_sa + (r_s + r1_s)*s_sb - beta_s*(r1_s - e(i)));
            Apsi1 = Ap_sa + l_s*(r_s + e(i))*Ap_sb + (((r_s + e(i))^2) / 2)*Ap_sc ...
                + (((r1_s - e(i))^2) / 2)*Ap_sd;
            Ap_s(i) = 2*n_s*Apsi1;
            s_c(i) = pi*(dc + 2*e(i));
            Ap_c(i) = pi*(dc + 2*e(i))^2 / 4;
        case{3}     %星孔第二阶段，圆柱第一阶段
            s_s(i) = 2*n_s*(l_s*s_sa + (r_s + e(i))*s_sb);
            Apsi2 = Ap_sa + l_s*(r_s + e(i))*Ap_sb + (((r_s + e(i))^2) / 2)*Ap_sc;
            Ap_s(i) = 2*n_s*Apsi2;
            s_c(i) = pi*(dc + 2*e(i));
            Ap_c(i) = pi*(dc + 2*e(i))^2 / 4;
        case{7}     %星孔第三阶段，圆柱第一阶段
            s_s(i) = 2*n_s*(l_s*s_sc + (r_s + e(i))*( beta_s + asin(s_sd / (r_s + e(i))) ) );
            Apsi3 = ( ((l_s + r_s + e(i))^2)*Ap_se ...
                + ((r_s + e(i))^2)*(Ap_sg + asin(Ap_sf / (r_s + e(i)))) ...
                + Ap_sf*(sqrt((r_s + e(i))^2 - Ap_sf^2)  + l_s*cos(Ap_sg)) ) / 2;
            Ap_s(i) = 2*n_s*Apsi3;
            s_c(i) = pi*(dc + 2*e(i));
            Ap_c(i) = pi*(dc + 2*e(i))^2 / 4;
        case{15}     %星孔结束，圆柱第一阶段
            s_s(i) = 0;
            if(Ap_s(i - 1 ) < (pi*Dr^2 / 4) )
                Ap_s(i) = Ap_s(i - 1) +0.0008*Apsi3;
            else
                Ap_s(i) = pi*Dr^2 / 4;
            end
            s_c(i) = pi*(dc + 2*e(i));
            Ap_c(i) = pi*(dc + 2*e(i))^2 / 4;
        case{17}     %星孔第一阶段，圆柱结束
            s_s(i) = 2*n_s*(l_s*s_sa + (r_s + r1_s)*s_sb - beta_s*(r1_s - e(i)));
            Apsi1 = Ap_sa + l_s*(r_s + e(i))*Ap_sb + (((r_s + e(i))^2) / 2)*Ap_sc ...
                + (((r1_s - e(i))^2) / 2)*Ap_sd;
            Ap_s(i) = 2*n_s*Apsi1;
            s_c(i) = 0;
            Ap_c(i) = pi*Dr^2 / 4;
        case{19}     %星孔第二阶段，圆柱结束
            s_s(i) = 2*n_s*(l_s*s_sa + (r_s + e(i))*s_sb);
            Apsi2 = Ap_sa + l_s*(r_s + e(i))*Ap_sb + (((r_s + e(i))^2) / 2)*Ap_sc;
            Ap_s(i) = 2*n_s*Apsi2;
            s_c(i) = 0;
            Ap_c(i) = pi*Dr^2 / 4;
        case{23}     %星孔第三阶段，圆柱结束
            s_s(i) = 2*n_s*(l_s*s_sc + (r_s + e(i))*( beta_s + asin(s_sd / (r_s + e(i))) ) );
            Apsi3 = ( ((l_s + r_s + e(i))^2)*Ap_se ...
                + ((r_s + e(i))^2)*(Ap_sg + asin(Ap_sf / (r_s + e(i)))) ...
                + Ap_sf*(sqrt((r_s + e(i))^2 - Ap_sf^2)  + l_s*cos(Ap_sg)) ) / 2;
            Ap_s(i) = 2*n_s*Apsi3;
            s_c(i) = 0;
            Ap_c(i) = pi*Dr^2 / 4;
        otherwise     %其它阶段
            s_s(i) = s_s(i - 1);
            Ap_s(i) = Ap_s(i - 1);
            s_c(i) = s_c(i - 1);
            Ap_c(i) = Ap_c(i - 1);
    end
    %计算其他参数
    Ab(i) = s_s(i)*Lp_s + s_c(i)*Lp_c;
    Vg(i) = Ap_s(i)*Lp_s + Ap_c(i)*Lp_c;
    m_b(i) = rho_p*Ab(i)*rb(i);
    m_p(i) = phi_m*p(i)*At / c;
    F(i) = 2000*m_p(i) + 4.4*At*(p(i)-p0);
    m_z(i) = m_c + (Vc - Vg(i))*rho_p;
    a(i) = (F(i) - f(i - 1)) / m_z(i) ;
    V(i) = V(i - 1) + a(i)*dt;
    
    %计算阻力系数
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
    
    %记录最大值
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
    
    %记录最小值
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
    
    %数组长度约束
    if i >= long
        break;
    end
    
end
%循环后处理

i = i + 1;
n(i:1:long) = [];

F(i:1:long) = [];
V(i:1:long) = [];

t = dt*n;
t_max = dt*(i - 1);     %燃烧总时间
prx = 1.05;
pri = -0.05;



str = [' 燃烧总时间： ',num2str(t_max),'s'];
disp(str);

%数据暂存

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

%数据处理
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


%数据输出


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
% title('速度');
xlabel('t/s');
ylabel('V/(m/s)');
legend('case7','case8');


disp('done');
