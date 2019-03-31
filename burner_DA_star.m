%181016
%Dashaway
%燃烧参数计算
%星型装药，等面燃烧

%181114
%加入燃气生成率和质量流率曲线
%190117
%增加注释

clear;
close all;
%初始值设定
%常量
%辅助常量
long = 50000;        %数组长度
n = 1:1:long;       %绘图横坐标
dt = 5e-5;      %步进值
%已知常量
p0 = 1.02e5;    %初始压强(Pa)
D = 0.132;       %外径(m)
ep = 0.04;     %总肉厚(m)
r_s = 0.003;        %星尖圆弧半径(m)
r1_s = 0.003;        %星根过渡圆弧半径(m)
l_s = 0.023;        %特征长度(m)
n_s = 7;        %星角数
beta_s = pi / n_s;        %等分角(rad)
epsilon_s = 1;      %角分数
theta_s = 0.620465;      %星根半角(rad)
Lp = 0.22;      %装药长度(m)
gamma = 1.2;        %比热比
rho_p = 1730;       %密度(kg/m^3)
n_p = 0.302;      %压强指数
c = 1600;       %特征速度(m/s)
At = 5e-4;     %喷管喉部面积(m^2)
r0 = 5e-3;      %?初始燃速(m/s)
alpha_r = 1.7e-4;        %?燃速系数
phi_alpha = 1;       %?侵蚀函数
phi_m = 1;      %?

%计算得常量
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数
%周长参数
s_a = (sin(epsilon_s*beta_s)) / sin(theta_s) + (1 - epsilon_s)*beta_s;
s_b = (pi / 2 + beta_s - theta_s - cot(theta_s));
s_c = (1 - epsilon_s)*beta_s;
s_d =  (l_s*sin(epsilon_s*beta_s));
s0 = 2*n_s*(l_s*s_a + (r_s + r1_s)*s_b  - beta_s*r1_s);      %燃烧面初始边长
%si1 = l_s*s_a + (r_s + r1_s)*s_b - beta_s*(r1_s - e);
%si2 = l_s*s_a + (r_s + e)*s_b;
%si3 = l_s*s_c + (r_s + e)*( beta_s + asin(s_d / (r_s + e)) );
%通气面积参数
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

Vp0 = (pi*(D^2) / 4 - Ap0)*Lp;       %药柱体积(m^3)
mp = Vp0*rho_p;     %药柱质量(kg)


%变量
%计算用变量（已赋初值）
p = p0*ones(1,long);        %实际压强(Pa)
r = r0*ones(1,long);     %燃速(m/s)
e = 0*ones(1,long);     %已烧去肉厚(m)
s = s0*ones(1,long);     %燃烧面边长(m)
m_b = rho_p*s0*Lp*r0*ones(1,long);     %燃气生成率(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %质量流率(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %推力(?)
Ab = s0*Lp*ones(1,long);     %燃烧面积(m^2)
Ap = Ap0*ones(1,long);     %通气面积(m^2)
Vg = Ap*Lp;     %自由体积(m^3)
%循环变量
i = 1;


%数值计算
%循环前常参数计算
%压强项中参数
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;

e_a = r1_s;
e_b = l_s*(sin(epsilon_s*beta_s)) / cos(theta_s) - r_s;

%循环计算
while (e(i) <= ep)
    %龙格库塔逐步计算，计算出新的压强值
    %dp = p_a*(Ab(i) / Vg(i))*p(i)^n_p  - p_b*p(i) / Vg(i);
    k1 = p_a*(Ab(i) / Vg(i))*( p(i)^n_p ) - p_b*p(i) / Vg(i);
    k2 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k1/2)^n_p ) - p_b*(p(i) + dt*k1/2) / Vg(i);
    k3 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k2/2)^n_p ) - p_b*(p(i) + dt*k2/2) / Vg(i);
    k4 = p_a*(Ab(i) / Vg(i))*( (p(i) + dt*k3)^n_p )  - p_b*(p(i) + dt*k3) / Vg(i);
    dp = (k1 + 2*k2 + 2*k3 + k4)/6;
    i = i + 1;
    p(i) = p(i - 1) + dp*dt;

    %参数修正
    %逐项计算其他参数的新值
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

%循环后处理
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

%数据输出

figure;
hold on;
subplot(2,1,1);
plot(t,p);
grid on,axis tight ;
title('燃烧室压力')

subplot(2,1,2);
plot(t,e);
grid on,axis tight ;
title('烧去肉厚')

figure;
hold on;
subplot(2,1,1);
plot(t,Ab);
grid on,axis tight ;
title('燃烧面积')

subplot(2,1,2);
plot(t,Ap);
grid on,axis tight ;
title('通气面积')

figure;
hold on;
subplot(2,1,1);
plot(t,m_b);
grid on,axis tight ;
title('燃气生成率')

subplot(2,1,2);
plot(t,m_p);
grid on,axis tight ;
title('质量流率')

figure;
hold on;
plot(t,F);
grid on,axis tight ;
title('推力')

figure;
hold on;
plot(t,p);
grid on,axis tight ;
title('燃烧室压力')

