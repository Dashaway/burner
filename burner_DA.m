%181016
%Dashaway
%燃烧参数计算
%圆柱装药，等面燃烧

%190117
%增加注释
%换为圆柱药柱，药柱数量增加

clear;
close all;
%初始值设定
%常量
%辅助常量
long = 100000;        %数组长度
n = 1:1:long;       %绘图横坐标
dt = 5e-5;      %步进值
%已知常量
p0 = 1.02e5;    %初始压强(Pa)
Dr = 0.132;       %燃烧室外径(m)
D = 0.0265;         %药柱初始外径(m)
d = 0.0125;       %药柱初始内径(m)
ep = (D - d)/2;     %总肉厚(m)
n_s = 17;    %药柱数量
Lp = 0.22;      %装药长度(m)
gamma = 1.2;        %比热比
rho_p = 1730;       %密度(kg/m^3)
n_p = 0.302;      %压强指数
c = 1600;       %特征速度(m/s)
At = 1.884785e-3;     %喷管喉部面积(m^2)
r0 = 5e-3;      %?初始燃速(m/s)
alpha_r = 1.7e-4;        %?燃速系数
phi_alpha = 1;       %?侵蚀函数
phi_m = 1;      %?

%计算得常量

Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数
%周长参数
s0 = n_s*pi*d;      %燃烧面初始边长
    %s =  n_s*pi*(d + e);
%通气面积参数
Ap0 = pi*Dr^2 / 4 - n_s*pi*D^2 / 4 + n_s*pi*d^2 / 4;
    %Ap = pi*Dr^2 / 4 - n_s*pi*D^2 / 4 + n_s*pi*(d + e)^2 / 4;


%变量
%计算用变量（已赋初值）
p = p0*ones(1,long);        %实际压强(Pa)
r = r0*ones(1,long);     %燃速(m/s)
e = 0*ones(1,long);     %已烧去肉厚(m)
s = s0*ones(1,long);     %燃烧面实际边长(m)
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


%循环计算
while e <= ep
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
    
    s(i) = n_s*pi*(d + e(i));
    Ap(i) = pi*Dr^2 / 4 - n_s*pi*D^2 / 4 + n_s*pi*(d + e(i))^2 / 4;
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

