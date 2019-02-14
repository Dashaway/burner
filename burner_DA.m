%181016
%Dashaway
%燃烧参数计算


%190212
%换为跑道型内孔
%更改参数

%190213
%更改圆弧角度


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
Dr = 0.150;       %燃烧室外径(m)
D = Dr;         %药柱外径(m)
r = 0.003;        %微圆弧(m)
l = 0.055;        %中心距(m)
R = 0.018;       %中心孔(m)
m_s = 2;        %孔药比
n_s = 4;    %药柱数量
Lp = 0.22;      %装药长度(m)
gamma = 1.2;        %比热比
rho_p = 1730;       %密度(kg/m^3)
n_p = 0.302;      %压强指数
c = 1600;       %特征速度(m/s)
At = 0.884785e-3;     %喷管喉部面积(m^2)
rb_0 = 5e-3;      %?初始燃速(m/s)
alpha_r = 1.7e-4;        %?燃速系数
phi_alpha = 1;       %?侵蚀函数
phi_m = 1;      %?

%计算得常量

Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数
trd = (2*pi) / (n_s*(m_s + 1));     %药圆心角
trr = (2*pi*m_s) / (n_s*(m_s + 1));     %弧圆心角
%判别条件
e_a = l*sin(trd/2) - r;       %第一阶段
e_b = l*sin(trr/2) - r;
e_c = D/2 - r - l;
e_d = (l - r - R)/2;
e_e = (D/2)*(sin(trd) / sin(trd + as2)) - r;

%计算中间值
as1 =  asin((D/2 - l) / (r + e));
as2 =  asin( (l*sin(pi / (2*n_s))) / (r + e) );
ac1 = acos( (4*l^2 + D^2 - 4*(r + e)^2) / 4*l*D ); 
ac2 = acos( (4*(r + e)^2 + 4*l^2 - D^2) / (8*l*(r + e)) ); 
%         as =  asin( (l*sin(pi / (2*n_s))) / (r + e(i)) );
%         ac1 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) ); 
%         ac2 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / 4*l*D ); 

%周长参数
s_a0 = (4*pi*l*m_s) / (m_s + 1) + 2*pi*(R + n_s*r);      %第一阶段燃烧面初始边长
    s_a =  (4*pi*l*m_s) / (m_s + 1) + 2*pi*((R + e) + n_s*(r + e));    %第一阶段燃烧周长公式
    s_b = 2*n_s*(as1 + pi/2)*(r + e) + 2*n_s*(trr/2)*(l + R - r) ;  %第二阶段燃烧周长公式
    s_c = 2*n_s*(r + e)*(ac2 + 2*as2 - pi) + 2*pi*R + 2*n_s*trr*(l - r);  %第三阶段燃烧周长公式
    s_d = 2*n_s*(r + e)*(as2 - trd/2) + 2*pi*(R + e) + 2*n_s*trr*(l - r - e);  %第四阶段燃烧周长公式
%通气面积参数
Ap_a0 = (4*pi*l*m_s*r) / (m_s + 1) + n_s*pi*r^2 + pi*R^2;     %第一阶段初始通气面积
    Ap_a = (4*pi*l*m_s*(r + e)) / (m_s + 1) ... %第一阶段通气面积
        + n_s*pi*(r + e)^2 + pi*(R + e)^2;   
    Ap_b = n_s*(trr/2)*(R + e)^2 + n_s*(as1 + pi/2)*(r + e)^2 ...  %第二阶段通气面积
        + n_s*D*ac1 - n_s*l*(D/2) + 2*n_s*l*(r + e)*trr;
    %Ap_c = pi*(D^2)/4 - pi*(R + e)^2 - n_s*D*(pi/2 -ac2) ...%第三阶段通气面积
        %*(  D/2 - (r + e)*sin(as + pi / (2*n_s)) / sin(pi / (2*n_s))  );
    Ap_c = pi*(D^2)/4 - pi*(R + e)^2;  %第三阶段通气面积
    Ap_d = pi*(D^2)/4 - pi*(R + e)^2;  %第四阶段通气面积
    
        

%变量
%计算用变量（已赋初值）
p = p0*ones(1,long);        %实际压强(Pa)
rb = rb_0*ones(1,long);     %燃速(m/s)
e = 0*ones(1,long);     %已烧去肉厚(m)
s = s_a0*ones(1,long);     %燃烧面实际边长(m)
m_b = rho_p*s_a0*Lp*rb_0*ones(1,long);     %燃气生成率(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %质量流率(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %推力(?)
Ab = s_a0*Lp*ones(1,long);     %燃烧面积(m^2)
Ap = Ap_a0*ones(1,long);     %通气面积(m^2)
Vg = Ap*Lp;     %自由体积(m^3)
%循环变量
i = 1;
%状态变量
SW = 0;     %阶段选择



%数值计算
%循环前常参数计算
%压强项中参数
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;


%循环计算
while e <= e_c
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
    rb(i) = alpha_r*p(i)^n_p;
    e(i) = e(i - 1) + rb(i)*dt;
    %分阶段计算燃烧周长和通气面积
    if(e < e_a)
        s(i) = 2*pi*(R + l + e(i)) + 2*pi*n_s*(r + e(i));
        Ap(i) = 2*pi*l*(r + e(i)) + n_s*pi*(r + e(i))^2 + pi*(R + e(i))^2;
    elseif(e < e_b)
        as =  asin( (l*sin(pi / (2*n_s))) / (r + e(i)) );
        s(i) = 4*n_s*(r + e(i))*as + 2*pi*(l + R + e(i));
        Ap(i) = 2*n_s*((r + e(i))^2)*as + 2*pi*l*(r + e(i)) + pi*(R + e(i))^2 ...
            + n_s*((r + e(i))^2)*sin(2*as) ;
    else
        as =  asin( (l*sin(pi / (2*n_s))) / (r + e(i)) );
        ac1 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) ); 
        ac2 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / 4*l*D ); 
        s(i) = 2*n_s*(r + e(i))*(ac1 + 2*as - pi) + 2*pi*(l + R - r);
        Ap(i) = pi*(D^2)/4 - pi*(R + e(i))^2 - n_s*D*(pi/2 -ac2) ...
            *( D/2 - (r + e(i))*sin(as + pi / (2*n_s)) / sin(pi / (2*n_s))  );
    end
    
    Ab(i) = s(i)*Lp;
    Vg(i) = Ap(i)*Lp;
    m_b(i) = rho_p*Ab(i)*rb(i);
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

