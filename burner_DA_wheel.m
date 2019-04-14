%181016
%Dashaway
%燃烧参数计算


%190327
%轮式装药


clear;
close all;
%初始值设定


%常量

%计算辅助常量
long = 150000;        %数组长度
n = 1:1:long;       %绘图横坐标
dt = 5e-5;      %步进值
%燃烧室参数
Dr = 0.240;       %燃烧室外径(m)
% At = 1.884785e-3;     %喷管喉部面积(m^2)
At = 7.853982e-3;     %喷管喉部面积(m^2)
%装药参数
%轮式装药参数
Dc = 0.238;     %外径(m)
dc = 0.0868;      %内径(m)
e_s1 = 0.0155;       %厚肉厚(m)
e_s2 = 0.0033;       %薄肉厚(m)
n_s = 12;       %轮孔数
l_s1 = 0.0995;      %内圆弧特征长度(m)
l_s2 = 0.053;       %外圆弧特征长度(m)
r_s1 = 0.004;       %上过渡圆弧半径(m)
r_s2 = 0.003;       %下过渡圆弧半径(m)
Lp = 0.345;      %药柱长度(m)
%推进剂参数
n_p = 0.302;      %压强指数
rho_p = 1730;       %密度(kg/m^3)
c = 1600;       %特征速度(m/s)
rb_0 = 5e-3;      %?初始燃速(m/s)
alpha_r = 1.7e-4;        %?燃速系数(m/s*MPa)

%已知常量
p0 = 1.02e5;    %初始压强(Pa)
gamma = 1.2;        %比热比
phi_alpha = 1;       %?侵蚀函数
phi_m = 1;      %?

%计算得常量
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数
beta_s = 2*pi / n_s;    %轮孔角(rad)

%计算中间值
as1 = asin((e_s1 + r_s1) / l_s1);       %厚轮辐上部角分数
as2 = asin((e_s1 + r_s2) / l_s2);       %厚轮辐下部角分数
as3 = asin((e_s2 + r_s1) / l_s1);       %薄轮辐上部角分数
as4 = asin((e_s2 + r_s2) / l_s2);       %薄轮辐下部角分数
h_s1 = l_s1*cos(as1) - l_s2*cos(as2);       %厚轮辐高度
h_s2 = l_s1*cos(as3) - l_s2*cos(as4);       %薄轮辐高度
%压强项中参数
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;
%各阶段燃烧周长
% si1 = (l_s1 + r_s1 + e)*(beta_s - as1 - as3) + (r_s1 + e)*(pi + as1 + as2) ...
%     + h_s1 + h_s2 + (r_s2 + e)*(pi - as2 - as4) ...
%     + (l_s2 - r_s2 - e)*(beta_s - as2 - as4) + (e + dc/2)*beta_s;
% si2 = (l_s1 + r_s1 + e)*(beta_s - as1 - as3) + (r_s1 + e)*(pi/2 + as1) ...
%     + l_s1*cos(as1) - (e + dc/2) ...
%     + (e + dc/2)*asin((l_s2*sin(as2) - r_s2 - e) / (e + dc/2)) ...
%     + (r_s1 + e)*(as3 + asin( l_s1*sin(as3) / (r_s1 + e) ));
% s = n_s*si*Lp;
%各阶段通气面积
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


%初始参数
%圆柱装药
%周长参数
si0 = (l_s1 + r_s1)*(beta_s - as1 - as3) + r_s1*(pi + as1 + as2) ...
    + h_s1 + h_s2 + r_s2*(pi - as2 - as4) ...
    + (l_s2 - r_s2)*(beta_s - as2 - as4) + (dc/2)*beta_s;
s_0 = n_s*si0*Lp;
%通气面积参数
Api0 = (l_s1^2)*((sin(2*as1) + sin(2*as2)) / 4) ...
    - (l_s2^2)*((sin(2*as2) + sin(2*as4)) / 4) ...
    + ((l_s1 + r_s1)^2)*((beta_s - as1 - as3) / 2) ...
    - ((l_s2 - r_s2)^2)*((beta_s - as2 - as4) / 2) ...
    + (r_s1^2)*((pi + as1 + as3) / 2) ...
    + (r_s2^2)*((pi - as2 - as4) / 2) ...
    - h_s1*e_s1 - h_s2*e_s2 + ((dc/2)^2)*beta_s/2;
Ap_0 = n_s*Api0;

Vp0 = (pi*(Dr^2) / 4 - Ap_0)*Lp;       %药柱体积(m^3)
mp = Vp0*rho_p;     %药柱质量(kg)
%压强项


%约束条件
%参数约束条件

%分阶段判别条件


%各阶段对应条件
%     大推力阶段
%     e1 = ((e <= e_s1) & (e <= e_s2));
%     sw = 1;
%     小推力阶段
%     e2 = ((e <= e_s1) & (e > e_s2));
%     sw = 3;
%程序结束条件
ep = e_s1;
%     e >= ep;


%变量
%计算用变量数组（已赋初值）
p = p0*ones(1,long);        %实际压强(Pa)
rb = rb_0*ones(1,long);     %燃速(m/s)
e = 0*ones(1,long);     %已烧去肉厚(m)
s = s_0*ones(1,long);     %燃烧面实际周长(m)
m_b = rho_p*s_0*Lp*rb_0*ones(1,long);     %燃气生成率(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %质量流率(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %推力(N)
Ab = s_0*Lp*ones(1,long);     %燃烧面积(m^2)
Ap = Ap_0*ones(1,long);     %通气面积(m^2)
Vg = Ap*Lp;     %自由体积(m^3)
%最大值变量
p_max = p0;        %最大压强(Pa)
e_max = 0;     %最大已烧去肉厚(m)
s_max = s_0;     %燃烧面最大周长(m)
m_b_max = rho_p*s_0*Lp*rb_0;     %最大燃气生成率(kg/s)
m_p_max = (phi_m*p0*At / c);     %最大质量流率(kg/s)
F_max = 2000*(p0*At / c);     %最大推力(N)
Ab_max = s_0*Lp;     %最大燃烧面积(m^2)
Ap_max = Ap_0;     %最大通气面积(m^2)
Vg_max = Ap_0*Lp;     %最大自由体积(m^3)
%最小值变量
p_min = p0;        %最小压强(Pa)
e_min = 0;     %最小已烧去肉厚(m)
s_min = s_0;     %燃烧面最小周长(m)
m_b_min = rho_p*s_0*Lp*rb_0;     %最小燃气生成率(kg/s)
m_p_min = (phi_m*p0*At / c);     %最小质量流率(kg/s)
F_min = 2000*(p0*At / c);     %最小推力(N)
Ab_min = s_0*Lp;     %最小燃烧面积(m^2)
Ap_min = Ap_0;     %圆柱最小通气面积(m^2)
Vg_min = Ap_0*Lp;     %最小自由体积(m^3)
%循环变量
i = 1;
j = 1;
sw = zeros(1,long);        %阶段选择
swc = zeros(2,100);        %阶段变化点

%计算中间值

%各阶段燃烧周长
%各阶段通气面积


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
    if(e(i) <= e_s2)
        sw(i) = sw(i) + 0*2;
    else
        sw(i) = sw(i) + 1*2;
    end
    %记录阶段改变时的循环数值
    if(sw(i) ~= sw(i - 1))
        swc(1,j) = i;
        swc(2,j) = sw(i);
        j = j + 1;
    end
    %分阶段计算燃烧周长和通气面积
    switch sw(i)
        case{1}     %大推力阶段
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
        case{3}     %小推力阶段
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
    
    %数组长度约束
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
sw(i:1:long) = [];
swc(:,j:100) = [];

t = dt*n;
t_max = dt*(i - 1);     %燃烧总时间
prx = 1.05;
pri = -0.05;

str = [' 燃烧总时间： ',num2str(t_max),'s'];
disp(str);

%数据输出
%写入文档


%绘制图像
figure;
hold on;
subplot(2,1,1);
plot(t,p);
axis ([pri*t_max,prx*t_max,(pri*(p_max - p_min) + p_min), ...
    (prx*(p_max - p_min) + p_min)]);
title('燃烧室压力');
xlabel('时间(s)');
ylabel('压强(Pa)');
% legend('算例2');

subplot(2,1,2);
plot(t,e);
axis ([pri*t_max,prx*t_max,(pri*(e_max - e_min) + e_min), ...
    (prx*(e_max - e_min) + e_min)]);
title('烧去肉厚');
xlabel('时间(s)');
ylabel('长度(m)');
% legend('算例2');

figure;
hold on;
subplot(2,1,1);
plot(t,Ab);
axis ([pri*t_max,prx*t_max,(pri*(Ab_max - Ab_min) + Ab_min), ...
    (prx*(Ab_max - Ab_min) + Ab_min)]);
title('燃烧面积');
xlabel('时间(s)');
ylabel('面积(m^2)');
% legend('算例2');

subplot(2,1,2);
plot(t,Ap);
axis ([pri*t_max,prx*t_max,(pri*(Ap_max - Ap_min) + Ap_min), ...
    (prx*(Ap_max - Ap_min) + Ap_min)]);
title('通气面积');
xlabel('时间(s)');
ylabel('面积(m^2)');
% legend('算例2');

figure;
hold on;
subplot(2,1,1);
plot(t,m_b);
axis ([pri*t_max,prx*t_max,(pri*(m_b_max - m_b_min) + m_b_min), ...
    (prx*(m_b_max - m_b_min) + m_b_min)]);
title('燃气生成率');
xlabel('时间(s)');
ylabel('质量流率(kg/s)');
% legend('算例2');

subplot(2,1,2);
plot(t,m_p);
axis ([pri*t_max,prx*t_max,(pri*(m_p_max - m_p_min) + m_p_min), ...
    (prx*(m_p_max - m_p_min) + m_p_min)]);
title('质量流率');
xlabel('时间(s)');
ylabel('质量流率(kg/s)');
% legend('算例2');

figure;
hold on;
plot(t,F);
axis ([pri*t_max,prx*t_max,(pri*(F_max - F_min) + F_min), ...
    (prx*(F_max - F_min) + F_min)]);
title('推力');
xlabel('时间(s)');
ylabel('力(N)');
legend('算例2');

figure;
hold on;
plot(t,p);
axis ([pri*t_max,prx*t_max,(pri*(p_max - p_min) + p_min), ...
    (prx*(p_max - p_min) + p_min)]);
title('燃烧室压力');
xlabel('时间(s)');
ylabel('压强(Pa)');
legend('算例2');


%结束
