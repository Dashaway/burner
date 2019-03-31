%181016
%Dashaway
%燃烧参数计算


%190219
%整理公式
%完善注释

%190321
%分段燃烧


clear;
close all;
%初始值设定


%常量

%计算辅助常量
long = 150000;        %数组长度
n = 1:1:long;       %绘图横坐标
dt = 2e-5;      %步进值
%燃烧室参数
Dr = 0.082;       %燃烧室外径(m)
% At = 1.884785e-3;     %喷管喉部面积(m^2)
At = 1.103e-3;     %喷管喉部面积(m^2)
%装药参数
%星孔装药参数
Ds = 0.082;       %外径(m)
ep_s = 0.018;     %肉厚(m)
r_s = 0.0012;        %星尖圆弧半径(m)
r1_s = 0.001;        %星根过渡圆弧半径(m)
l_s = 0.022;        %特征长度(m)
n_s = 6;        %星角数
epsilon_s = 0.7;      %角分数
theta_s = pi*25/180;      %星根半角(rad)
%圆柱装药参数
Dc = 0.082;     %外径(m)
dc = 0.01;      %内径(m)
ep_c = (Dc - dc) / 2;       %肉厚(m)
%分段参数
Lp = 0.49;      %药柱长度(m)
pers = 0.7;     %星孔段占比
Lp_s = Lp*pers;     %星孔段长度
Lp_c = Lp - Lp_s;       %圆柱段长度
%已知常量
p0 = 1.02e5;    %初始压强(Pa)
gamma = 1.2;        %比热比
rho_p = 1700;       %密度(kg/m^3)
n_p = 0.4;      %压强指数
c = 1584;       %特征速度(m/s)
rb_0 = 9e-3;      %?初始燃速(m/s)
alpha_r = 9e-5;        %?燃速系数
phi_alpha = 1;       %?侵蚀函数
phi_m = 1;      %?

%计算得常量
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数
beta_s = pi / n_s;        %等分角(rad)

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
%循环变量
i = 1;
j = 1;
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
            Ap_s(i) = pi*Dr^2 / 4;
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
s_s(i:1:long) = [];
s_c(i:1:long) = [];
m_b(i:1:long) = [];
m_p(i:1:long) = [];
F(i:1:long) = [];
Ab(i:1:long) = [];
Ap_s(i:1:long) = [];
Ap_c(i:1:long) = [];
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


subplot(2,2,3);
plot(t,Ap_s);
axis ([pri*t_max,prx*t_max,(pri*(Ap_s_max - Ap_s_min) + Ap_s_min), ...
    (prx*(Ap_s_max - Ap_s_min) + Ap_s_min)]);
title('星孔通气面积');
xlabel('时间(s)');
ylabel('面积(m^2)');
subplot(2,2,4);
plot(t,Ap_c);
axis ([pri*t_max,prx*t_max,(pri*(Ap_c_max - Ap_c_min) + Ap_c_min), ...
    (prx*(Ap_c_max - Ap_c_min) + Ap_c_min)]);
title('圆柱通气面积');
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
% legend('算例2');

figure;
hold on;
plot(t,p);
axis ([pri*t_max,prx*t_max,(pri*(p_max - p_min) + p_min), ...
    (prx*(p_max - p_min) + p_min)]);
title('燃烧室压力');
xlabel('时间(s)');
ylabel('压强(Pa)');
% legend('算例2');


%双Y轴图像
%打开图，设置左右y轴属性
left_color = [0 0 0];
right_color = [0 0 0];
set(figure,'defaultAxesColorOrder',[left_color;right_color]);
hold on;
%激活左侧
yyaxis left;
plot(t,p,'-','LineWidth',1,'color','k');   
ylabel('压强(Pa)')
%设置刻度
per = 0.6;      %定比例
cou = p_max;        %选参数
mag = floor(log10(cou / per));        %求数量级
fir = floor((cou / per) / (10^mag));      %取第一位
sec = floor((cou / per - fir*10^mag) / (10^(mag - 1)));       %取第二位
ran = fir*10^mag + sec*10^(mag - 1);        %定量程
cal = fir*10^(mag - 1);     %定刻度
sta = -0.5*fir*10^(mag - 1);        %定起点
axis([pri*t_max,prx*t_max,sta,ran]);
set(gca,'YTick',0:cal:ran);
%激活右侧
yyaxis right;
plot(t,F,'-.','LineWidth',1,'color','k');
ylabel('力(N)');
%设置刻度
per = 0.8;      %定比例
cou = F_max;        %选参数
mag = floor(log10(cou / per));        %求数量级
fir = floor((cou / per) / (10^mag));      %取第一位
sec = floor((cou / per - fir*10^mag) / (10^(mag - 1)));       %取第二位
ran = fir*10^mag + sec*10^(mag - 1);        %定量程
cal = fir*10^(mag - 1);     %定刻度
sta = -0.5*fir*10^(mag - 1);        %定起点
axis([pri*t_max,prx*t_max,sta,ran]);
set(gca,'YTick',0:cal:ran);
%设置X轴和标题
xlabel('时间(s)');
title('压强与推力');
legend('压强', '推力');

%结束
