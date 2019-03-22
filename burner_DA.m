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
long = 100000;        %数组长度
n = 1:1:long;       %绘图横坐标
dt = 2e-5;      %步进值
%燃烧室参数
Dr = 0.132;       %燃烧室外径(m)
% At = 1.884785e-3;     %喷管喉部面积(m^2)
At = 5e-4;     %喷管喉部面积(m^2)
%装药参数
%星孔装药参数
Ds = 0.132;       %外径(m)
ep_s = 0.04;     %肉厚(m)
r_s = 0.003;        %星尖圆弧半径(m)
r1_s = 0.003;        %星根过渡圆弧半径(m)
l_s = 0.023;        %特征长度(m)
n_s = 7;        %星角数
epsilon_s = 1;      %角分数
theta_s = 0.620465;      %星根半角(rad)
%圆柱装药参数
Dc = 0.132;     %外径(m)
rc = 0.05;      %内径(m)
ep_c = (Dc - rc) / 2;       %肉厚(m)
%分段参数
Lp = 0.22;      %药柱长度(m)
pers = 0.5;     %星孔段占比
Lp_s = Lp*pers;     %星孔段长度
Lp_c = Lp - Lp_s;       %圆柱段长度
%已知常量
p0 = 1.02e5;    %初始压强(Pa)
gamma = 1.2;        %比热比
rho_p = 1730;       %密度(kg/m^3)
n_p = 0.302;      %压强指数
c = 1600;       %特征速度(m/s)
rb_0 = 5e-3;      %?初始燃速(m/s)
alpha_r = 1.7e-4;        %?燃速系数
phi_alpha = 1;       %?侵蚀函数
phi_m = 1;      %?

%计算得常量
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数
beta_s = pi / n_s;        %等分角(rad)

%计算中间值
%星孔装药
%燃烧周长项中参数
s_a = (sin(epsilon_s*beta_s)) / sin(theta_s) + (1 - epsilon_s)*beta_s;
s_b = (pi / 2 + beta_s - theta_s - cot(theta_s));
s_c = (1 - epsilon_s)*beta_s;
s_d =  (l_s*sin(epsilon_s*beta_s));
%通气面积项中参数
Ap_a = ((l_s^2) / 2)*( (1 - epsilon_s)*beta_s + (sin(epsilon_s*beta_s)) ...
    *( (cos(epsilon_s*beta_s)) - (sin(epsilon_s*beta_s))*cot(theta_s) ) );
Ap_b = s_a;
Ap_c = s_b;
Ap_d = cot(theta_s) + theta_s - pi / 2;
Ap_e = s_c;
Ap_f = s_d;
Ap_g = epsilon_s*beta_s;
%压强项中参数
p_a = rho_p*alpha_r*phi_alpha*Gamma^2*c^2;
p_b = phi_m*Gamma^2*c*At;

%初始参数
%星孔装药
%周长参数
s_s_a0 = 2*n_s*(l_s*s_a + (r_s + r1_s)*s_b  - beta_s*r1_s);      
%通气面积参数
Ap_si0 = Ap_a + l_s*r_s*Ap_b + ((r_s^2) / 2)*Ap_c + ((r1_s^2) / 2)*Ap_d;
Ap_s_a0 = 2*n_s*Ap_si0;
%圆柱装药
s_c_0 = 2*pi*rc;
Ap_c_0 = pi*(rc^2);


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
s_s = s_s_a0*ones(1,long);     %星孔燃烧面实际边长(m)
s_c = s_c_0*ones(1,long);     %圆柱燃烧面实际边长(m)
m_b = rho_p*s_a0*Lp*rb_0*ones(1,long);     %燃气生成率(kg/s)
m_p = (phi_m*p0*At / c)*ones(1,long);     %质量流率(kg/s)
F = 2000*(p0*At / c)*ones(1,long);     %推力(N)
Ab = s_a0*Lp*ones(1,long);     %燃烧面积(m^2)
Ap_s = Ap_s_a0*ones(1,long);     %星孔通气面积(m^2)
Ap_c = Ap_c_0*ones(1,long);     %圆柱通气面积(m^2)
Vg = Ap*Lp;     %自由体积(m^3)
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
%循环变量
i = 1;
j = 1;
sw = 0*ones(1,long);        %阶段选择
swc = 0*ones(1,100);        %阶段变化点

%计算中间值

%各阶段燃烧周长
%     si1 = l_s*s_a + (r_s + r1_s)*s_b - beta_s*(r1_s - e);
%     si2 = l_s*s_a + (r_s + e)*s_b;
%     si3 = l_s*s_c + (r_s + e)*( beta_s + asin(s_d / (r_s + e)) );
%     s = 2*n_s*si;
%各阶段通气面积
%     Api1 = Ap_a + l_s*(r_s + e)*Ap_b + (((r_s + e)^2) / 2)*Ap_c ...
%        + (((r1_s - e)^2) / 2)*Ap_d;
%     Api2 = Ap_a + l_s*(r_s + e)*Ap_b + (((r_s + e)^2) / 2)*Ap_c;
%     Api3 = ( ((l_s + r_s + e)^2)*Ap_e ...
%        + ((r_s + e)^2)*(Ap_g + asin(Ap_f / (r_s + e))) ...
%        + Ap_f*(sqrt((r_s + e)^2 - Ap_f^2)  + l_s*cos(Ap_g)) ) / 2;
%     Ap = 2*n_s*Api;

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
        swc(j) = i;
        j = j + 1;
    end
    %分阶段计算燃烧周长和通气面积
    switch sw(i)
        case{1}     %星孔第一阶段，圆柱第一阶段
            s(i) = 2*n_s*(l_s*s_a + (r_s + r1_s)*s_b - beta_s*(r1_s - e(i)));
            Api1 = Ap_a + l_s*(r_s + e(i))*Ap_b + (((r_s + e(i))^2) / 2)*Ap_c ...
                + (((r1_s - e(i))^2) / 2)*Ap_d;
            Ap(i) = 2*n_s*Api1;
        case{3}     %星孔第二阶段，圆柱第一阶段
            s(i) = 2*n_s*(l_s*s_a + (r_s + e(i))*s_b);
            Api2 = Ap_a + l_s*(r_s + e(i))*Ap_b + (((r_s + e(i))^2) / 2)*Ap_c;
            Ap(i) = 2*n_s*Api2;
        case{7}     %星孔第三阶段，圆柱第一阶段
            s(i) = 2*n_s*(l_s*s_c + (r_s + e(i))*( beta_s + asin(s_d / (r_s + e(i))) ) );
            Api3 = ( ((l_s + r_s + e(i))^2)*Ap_e ...
                + ((r_s + e(i))^2)*(Ap_g + asin(Ap_f / (r_s + e(i)))) ...
                + Ap_f*(sqrt((r_s + e(i))^2 - Ap_f^2)  + l_s*cos(Ap_g)) ) / 2;
            Ap(i) = 2*n_s*Api3;
        case{15}     %星孔结束，圆柱第一阶段
            s(i) = 2*n_s*(l_s*s_c + (r_s + e(i))*( beta_s + asin(s_d / (r_s + e(i))) ) );
            Api3 = ( ((l_s + r_s + e(i))^2)*Ap_e ...
                + ((r_s + e(i))^2)*(Ap_g + asin(Ap_f / (r_s + e(i)))) ...
                + Ap_f*(sqrt((r_s + e(i))^2 - Ap_f^2)  + l_s*cos(Ap_g)) ) / 2;
            Ap(i) = 2*n_s*Api3;
        case{17}     %星孔第一阶段，圆柱结束
            s(i) = 2*n_s*(l_s*s_a + (r_s + r1_s)*s_b - beta_s*(r1_s - e(i)));
            Api1 = Ap_a + l_s*(r_s + e(i))*Ap_b + (((r_s + e(i))^2) / 2)*Ap_c ...
                + (((r1_s - e(i))^2) / 2)*Ap_d;
            Ap(i) = 2*n_s*Api1;
        case{19}     %星孔第二阶段，圆柱结束
            s(i) = 2*n_s*(l_s*s_a + (r_s + e(i))*s_b);
            Api2 = Ap_a + l_s*(r_s + e(i))*Ap_b + (((r_s + e(i))^2) / 2)*Ap_c;
            Ap(i) = 2*n_s*Api2;
        case{23}     %星孔第三阶段，圆柱结束
            s(i) = 2*n_s*(l_s*s_c + (r_s + e(i))*( beta_s + asin(s_d / (r_s + e(i))) ) );
            Api3 = ( ((l_s + r_s + e(i))^2)*Ap_e ...
                + ((r_s + e(i))^2)*(Ap_g + asin(Ap_f / (r_s + e(i)))) ...
                + Ap_f*(sqrt((r_s + e(i))^2 - Ap_f^2)  + l_s*cos(Ap_g)) ) / 2;
            Ap(i) = 2*n_s*Api3;
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
