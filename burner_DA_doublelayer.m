%181016
%Dashaway
%燃烧参数计算


%190321
%分段燃烧

%190326
%内外层不同装药


clear;
close all;
%初始值设定


%常量

%计算辅助常量
long = 150000;        %数组长度
n = 1:1:long;       %绘图横坐标
dt = 1e-4;      %步进值
%燃烧室参数
Dr = 0.300;       %燃烧室外径(m)
At = 5.026548e-3;     %喷管喉部面积(m^2)
%装药参数
%内外圆柱装药参数
Dc = 0.300;     %外径(m)
dm = 0.100;     %中径(m)
dc = 0.050;      %内径(m)
ep_m = (dm - dc) / 2;       %内肉厚(m)
ep_c = (Dc - dc) / 2;       %总肉厚(m)
Lp = 0.500;      %药柱长度(m)
%内推进剂参数
n_p1 = 0.4;      %压强指数
rho_p1 = 1700;       %密度(kg/m^3)
c1 = 1584;       %特征速度(m/s)
rb1_0 = 8.77e-5;      %?初始燃速(m/s)
alpha_r1 = 3e-5;        %?燃速系数
%外推进剂参数
n_p2 = 0.302;      %压强指数
rho_p2 = 1730;       %密度(kg/m^3)
c2 = 1600;       %特征速度(m/s)
rb2_0 = 5e-3;      %?初始燃速(m/s)
alpha_r2 = 1.71e-4;        %?燃速系数
%已知常量
p0 = 1.02e5;    %初始压强(Pa)
gamma = 1.2;        %比热比
phi_alpha = 1;       %?侵蚀函数
phi_m = 1;      %?

%计算得常量
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数

%计算中间值

%压强项中参数
p_a1 = rho_p1*alpha_r1*phi_alpha*Gamma^2*c1^2;
p_b1 = phi_m*Gamma^2*c1*At;
p_a2 = rho_p2*alpha_r2*phi_alpha*Gamma^2*c2^2;
p_b2 = phi_m*Gamma^2*c2*At;

%初始参数
%圆柱装药
%周长参数
s_0 = pi*dc;
%通气面积参数
Ap_0 = pi*(dc^2) / 4;


Vp10 = (pi*(dm^2) / 4 - pi*(dc^2) / 4)*Lp;       %内药柱体积(m^3)
Vp20 = (pi*(Dr^2) / 4 - pi*(dm^2) / 4)*Lp;       %外药柱体积(m^3)
Vp0 = Vp10 + Vp20;      %总药柱体积(m^3)
mp1 = Vp10*rho_p1;     %内药柱质量(kg)
mp2 = Vp20*rho_p2;     %外药柱质量(kg)
mp = mp1 + mp2;     %总药柱质量(kg)

%压强项
p_a = p_a1;
p_b = p_b1;
n_p = n_p1;
alpha_r = alpha_r1;

%约束条件
%参数约束条件

%分阶段判别条件


%各阶段对应条件
%     星孔第一阶段，圆柱第一阶段
%     e1 = ((e <= ep_m) & (e <= ep_c));
%     sw = 1;
%     星孔第二阶段，圆柱第一阶段
%     e2 = ((e > ep_m) & (e <= ep_c));
%     sw = 3;
%程序结束条件
ep = ep_c;
%     e >= ep;


%变量
%计算用变量数组（已赋初值）
p = p0*ones(1,long);        %实际压强(Pa)
rb = rb1_0*ones(1,long);     %燃速(m/s)
e = 0*ones(1,long);     %已烧去肉厚(m)
s = s_0*ones(1,long);     %燃烧面实际周长(m)
m_b = rho_p1*s_0*Lp*rb1_0*ones(1,long);     %燃气生成率(kg/s)
m_p = (phi_m*p0*At / c1)*ones(1,long);     %质量流率(kg/s)
F = 2000*(p0*At / c1)*ones(1,long);     %推力(N)
Ab = s_0*Lp*ones(1,long);     %燃烧面积(m^2)
Ap = Ap_0*ones(1,long);     %通气面积(m^2)
Vg = Ap*Lp;     %自由体积(m^3)
%最大值变量
p_max = p0;        %最大压强(Pa)
e_max = 0;     %最大已烧去肉厚(m)
s_max = s_0;     %燃烧面最大周长(m)
m_b_max = rho_p1*s_0*Lp*rb1_0;     %最大燃气生成率(kg/s)
m_p_max = (phi_m*p0*At / c1);     %最大质量流率(kg/s)
F_max = 2000*(p0*At / c1);     %最大推力(N)
Ab_max = s_0*Lp;     %最大燃烧面积(m^2)
Ap_max = Ap_0;     %最大通气面积(m^2)
Vg_max = Ap_0*Lp;     %最大自由体积(m^3)
%最小值变量
p_min = p0;        %最小压强(Pa)
e_min = 0;     %最小已烧去肉厚(m)
s_min = s_0;     %燃烧面最小周长(m)
m_b_min = rho_p1*s_0*Lp*rb1_0;     %最小燃气生成率(kg/s)
m_p_min = (phi_m*p0*At / c1);     %最小质量流率(kg/s)
F_min = 2000*(p0*At / c1);     %最小推力(N)
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
    s(i) = pi*(dc + 2*e(i));
    Ap(i) = pi*(dc + 2*e(i))^2 / 4;
    %判断此时处于哪个阶段
    sw(i) = 1;
    if(e(i) <= ep_m)
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
        case{1}     %内装药阶段
            p_a = p_a1;
            p_b = p_b1;
            n_p = n_p1;
            alpha_r = alpha_r1;
            m_b(i) = rho_p1*Ab(i)*rb(i);
            m_p(i) = phi_m*p(i)*At / c1;
        case{3}     %外装药阶段
            p_a = p_a2;
            p_b = p_b2;
            n_p = n_p2;
            alpha_r = alpha_r2;
            m_b(i) = rho_p2*Ab(i)*rb(i);
            m_p(i) = phi_m*p(i)*At / c2;
        otherwise     %其它阶段
            m_b(i) = m_b(i - 1);
            m_p(i) = m_p(i - 1);
    end
    %计算其他参数
    Ab(i) = s(i)*Lp;
    Vg(i) = Ap(i)*Lp;
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


%双y轴2014版
figure;
hold on;
%AX(1)和 AX(2)分别是左右 axes 的句柄，可以用 set()函数处理
[AX,H1,H2] = plotyy(t,p,t,F); 
set(AX,'Xlim',[pri*t_max,prx*t_max]);
xlabel('时间(s)');
title('压强与推力');
set(AX(:),'Ycolor','k');

set(get(AX(1),'Ylabel'),'string','压强(Pa)','color','k','linewidth',1.0); 
set(get(AX(2),'Ylabel'),'string','力(N)','color','k','linewidth',1.0); 

per = 0.6;      %定比例
cou = p_max;        %选参数
mag = floor(log10(cou / per));        %求数量级
fir = floor((cou / per) / (10^mag));      %取第一位
sec = floor((cou / per - fir*10^mag) / (10^(mag - 1)));       %取第二位
ran = fir*10^mag + sec*10^(mag - 1);        %定量程
cal = fir*10^(mag - 1);     %定刻度
sta = -0.5*fir*10^(mag - 1);        %定起点
set(AX(1),'Ylim',[sta,ran],'yTick',0:cal:ran);

per = 0.8;      %定比例
cou = F_max;        %选参数
mag = floor(log10(cou / per));        %求数量级
fir = floor((cou / per) / (10^mag));      %取第一位
sec = floor((cou / per - fir*10^mag) / (10^(mag - 1)));       %取第二位
ran = fir*10^mag + sec*10^(mag - 1);        %定量程
cal = fir*10^(mag - 1);     %定刻度
sta = -0.5*fir*10^(mag - 1);        %定起点
set(AX(2),'Ylim',[sta,ran],'yTick',0:cal:ran);

box off;
set(H1,'LineStyle','-','color','k','linewidth',1.0);
set(H2,'LineStyle','-.','color','k','linewidth',1.0);
legend('压强','推力');

%结束
