%180916
%Dashaway
%燃烧参数计算

%dp计算公式有误
%部分参数无合理来源

clear;
close all;
%初始值设定
%常量
%辅助常量
long = 5000;        %数组长度
n = 1:1:long;       %绘图横坐标
t = 0;      %时间
dt = 0.01;      %步进值
%已知常量
p0 = 1.02e5;    %初始压强
D = 0.132;       %外径
d = 0.06;       %初始内径
ep = (D - d)/2;     %总肉厚
Lp = 0.22;      %装药长度
gamma = 1.2;        %比热比
rho_p = 1730;       %密度
np = 0.302;      %压强指数
c = 1600;       %特征速度
alpha = 1.41e-3;
phi = 1; 
%计算得常量
Gamma = ( (2/(gamma + 1) )^( (gamma + 1)/(2*(gamma - 1) ) ) )*sqrt(gamma);      %比热比函数
%未确定常量
phi_m = 1;
At = 0.002;

%变量
%计算用变量（已赋初值）
e = 0*ones(1,long);     %已烧去肉厚
p = p0*ones(1,long);        %实际压强
Ab = pi*d*ones(1,long);     %燃烧面积
Ap = (pi/4)*(d^2)*ones(1,long);     %通气面积
Vg = Ap*Lp;     %自由体积
%循环变量
i = 1;


%数值计算
%循环前计算
a = rho_p*1.41e-3*Gamma^2*c^2*phi;
b = 1*Gamma^2*c*At;

%循环计算
while e <= ep
    dp = 5e4;
    %dp = -1e4;
    %dp = a*Ab(i)*p(i)^np/Vg(i)-b*p(i)/Vg(i);
    i = i + 1;
    p(i) = p(i - 1) + Runge_Kutta_4th(dp,dt);
    e(i) = e(i - 1) + alpha*p(i)^np*dt;
    %p(i) = p(i - 1) + dp;
    %e(i) = e(i - 1) + 1e-3;
    %p(i) = p(i - 1) + Runge_Kutta_4th(dp,dt);
    %e(i) = e(i - 1) + alpha*p(i)^np*dt;
    
    %参数修正
    Ab(i) = pi*(d+2*e(i) )*Lp;
    Ap(i) = (pi/4)*(d+2*e(i) )^2;
    Vg(i) = Ap(i)*Lp;
    if i >= long
        break;
    end
end

%循环后处理
i = i + 1;
n(i:1:long) = [];
e(i:1:long) = [];
p(i:1:long) = [];
Ab(i:1:long) = [];
Ap(i:1:long) = [];
Vg(i:1:long) = [];


%数据输出

figure;
hold on;
subplot(2,1,1);
plot(n,p);
grid on,axis tight ;
title('燃烧室压力')

subplot(2,1,2);
plot(n,e);
grid on,axis tight ;
title('烧去肉厚')

figure;
hold on;
subplot(2,1,1);
plot(n,Ab);
grid on,axis tight ;
title('燃烧面积')

subplot(2,1,2);
plot(n,Ap);
grid on,axis tight ;
title('通气面积')



