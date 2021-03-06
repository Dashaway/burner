%181016
%Dashaway
%燃烧参数计算


%190329
%双Y轴按比例显示

%190414
%参数外部设置
%多重图像



%初始值设定


%常量

%计算辅助常量
long = 150000;        %数组长度
n = 1:1:long;       %绘图横坐标
dt = 5e-5;      %步进值
%燃烧室参数
Dr = 0.300;       %燃烧室外径(m)
At = 5.026548e-3;     %喷管喉部面积(m^2)

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


%计算得常量
Gamma = ( (2 / (gamma + 1))^( (gamma + 1) / (2*(gamma - 1)) ) ) ...
    *sqrt(gamma);      %比热比函数
trd = (2*pi) / (n_s*(m_s + 1));     %药圆心角
trr = (2*pi*m_s) / (n_s*(m_s + 1));     %弧圆心角

s_a0 = 2*n_s*trr*l + 2*pi*(R + n_s*r);      %燃烧面初始周长
Ap_a0 = 2*n_s*trr*l*r + n_s*pi*r^2 + pi*R^2;     %初始通气面积

Vp0 = (pi*(Dr^2) / 4 - Ap_a0)*Lp;       %药柱体积(m^3)
mp = Vp0*rho_p;     %药柱质量(kg)

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
df = zeros(1,long);
ddf = zeros(1,long);
st = zeros(2,100);
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
sw = zeros(1,long);        %阶段选择
swc = zeros(2,100);        %阶段变化点

%计算中间值
%     as1 =  asin( (l*sin(trd/2)) / (r + e) );
%     ac1 = acos( (4*l^2 + D^2 - 4*(r + e)^2) / (4*l*D) ); 
%     ac2 = acos( (4*(r + e)^2 + 4*l^2 - D^2) / (8*l*(r + e)) ); 
%     ac3 = acos( (l^2 + (r + e)^2 - (R + e)^2) / (2*l*(r + e)) );
%     ac4 = acos( (l^2 + (R + e)^2 - (r + e)^2) / (2*l*(R + e)) );
%各阶段燃烧周长
%     s1 =  2*n_s*trr*l + 2*pi*(n_s*(r + e) + (R + e));
%     s2a = 2*n_s*ac2*(r + e) + 2*n_s*(trr/2)*(l - r - e) + 2*pi*(R + e);
%     s2b =  2*n_s*trr*l + 2*pi*(R + e) + 4*n_s*as1*(r + e);
%     s3 = 2*n_s*(r + e)*(ac2 + 2*as1 - pi) + 2*pi*(R + e) + n_s*trr*(l - r - e);
%     s4 = 2*n_s*(r + e)*(as1 - trd/2) + 2*pi*(R + e) + n_s*trr*(l - r - e);
%     s5 = 2*n_s*(r + e)*(as1 - trd/2 - ac3) + 2*n_s*(R + e)*(trd/2 - ac4);  
%各阶段通气面积
%     Ap1 = 2*n_s*trr*l*(r + e) + n_s*pi*(r + e)^2 + pi*(R + e)^2;
%     Ap2a = pi*(R + e)^2 + n_s*ac2*(r + e)^2 ...  
%         + n_s*(((D/2)^2)*ac1 - n_s*l*(D/2)*sin(ac1) ...
%         + (n_s/2)*trr*( (D/2)^2 - (l - r - e)^2 );
%     Ap2b = 2*n_s*trr*l*(r + e) + pi*(R + e)^2 ...
%         + n_s*((r + e)^2)*(2*as1 + sin(2*as1));
%     Ap3 = pi*(R + e)^2 + n_s*ac2*(r + e)^2 ...  
%         + n_s*((D/2)^2)*ac1 - n_s*l*(D/2)*sin(ac1) ...
%         + (n_s/2)*trr*( (D/2)^2 - (l - r - e)^2 ) ....
%         - n_s*((r + e)^2)*(pi - 2*as1 - sin(2*as1));
%     Ap4 = pi*((D/2)^2) + pi*(R + e)^2 - (trr/2)*n_s*((l - r - e)^2);
%         - n_s*l*(r + e)*sin(as1 - trd/2) + n_s*((r + e)^2)*(as1 - trd/2);
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
            s(i) =  2*n_s*trr*l + 2*pi*(n_s*(r + e(i)) + (R + e(i)));
            Ap(i) = 2*n_s*trr*l*(r + e(i)) + n_s*pi*(r + e(i))^2 + pi*(R + e(i))^2;
        case{3}     %第二A阶段
            ac1 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / (4*l*D) );
            ac2 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) );
            s(i) = 2*n_s*ac2*(r + e(i)) + 2*n_s*(trr/2)*(l - r - e(i)) + 2*pi*(R + e(i));
            Ap(i) = pi*(R + e(i))^2 + n_s*ac2*(r + e(i))^2 ...
                + n_s*((D/2)^2)*ac1 - n_s*l*(D/2)*sin(ac1) ...
                + (n_s/2)*trr*( (D/2)^2 - (l - r - e(i))^2 );
        case{5}     %第二B阶段
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            s(i) =  2*n_s*trr*l + 2*pi*(R + e(i)) + 4*n_s*as1*(r + e(i));
            Ap(i) = 2*n_s*trr*l*(r + e(i))  + pi*(R + e(i))^2 ...
                + n_s*((r + e(i))^2)*(2*as1 + sin(2*as1));
        case{7}     %第三阶段
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            ac1 = acos( (4*l^2 + D^2 - 4*(r + e(i))^2) / (4*l*D) );
            ac2 = acos( (4*(r + e(i))^2 + 4*l^2 - D^2) / (8*l*(r + e(i))) );
            s(i) = 2*n_s*(r + e(i))*(ac2 + 2*as1 - pi) + 2*pi*(R + e(i)) ...
                + n_s*trr*(l - r - e(i));
            Ap(i) = pi*(R + e(i))^2 + n_s*ac2*(r + e(i))^2 ...
                + n_s*((D/2)^2)*ac1 - n_s*l*(D/2)*sin(ac1) ...
                + (n_s/2)*trr*( (D/2)^2 - (l - r - e(i))^2 ) ....
                - n_s*((r + e(i))^2)*(pi - 2*as1 - sin(2*as1));
        case{15}     %第四阶段
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            s(i) = 2*n_s*(r + e(i))*(as1 - trd/2) + 2*pi*(R + e(i)) + n_s*trr*(l - r - e(i));
            Ap(i) = pi*((D/2)^2) + pi*((R + e(i))^2) - (trr/2)*n_s*((l - r - e(i))^2) ...
                - n_s*l*(r + e(i))*sin(as1 - trd/2) + n_s*((r + e(i))^2)*(as1 - trd/2);
        case{31}     %第五阶段
            as1 =  asin( (l*sin(trd/2)) / (r + e(i)) );
            ac3 = acos( (l^2 + (r + e(i))^2 - (R + e(i))^2) / (2*l*(r + e(i))) );
            ac4 = acos( (l^2 + (R + e(i))^2 - (r + e(i))^2) / (2*l*(R + e(i))) );
            s(i) = 2*n_s*(r + e(i))*(as1 - trd/2 - ac3) + 2*n_s*(R + e(i))*(trd/2 - ac4);
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
    if(F(i) < 0 )
        F(i) = 0;
    end
    if (i  <= 1)
        df(i) = 0;
    else
        df(i) = F(i) - F(i - 1);
    end
    if (i  <= 1)
        ddf(i) = 0;
    else
        ddf(i) = df(i) - df(i - 1);
    end
    
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
df(i:1:long) = [];
ddf(i:1:long) = [];
t = dt*n;
t_max = dt*(i - 1);     %燃烧总时间
prx = 1.05;
pri = -0.05;



str = [' t_max = ',num2str(t_max)];
disp(str);

%数据输出






%结束
