clear; clc;
%% 参数设置
lambda = 1064e-9;  % 激光波长（单位：米）
R_75 = 75e-9;         % 月尘半径（单位：米）
PCE = 0.51;        % 激光电-光功率转换效率
EHCE = 0.264;      % 能量收集转换效率
dr = 2.1;          % 接收器孔径直径
phit=2e-6;
phir=2e-6;
nm_real = 1.733;        % 月尘的折射率
nm_imag = 0.01;
nm = nm_real + 1i * nm_imag;
%% 计算尺寸参数
k = 2 * pi / lambda;  % 波数
r = linspace(0.0, 200e-9, 1000);  % 粒子半径范围（从1nm到200nm）
x = k * r;  % 尺寸参数
mx = nm * x;  % 尺寸参数乘折射率
%% 衰减系数与密度和粒子半径的关系
n_max = 50;  % 设置最大阶数
Cext = zeros(1, length(r));
for n = 1:n_max
    % 计算贝塞尔函数 j_n(x) 和 j_n(mx)
    Jx = besselj(n + 0.5, x);  % 贝塞尔函数 j_n(x)
    Jmx = besselj(n + 0.5, mx); % 贝塞尔函数 j_n(mx)
    phix = sqrt(pi * x / 2) .* Jx;
    phimx = sqrt(pi * mx / 2) .* Jmx;
    % 计算导数
    phix_ = gradient(phix, x);  % 计算phix的导数
    phimx_ = gradient(phimx, mx);
    % 计算Hankel函数
    Hx = besselh(n + 0.5, 1, x);  % Hankel函数
    epsilonx = sqrt(pi * x / 2) .* Hx;
    epsilonx_ = gradient(epsilonx, x);  
    % 计算 an 和 bn
    an = (nm * phimx .* phix_ - phix .* phimx_) ./ (nm * phimx .* epsilonx_ - epsilonx .* phimx_);
    bn = (phimx .* phix_ - nm * phix .* phimx_) ./ (phimx .* epsilonx_ - nm * epsilonx .* phimx_);   
    % 计算消光截面并累加
    Cext = Cext + (2 * n + 1) * real(an + bn);
end
% 消光截面
Cext_total = (2 * pi / k^2) * Cext;
%计算衰减系数 alpha
N=linspace(0,20000,1000)*1e6;  %设置月尘密度计算范围
alpha = zeros(length(N), length(r));
for i = 1:length(N)
    alpha(i, :) = N(i) * Cext_total;  % 对应每个粒子密度计算衰减系数
end
%%绘图
[R, N_grid] = meshgrid(r, N);  
figure(1)
surf(R*1e9, N_grid*1e-6, alpha);
xlabel('粒子半径 (nm)');
ylabel('粒子密度(cm^3)');
zlabel('衰减系数(m^{-1})');
title('衰减系数与粒子密度和粒子半径的关系');
colorbar;  %颜色条
shading interp;%去除网格线  
grid on;

%% 光照区Ldust与传输距离和月球表面以上高度的函数
%% 计算75 nm半径的消光截面
x_75=k*R_75;  % 尺寸参数
mx_75=nm*x_75;  % 尺寸参数乘以折射率
n_max=50;  % 最大阶数
Cext_75nm=0;  % 初始化消光截面值
for n = 1:n_max
    % 计算贝塞尔函数 j_n(x) 和 j_n(mx)
    Jx_75=besselj(n+0.5,x_75); % 贝塞尔函数 j_n(x)
    Jmx_75=besselj(n+0.5,mx_75); % 贝塞尔函数 j_n(mx)
    phix_75=sqrt(pi*x_75/2).*Jx_75;
    phimx_75=sqrt(pi*mx_75/2).*Jmx_75;    
    % 计算导数
    phix__75=(sqrt(pi*(x_75+1e-12)/2).*besselj(n+0.5,x_75+1e-12)-phix_75)/1e-12;
    phimx__75=(sqrt(pi*(mx_75+1e-12)/2).*besselj(n+0.5,mx_75+1e-12)-phimx_75)/1e-12;   
    % 计算Hankel函数
    Hx_75=besselh(n+0.5,1,x_75);  
    epsilonx_75=sqrt(pi*x_75/2).*Hx_75;
    epsilonx__75=(sqrt(pi*(x_75+1e-12)/2).*besselh(n+0.5,1,x_75+1e-12)-epsilonx_75)/1e-12;    
    % 计算an 和 bn
    an_75=(nm*phimx_75.*phix__75-phix_75.*phimx__75)./(nm*phimx_75.*epsilonx__75-epsilonx_75.*phimx__75);
    bn_75=(phimx_75.*phix__75-nm*phix_75.*phimx__75)./(phimx_75.*epsilonx__75-nm*epsilonx_75.*phimx__75);
    % 累加消光截面
    Cext_75nm=Cext_75nm+(2*n+1)*real(an_75+bn_75);
end
% 最终消光截面（单位：m^2）
Cext_75nm=(2*pi/k^2)*Cext_75nm;

%% 计算对数分布75nm的消光截面
mu = log(75e-9);  % 均值
sigma = 0.2;  % 标准差
R = linspace(10e-9, 300e-9, 10000);  % 粒径范围从10nm到300nm
Cext = zeros(size(R));
% 计算每个粒径下的消光截面
for i = 1:length(R)
    r = R(i);
    x = k * r;  % 尺寸参数
    mx = nm * x;  % 尺寸参数乘折射率
    Cext_r = 0;  % 初始消光截面值
    for n = 1:n_max
        % 计算贝塞尔函数 j_n(x) 和 j_n(mx)
        Jx = besselj(n+0.5, x);  
        Jmx = besselj(n+0.5, mx);  
        phix = sqrt(pi * x / 2) .* Jx;
        phimx = sqrt(pi * mx / 2) .* Jmx;    
        % 计算导数
        phix_ = (sqrt(pi * (x + 1e-12) / 2) .* besselj(n+0.5, x + 1e-12) - phix) / 1e-12;
        phimx_ = (sqrt(pi * (mx + 1e-12) / 2) .* besselj(n+0.5, mx + 1e-12) - phimx) / 1e-12;   
        % 计算Hankel函数
        Hx = besselh(n+0.5, 1, x);  
        epsilonx = sqrt(pi * x / 2) .* Hx;
        epsilonx_ = (sqrt(pi * (x + 1e-12) / 2) .* besselh(n+0.5, 1, x + 1e-12) - epsilonx) / 1e-12;    
        % 计算an 和 bn
        an = (nm * phimx .* phix_ - phix .* phimx_) ./ (nm * phimx .* epsilonx_ - epsilonx .* phimx_);
        bn = (phimx .* phix_ - nm * phix .* phimx_) ./ (phimx .* epsilonx_ - nm * epsilonx .* phimx_);
        % 累加消光截面
        Cext_r = Cext_r + (2 * n + 1) * real(an + bn);
    end
    % 最终消光截面（单位：m^2）
    Cext(i) = (2 * pi / k^2) * Cext_r;
end
% 定义对数高斯分布的概率密度函数
lognormal_pdf = @(r) (1 ./ (r * sigma * sqrt(2 * pi))) .* exp(-((log(r) - mu).^2) / (2 * sigma^2));
% 计算平均消光截面
Cext_avg = trapz(R, Cext .* lognormal_pdf(R)) / trapz(R, lognormal_pdf(R));

%% 光照区的Ldust与传输距离和月球表面以上高度的函数
d_max = 100e3;   % 最大传输距离 (单位: m)
h_max = 2;    % 最大高度 (单位: m)
d = linspace(50, d_max, 1000);  % 传输距离 (从 50m 到 100 km)
h = linspace(0, h_max, 200);  % 高度 (从 0 到 200 cm)
[D, H] = meshgrid(d, h);      % 创建传输距离和高度的网格
%%衰减系数计算
alpha_h = Cext_75nm*(-9.5519e9*H.^5+ 5.8237e10*H.^4-1.382e11*H.^3+1.59987e11*H.^2-9.307052e10*H+2.5539e10);  % 衰减系数随高度的变化
Ldust = exp(-alpha_h .* D);  % Ldust 随高度和传输距离变化
alpha_havg=Cext_avg*(-9.5519e9*H.^5+ 5.8237e10*H.^4-1.382e11*H.^3+1.59987e11*H.^2-9.307052e10*H+2.5539e10);
Ldust_avg = exp(-alpha_havg .* D);
%%绘制三维图
figure(2);
surf(H, D, Ldust_avg);  % 绘制三维曲面图
xlabel('月面以上高度 H (m)');
ylabel('传输距离 R (m)');
zlabel('L_{dust}');
title('光照区L_{dust} 与月面高度H和传输距离R的关系');
colorbar;  % 添加颜色条
shading interp; %去除网格线  
grid on;
%% 阴影区的Ldust与传输距离和月球表面以上高度的函数
%%衰减系数计算
alpha_h_shadow = Cext_75nm*((-9.5519e9*H.^5+ 5.8237e10*H.^4-1.382e11*H.^3+1.59987e11*H.^2-9.307052e10*H+2.5539e10)*1e-4);  % 衰减系数随高度的变化
Ldust_shadow = exp(-alpha_h_shadow .* D);  % Ldust 随高度和传输距离变化
alpha_h_shadowavg = Cext_avg*((-9.5519e9*H.^5+ 5.8237e10*H.^4-1.382e11*H.^3+1.59987e11*H.^2-9.307052e10*H+2.5539e10)*1e-4);  % 衰减系数随高度的变化
Ldust_shadowavg = exp(-alpha_h_shadowavg .* D);
%%绘制三维图
figure(3);
surf(H, D, Ldust_shadowavg);  % 绘制三维曲面图
zlim([0,1]);
xlabel('月面以上高度 H (m)');
ylabel('传输距离 R (m)');
zlabel('L_{dust}');
title('阴影区L_{dust} 与 月面高度 H 和传输距离 R 的关系');
colorbar;  % 添加颜色条
shading interp; %去除网格线  
grid on;
%% 光照区 不同高度的平均收获功率随传输功率和距离的关系
Hgiven=[0.2;0.5;1];
alpha_Hgiven=Cext_75nm * (-9.5519e9*Hgiven.^5+ 5.8237e10*Hgiven.^4-1.382e11*Hgiven.^3+1.59987e11*Hgiven.^2-9.307052e10*Hgiven+2.5539e10);
alpha_Hgivenavg=Cext_avg * (-9.5519e9*Hgiven.^5+ 5.8237e10*Hgiven.^4-1.382e11*Hgiven.^3+1.59987e11*Hgiven.^2-9.307052e10*Hgiven+2.5539e10);
%%参数设置
d_max = 1e5;        % 最大传输距离 (单位: m)
p_max = 1e4;        % 最大功率 (单位: w)
d = linspace(50, d_max, 1000);  % 传输距离 (从 50 到 100 km)
pt = linspace(0, p_max, 1000);  % 功率 (从 0 到 10kw)
[D, Pt] = meshgrid(d, pt);      
theta_=dr./D;%发散角
dt=1.22*lambda./theta_;  % 计算发射器孔径 (随距离变化)
dr = repmat(dr, size(dt));
Gt=(pi*dt./lambda).^2;
Gr=(pi*dr./lambda).^2;
Lt=exp(-Gt*phit^2);
Lr=exp(-Gt*phir^2);
KS=lambda./(4*pi*D);
KS2=KS.^2;
%20cm

Ldust_20cm=exp(-alpha_Hgiven(1).*D);
Ph_20cm=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_20cm;
%50cm
Ldust_50cm=exp(-alpha_Hgiven(2).*D);
Ph_50cm=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_50cm;
%100cm
Ldust_100cm=exp(-alpha_Hgiven(3).*D);
Ph_100cm=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_100cm;
%20cm
Ldust_20cmavg=exp(-alpha_Hgivenavg(1).*D);
Ph_20cmavg=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_20cmavg;
%50cm
Ldust_50cmavg=exp(-alpha_Hgivenavg(2).*D);
Ph_50cmavg=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_50cmavg;
%100cm
Ldust_100cmavg=exp(-alpha_Hgivenavg(3).*D);
Ph_100cmavg=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_100cmavg;
% 绘制三维图
figure(4);
surf(D, Pt, Ph_20cmavg, 'FaceColor', 'r','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
hold on;
surf(D, Pt, Ph_50cmavg, 'FaceColor', 'g','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
surf(D, Pt, Ph_100cmavg, 'FaceColor', 'b','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
grid on;
xlabel('传输距离 d (m)');
ylabel('传输功率 Pt (W)');
zlabel('平均收获功率 Ph (W)');
title('光照区不同高度的平均收获功率随传输功率和距离的关系');
legend('20 cm', '50 cm', '100 cm');
dr = 2.1;          % 接收器孔径直径
%% 阴影区 不同高度的平均收获功率随传输功率和距离的关系
alpha_Hgiven_shadow=Cext_75nm*((-9.5519e9*Hgiven.^5+ 5.8237e10*Hgiven.^4-1.382e11*Hgiven.^3+1.59987e11*Hgiven.^2-9.307052e10*Hgiven+2.5539e10)*1e-4);
alpha_Hgiven_shadowavg=Cext_avg*((-9.5519e9*Hgiven.^5+ 5.8237e10*Hgiven.^4-1.382e11*Hgiven.^3+1.59987e11*Hgiven.^2-9.307052e10*Hgiven+2.5539e10)*1e-4);
%20cm
Ldust_20cm_shadow=exp(-alpha_Hgiven_shadow(1)*D);
Ph_20cm_shadow=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_20cm_shadow;
%50cm
Ldust_50cm_shadow=exp(-alpha_Hgiven_shadow(2)*D);
Ph_50cm_shadow=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_50cm_shadow;
%100cm
Ldust_100cm_shadow=exp(-alpha_Hgiven_shadow(3)*D);
Ph_100cm_shadow=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_100cm_shadow;
%20cm
Ldust_20cm_shadowavg=exp(-alpha_Hgiven_shadowavg(1)*D);
Ph_20cm_shadowavg=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_20cm_shadowavg;
%50cm
Ldust_50cm_shadowavg=exp(-alpha_Hgiven_shadowavg(2)*D);
Ph_50cm_shadowavg=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_50cm_shadowavg;
%100cm
Ldust_100cm_shadowavg=exp(-alpha_Hgiven_shadowavg(3)*D);
Ph_100cm_shadowavg=Pt.*KS2.*(PCE*EHCE.*Lt.*Gt.*Lr.*Gr).*Ldust_100cm_shadowavg;
% 绘制三维图
figure(5);
surf(D, Pt, Ph_20cm_shadowavg, 'FaceColor', 'r','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
hold on;
surf(D, Pt, Ph_50cm_shadowavg, 'FaceColor', 'g','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
surf(D, Pt, Ph_100cm_shadowavg, 'FaceColor', 'b','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
grid on;
xlabel('传输距离 d (m)');
ylabel('传输功率 Pt (W)');
zlabel('平均收获功率 Ph (W)');
title('阴影区不同高度的平均收获功率随传输功率和距离的关系');
legend('20 cm', '50 cm', '100 cm');

%% 倾斜光路的月尘衰减系数计算
%%Ldust与月面高度H和水平距离R的关系
%定义参数
H_vdmax = 2;  % 垂直高度
H_vdmin = 1e-2;  
H_vdstep = 1e-2;  
R_vdmin = 50;  % 水平距离
R_vdmax = 1e5;  
R_vdstep = 100; 
%H和R的范围
H_vd = H_vdmin:H_vdstep:H_vdmax;
R_vd = R_vdmin:R_vdstep:R_vdmax;
[H_grid, R_grid] = meshgrid(H_vd, R_vd);
%计算每个(H, R)对应的角度 theta
theta_grid = atan(H_grid ./ R_grid);
%月尘密度函数 N(h)
Nvd = @(h) -9.5519e9*h.^5+ 5.8237e10*h.^4-1.382e11*h.^3+1.59987e11*h.^2-9.307052e10*h+2.5539e10;
Nvd_shadow=@(h) (-9.5519e9*h.^5+ 5.8237e10*h.^4-1.382e11*h.^3+1.59987e11*h.^2-9.307052e10*h+2.5539e10)*1e-4;
%消光截面 Cext
Cext_vd = Cext_avg;
L_dustvd = zeros(size(H_grid));
L_dustvd_shadow = zeros(size(H_grid));
%计算L_dust
for i = 1:numel(H_grid)
    integrand = @(h) Nvd(h) * Cext_vd / sin(theta_grid(i));
    integral_value = integral(integrand, 0, H_grid(i));
    L_dustvd(i) = exp(-integral_value);
    
    integrand_shadow = @(h) Nvd_shadow(h) * Cext_vd / sin(theta_grid(i));
    integral_value_shadow = integral(integrand_shadow, 0, H_grid(i));
    L_dustvd_shadow(i) = exp(-integral_value_shadow);
end
%绘制三维图
figure(6);
surf(H_grid, R_grid, L_dustvd);
xlabel('垂直高度 H (m)');
ylabel('水平距离 R (m)');
zlabel('L_{dust}');
title('光照区 L_{dust} 随垂直高度 H 和水平距离 R 的变化');
colorbar;  % 添加颜色条
shading interp; %去除网格线  
figure(7)
surf(H_grid, R_grid, L_dustvd_shadow);
zlim([0,1]);
xlabel('垂直高度 H (m)');
ylabel('水平距离 R (m)');
zlabel('L_{dust}');
title('阴影区 L_{dust} 随垂直高度 H 和水平距离 R 的变化');
colorbar;  % 添加颜色条
shading interp; %去除网格线  

%% 倾斜光路下 激光发射器在月表不同高度H 处的 探测器平均收获功率 随 传输功率Pt 和 收发装置水平距离R 的关系
Hgiven=[0.2;0.5;1];
%定义参数
Pt_vdmax = 1e4;  % 发射功率
Pt_vdmin = 0;  
R_vdmin = 50;  % 水平距离
R_vdmax = 1e5;  
%R和Pt的范围
R_vd = linspace(R_vdmin, R_vdmax, 1000);  % 传输距离 (从 50 到 100 km)
Pt_vd = linspace(Pt_vdmin, Pt_vdmax, 1000);  % 功率 (从 0 到 10kw)
[R_grid , Pt_grid] = meshgrid(R_vd , Pt_vd);

%20cm
%计算每个R对应的角度 theta
thetavd_20cm = atan(Hgiven(1) ./ R_grid);
%月尘密度函数 N(h)
Nvd = @(h) -9.5519e9*h.^5+ 5.8237e10*h.^4-1.382e11*h.^3+1.59987e11*h.^2-9.307052e10*h+2.5539e10;
Nvd_shadow = @(h) (-9.5519e9*h.^5+ 5.8237e10*h.^4-1.382e11*h.^3+1.59987e11*h.^2-9.307052e10*h+2.5539e10)*1e-4;
H_min_20cm=1e-2;
H_max_20cm=Hgiven(1);
H_grid_20cm=linspace(H_min_20cm, H_max_20cm, 1000);
for i = 1:numel(H_grid_20cm)
    %光照区
    integrand = @(h) Nvd(h) * Cext_vd / sin(thetavd_20cm(1,i));
    integral_value = integral(integrand, 0, H_grid_20cm(i));
    L_dustvd_20cm(i) = exp(-integral_value);
    %阴影区
    integrand_shadow = @(h) Nvd_shadow(h) * Cext_vd / sin(thetavd_20cm(1,i));
    integral_value_shadow = integral(integrand_shadow, 0, H_grid_20cm(i));
    L_dust_shadowvd_20cm(i) = exp(-integral_value_shadow);
end
D_20cm=sqrt(R_grid.^2 + Hgiven(1)^2);
theta__20cm=dr./D_20cm;%发散角
dt_20cm=1.22*lambda./theta__20cm;  % 计算发射器孔径 (随距离变化)
dr_20cm = repmat(dr, size(dt_20cm));
Gt_20cm=(pi*dt_20cm./lambda).^2;
Gr_20cm=(pi*dr_20cm./lambda).^2;
Lt_20cm=exp(-Gt_20cm*phit^2);
Lr_20cm=exp(-Gt_20cm*phir^2);
KS_20cm=lambda./(4*pi*D_20cm);
KS2_20cm=KS_20cm.^2;
%光照区
Ph_20cmavg=Pt_grid.*KS2_20cm.*(PCE*EHCE.*Lt_20cm.*Gt_20cm.*Lr_20cm.*Gr_20cm).*L_dustvd_20cm;
%阴影区
Ph_20cmavg_shadow=Pt_grid.*KS2_20cm.*(PCE*EHCE.*Lt_20cm.*Gt_20cm.*Lr_20cm.*Gr_20cm).*L_dust_shadowvd_20cm;

%50cm
%计算每个R对应的角度 theta
thetavd_50cm = atan(Hgiven(2) ./ R_grid);
%月尘密度函数 N(h)
Nvd = @(h) -9.5519e9*h.^5+ 5.8237e10*h.^4-1.382e11*h.^3+1.59987e11*h.^2-9.307052e10*h+2.5539e10;
H_min_50cm=1e-2;
H_max_50cm=Hgiven(2);
H_grid_50cm=linspace(H_min_50cm, H_max_50cm, 1000);
for i = 1:numel(H_grid_50cm)
    %光照区
    integrand = @(h) Nvd(h) * Cext_vd / sin(thetavd_50cm(1,i));
    integral_value = integral(integrand, 0, H_grid_50cm(i));
    L_dustvd_50cm(i) = exp(-integral_value);
    %阴影区
    integrand_shadow = @(h) Nvd_shadow(h) * Cext_vd / sin(thetavd_50cm(1,i));
    integral_value_shadow = integral(integrand_shadow, 0, H_grid_50cm(i));
    L_dust_shadowvd_50cm(i) = exp(-integral_value_shadow);
end
D_50cm=sqrt(R_grid.^2+Hgiven(2).^2);
theta__50cm=dr./D_50cm;%发散角
dt_50cm=1.22*lambda./theta__50cm;  % 计算发射器孔径 (随距离变化)
dr_50cm = repmat(dr, size(dt_50cm));
Gt_50cm=(pi*dt_50cm./lambda).^2;
Gr_50cm=(pi*dr_50cm./lambda).^2;
Lt_50cm=exp(-Gt_50cm*phit^2);
Lr_50cm=exp(-Gt_50cm*phir^2);
KS_50cm=lambda./(4*pi*D_50cm);
KS2_50cm=KS_50cm.^2;
%光照区
Ph_50cmavg=Pt_grid.*KS2_50cm.*(PCE*EHCE.*Lt_50cm.*Gt_50cm.*Lr_50cm.*Gr_50cm).*L_dustvd_50cm;
%阴影区
Ph_50cmavg_shadow=Pt_grid.*KS2_50cm.*(PCE*EHCE.*Lt_50cm.*Gt_50cm.*Lr_50cm.*Gr_50cm).*L_dust_shadowvd_50cm;

%100cm
%计算每个R对应的角度 theta
thetavd_100cm = atan(Hgiven(3) ./ R_grid);
%月尘密度函数 N(h)
Nvd = @(h) -9.5519e9*h.^5+ 5.8237e10*h.^4-1.382e11*h.^3+1.59987e11*h.^2-9.307052e10*h+2.5539e10;
H_min_100cm=1e-2;
H_max_100cm=Hgiven(3);
H_grid_100cm=linspace(H_min_100cm, H_max_100cm, 1000);
for i = 1:numel(H_grid_100cm)
    %光照区
    integrand = @(h) Nvd(h) * Cext_vd / sin(thetavd_100cm(1,i));
    integral_value = integral(integrand, 0, H_grid_100cm(i));
    L_dustvd_100cm(i) = exp(-integral_value);
    %阴影区
    integrand_shadow = @(h) Nvd_shadow(h) * Cext_vd / sin(thetavd_100cm(1,i));
    integral_value_shadow = integral(integrand_shadow, 0, H_grid_100cm(i));
    L_dust_shadowvd_100cm(i) = exp(-integral_value_shadow);
end
D_100cm=sqrt(R_grid.^2+Hgiven(3).^2);
theta__100cm=dr./D;%发散角
dt_100cm=1.22*lambda./theta__100cm;  % 计算发射器孔径 (随距离变化)
dr_100cm = repmat(dr, size(dt_100cm));
Gt_100cm=(pi*dt_100cm./lambda).^2;
Gr_100cm=(pi*dr_100cm./lambda).^2;
Lt_100cm=exp(-Gt_100cm*phit^2);
Lr_100cm=exp(-Gt_100cm*phir^2);
KS_100cm=lambda./(4*pi*D_100cm);
KS2_100cm=KS_100cm.^2;
%光照区
Ph_100cmavg=Pt_grid.*KS2_100cm.*(PCE*EHCE.*Lt_100cm.*Gt_100cm.*Lr_100cm.*Gr_100cm).*L_dustvd_100cm;
%阴影区
Ph_100cmavg_shadow=Pt_grid.*KS2_100cm.*(PCE*EHCE.*Lt_100cm.*Gt_100cm.*Lr_100cm.*Gr_100cm).*L_dust_shadowvd_100cm;

% 绘制三维图
figure(8);
surf(R_grid, Pt_grid, Ph_20cmavg, 'FaceColor', 'r','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
hold on;
surf(R_grid, Pt_grid, Ph_50cmavg, 'FaceColor', 'g','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
surf(R_grid, Pt_grid, Ph_100cmavg, 'FaceColor', 'b','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
grid on;
xlabel('水平距离 R (m)');
ylabel('传输功率 Pt (W)');
zlabel('平均收获功率 Ph (W)');
title('光照区不同高度的平均收获功率随传输功率和水平距离的关系');
legend('20 cm', '50 cm', '100 cm');
figure(9);
surf(R_grid, Pt_grid, Ph_20cmavg_shadow, 'FaceColor', 'r','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
hold on;
surf(R_grid, Pt_grid, Ph_50cmavg_shadow, 'FaceColor', 'g','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
surf(R_grid, Pt_grid, Ph_100cmavg_shadow, 'FaceColor', 'b','FaceAlpha', 0.5, 'EdgeColor', 'none');  % 绘制三维曲面图
grid on;
xlabel('水平距离 R (m)');
ylabel('传输功率 Pt (W)');
zlabel('平均收获功率 Ph (W)');
title('阴影区不同高度的平均收获功率Pt随传输功率Ph和水平距离R的关系');
legend('20 cm', '50 cm', '100 cm');

