%% Debug Data

clc;
clear;

b = 100; a0_t = 6.3; a0_r = 6.5; c_t = 8; c_r = 10; aero_t = 0; aero_r = -2*pi/180; geo_t = 5*pi/180; geo_r = 7*pi/180; N = 5;
%b = 100; a0_t = 2*pi; a0_r = 2*pi; c_t = 10; c_r = 10; aero_t = 0; aero_r = 0; geo_t = 5*pi/180; geo_r = 5*pi/180; N = 5;

[e, c_L, c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

% Print results
fprintf('e     = %.4f\n', e);
fprintf('c_L   = %.4f\n', c_L);
fprintf('c_Di  = %.4f\n', c_Di);

%% plot for delta as a function of taper ratio
clc;clear;

AR_lib = [4, 6, 8, 10];
taper_ratio = linspace(0.00000001, 1, 100);
c_r = 1;
a0_t = 2*pi; a0_r = 2*pi; aero_t = 0; aero_r = 0; geo_t = deg2rad(5); geo_r = geo_t; N = 50;


for k = 1:length(AR_lib)
    AR = AR_lib(k);

    for j = 1:length(taper_ratio)
        c_t = taper_ratio(j)*c_r;
        
        b = AR*(c_r+c_t)*0.5;
        [~,~,~,delta] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);
        delta_vec(j) = delta;
    end

    delta_matrix(k,:) = delta_vec;
end

figure; hold on;
xlim([0,1]);
ylim([0,0.2]);
for k = 1:length(AR_lib)
    plot(taper_ratio, delta_matrix(k,:), 'LineWidth', 2)
end

xlabel('Taper Ratio \lambda');
ylabel('\delta (Induced Drag Factor)');
legend('AR = '+string(AR_lib), 'Location','best');
grid on;
title('\delta vs \lambda with varying AR');



function [e,c_L,c_Di, delta_out] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)

for i = 1:1:N
    theta(i) = (i * pi) / (2 * N);
    y(i) = (-b / 2) * cos(theta(i));
    c(i) = ((c_r - c_t) / (b / 2)) * y(i) + c_r;
    a0(i) = ((a0_r - a0_t) / (b / 2)) * y(i) + a0_r;
    aero(i) = ((aero_r - aero_t) / (b / 2)) * y(i) + aero_r; % create vector of properties that vary with theta

    for j = 1:1:N
        M_ij(i,j) = ((4 .* b) ./ (a0(i) .* c(i))) .* sin((2.*j - 1) .* theta(i)) + (2.*j - 1) .* ((sin((2.*j - 1) .* theta(i))) ./ sin(theta(i)));
       
    end
    alpha_geo(i) = ((geo_r - geo_t)/(b/2))*y(i) + geo_r;%using geometric angle of attack
    alpha_L0(i)  = aero(i);   

d_i(i) = alpha_geo(i) - alpha_L0(i);
    
end
A = M_ij \ d_i';
S = (1/2) * (c_r + c_t) * b;
AR = (b^2 / S);
c_L = A(1) * pi * AR;
delta = 0;
for n = 2:N
    m = 2*n - 1;
    delta = delta + m*((A(n)/A(1))^2);
end
c_Di = ( c_L^2 / (pi * AR)) * (1 + delta);
e = 1 / (1 + delta);

delta_out = delta;


end
