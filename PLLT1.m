
clc;
clear all

 b = 100; a0_t = 6.3; a0_r = 6.5; c_t = 8; c_r = 10; aero_t = 0; aero_r = -2*pi/180; geo_t = 5*pi/180; geo_r = 7*pi/180; N = 5;

[e, c_L, c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

% Print results
fprintf('e     = %.4f\n', e);
fprintf('c_L   = %.4f\n', c_L);
fprintf('c_Di  = %.6f\n', c_Di);


function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)

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
n = 1:1:N;
delta = sum(n .* ((A(2:end) ./ A(1)).^2));
c_Di = ( c_L^2 / (pi * AR)) * (1 + delta);
e = 1 ./ (1+ delta);


end