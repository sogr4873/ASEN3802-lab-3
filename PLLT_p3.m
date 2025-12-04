function [e,c_L,c_Di, delta_out,cl_y] = PLLT_p3(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N,c_y)

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

% Compute spanwise sectional lift coefficient cl_y
cl_y = zeros(1,N);
for i = 1:N
    sum_n = 0;
    for n = 1:N
        m = 2*n - 1;
        sum_n = sum_n + A(n) * sin(m * theta(i));
    end
    cl_y(i) = (4.*b ./ c_y(i)) .* sum_n;
end



end
