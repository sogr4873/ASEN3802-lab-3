function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)

for i = 1:1:N
    theta(i) = (i * pi) / (2 * N);
    y(i) = (-b / 2) * cos(theta(i));
    c(i) = ((c_r - c_t) / (b / 2)) * y(i) + c_r;
    a0(i) = ((a0_r - a0_t) / (b / 2)) * y(i) + a0_r;
    aero(i) = ((aero_r - aero_t) / (b / 2)) * y(i) + aero_r; % create vector of properties that vary with theta

    for j = 1:1:N
        M_ij(i,j) = ((4 .* b) ./ (aero(i) .* c(i))) .* sin((2.*j - 1) .* theta(i)) + (2.*j - 1) .* ((sin((2.*j - 1) .* theta(i))) ./ sin(theta(i)));
       
    end
    d_i(i) = a0(i) - aero(i);
    
end
A = M_ij \ d_i';
S = (1/2) * (c_r + c_t) * b;
AR = (b^2 / S);
c_L = A(1) * pi * AR;
n = 1:1:N;
delta = sum(n .* ((A(2:end) ./ A(1)).^2));
c_Di = ( c_L^2 / (pi * AR)) * (1 + delta);
e = 1 / (1+ delta);
end