function [a_0l , a0] = TAT(mtot, ptot, n)
xoc_vec = linspace(0,1,n);
theta_vec = acos(1-2.*xoc_vec);
dzdx_conv = size(xoc_vec,2);

m = mtot / 100;
p = ptot / 10;

for i = 1:length(xoc_vec)
    if (xoc_vec(i) >= 0) && (xoc_vec(i) < p)
        dzdx_conv(i) = ( ( 2 * m ) / p^2 ) * ( p - (0.5.*(1-(cos(theta_vec(i))))));
    else 
        dzdx_conv(i) = ( ( 2 * m ) / ( 1 - p)^2 ) * ( p - (0.5.*(1-(cos(theta_vec(i))))));
    end
end

% find a_0l
integrand = dzdx_conv.*(cos(theta_vec) - 1);
a_0l = (-1/pi).*trapz(theta_vec,integrand);
a_0l = rad2deg(a_0l);
% a0 is 2pi, conv to per degree from per rad
a0 = (2*pi)*(pi/180);
end