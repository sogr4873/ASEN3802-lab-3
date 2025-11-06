function [x,y] = NACAgenerator(mtot, ptot, ttot, fignum)
m = mtot / 100;
p = ptot / 10;
t = ttot / 100;

xoc = linspace(0,1,100);
yt = ( t / 0.2) .* ( 0.2969 .* sqrt(xoc) - 0.1260 .* xoc - 0.3516 .* xoc.^2 + 0.2843 .* xoc.^3 - 0.1036 .* xoc.^4);

for i = 1: length(xoc)
    if (xoc(i) >= 0) && (xoc(i) < p)
        yc(i) = m .* (xoc(i) ./ p^2) * (2*p - xoc(i));
        dycdx(i) = ( ( 2 * m ) / p^2 ) * ( p - xoc(i) );
    else 
        yc(i) = m .* ( (1 - xoc(i) ) ./ (1 - p)^2 ) .* (1 + xoc(i) - 2*p);
        dycdx(i) = ( ( 2 * m ) / ( 1 - p)^2 ) * ( p - xoc(i));
    end
end

zeta = atan(dycdx);

xu = xoc - yt .* sin(zeta);
xl = xoc + yt .* sin(zeta);
yu = yc + yt .* cos(zeta);
yl = yc - yt .* cos(zeta);

x(:,1) = [xl(end:-1:1) xu];
y(:,1) = [yl(end:-1:1) yu];

figure(fignum);
plot(x,y);
hold on
plot(xoc, yc);
xlabel('x (% chord)');
ylabel('y % chord');
title(["Shape of NACA" + num2str(mtot) + num2str(ptot) + num2str(ttot) + " Airfoil"]);


end