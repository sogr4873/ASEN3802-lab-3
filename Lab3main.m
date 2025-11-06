clc
clear
close all;

%% Task 1
% NACA 0012
m0012 = 0;
p0012 = 0;
t0012 = 12;
[x0012,y0012] = NACAgenerator(m0012,p0012,t0012,0012, 8);

% NACA 0018
m0018 = 0;
p0018 = 0;
t0018 = 18;
[x0018,y0018] = NACAgenerator(m0018, p0018, t0018, 0018, 8);

% NACA 2418
m2418 = 2;
p2418 = 4;
t2418 = 18;
[x2418,y2418] = NACAgenerator(m2418, p2418, t2418, 2418, 8);

% Call Vortex Panel

Predict_Cl = Vortex_Panel(x0012,y0012,5);

TAT_Cl = @(a) 2*pi*(a*pi/180);
goal = TAT_Cl(5);
panels = 1:100;
panel_crit = 0;
for p = panels
    [testx0012,testy0012] = NACAgenerator(m0012,p0012,t0012,0012, p)
    test_cl = Vortex_Panel(testx0012,testy0012,5);
    if abs(test_cl - goal) < 0.01
        panel_crit = p;
        break;
    end
end



