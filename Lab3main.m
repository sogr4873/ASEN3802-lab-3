clc
clear
close all;



%% Task 1
% NACA 0018
m0018 = 0;
p0018 = 0;
t0018 = 18;
numpanels0018 = 120;
[x0018,y0018] = NACAgenerator(m0018, p0018, t0018, 0018, numpanels0018);

% NACA 2418
m2418 = 2;
p2418 = 4;
t2418 = 18;
numpanels2418 = 120;
[x2418,y2418] = NACAgenerator(m2418, p2418, t2418, 2418, numpanels2418);

%% Task 2
alpha = 5;
m0012 = 0;
p0012 = 0;
t0012 = 12;
%[x0012,y0012] = NACAgenerator(m0012,p0012,t0012,0012, 10);
%cl0012 = Vortex_Panel(x0012,y0012,alpha);

% finding the number of panels needed for less than 1% error
[exactx0012,exacty0012] = NACAgenerator(m0012,p0012,t0012,0012, 1000); 
exactcl0012 = Vortex_Panel(exactx0012,exacty0012,alpha); % finding "exact" coeff of lift (high panel #)
panels = 2:200;
for i = 1:length(panels)
    p = panels(i);
    [testx0012,testy0012] = NACAgenerator(m0012,p0012,t0012,0012, p);
    test_cl(i) = Vortex_Panel(testx0012,testy0012,alpha); % cl at various panel numbers
    if (abs(exactcl0012 - test_cl(i)) / exactcl0012) < 0.01 % relative error to exact coeff
        panel_crit = p; % break when error less than 1%
        % save # of panels needed to reach 1% error
        break;
    end
end

% deliverable 1: plotting cl vs number of panels

for j = 1:length(panels)
    [x0012,y0012] = NACAgenerator(m0012,p0012,t0012,0012, panels(j));
    estimatecl0012(j) = Vortex_Panel(x0012,y0012,alpha); % finding cl at all the panel #s
end

figure(1);
plot(panels, estimatecl0012);
hold on;
yline(exactcl0012);
yline(exactcl0012 - exactcl0012*0.01, 'Color','r');
xlabel('Number of Panels');
ylabel('Sectional Coefficient of Lift (c_L)');
ylim([0.5,0.7]);
title('Predicted c_L vs Number of Vortex Panels for Naca0012 Airfoil');
legend('c_L Predicted', '"Exact" c_L', '1% Error');

% Deliverable 2: Plot cL vs alpha for 3 airfoils
alphas = 0:10; % Define a range of angles of attack

m0006 = 0;
p0006 = 0;
t0006 = 6;

[x0006,y0006] = NACAgenerator(m0006,p0006,t0006,0006, 120);
[x0012,y0012] = NACAgenerator(m0012,p0012,t0012,0012, 120);

for k = 1:length(alphas)
    % NACA0006
    cl0006(k) = Vortex_Panel(x0006,y0006,alphas(k)); % finding cl at all the alpha values

    % NACA0012
    cl0012(k) = Vortex_Panel(x0012,y0012,alphas(k));

    % NACA0018
    cl0018(k) = Vortex_Panel(x0018,y0018,alphas(k));
end

figure(2);
plot(alphas, cl0006);
hold on;
plot(alphas, cl0012);
plot(alphas, cl0018);
xlabel('\alpha (deg)');
ylabel('c_L');
title('Coefficienct of Lift vs. Angle of Attack');
legend('NACA0006', 'NACA0012', 'NACA0018');

% find alpha(L=0) and a_0 based on plots

junk1 = polyfit(alphas, cl0006,1);
a0_0006 = junk1(1); % lift slope from polyfit
aL0_0006 = 0; % airfoils are symmetric so 0 lift at 0 AoA

junk2 = polyfit(alphas, cl0012,1);
a0_0012 = junk2(1);
aL0_0012 = 0;

junk3 = polyfit(alphas, cl0018,1);
a0_0018 = junk3(1);
aL0_0018 = 0;

%% find alpha(L=0) and a_0 based on TAT
% a_0l shoul be 0 for all since all symmetric
[a_0l0006_TAT, a00006_TAT] = TAT(m0006,p0006,1000);
[a_0l0012_TAT, a00012_TAT] = TAT(m0012,p0012,1000);
[a_0l0018_TAT, a00018_TAT] = TAT(m0018,p0018,1000);
