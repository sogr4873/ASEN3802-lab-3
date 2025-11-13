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

%% Task 3
alphas = -10:1:10;  %range I chose    
NP = 100;           
nTAT = 100;           % Integration points for TAT? might change 
airfoils = {'NACA 0012','NACA 2412','NACA 4412'};
geom = [0 0 12; 2 4 12; 4 4 12];  % [m p t]

Cl = zeros(numel(alphas),3);  % store Cl 
fits = struct([]);

for j = 1:3
    % Geometry
    m = geom(j,1); p = geom(j,2); t = geom(j,3);
    [x,y] = NACAgenerator(m,p,t,j,NP);

   
    for k = 1:numel(alphas)
        Cl(k,j) = Vortex_Panel(x,y,alphas(k));
    end

   
   coeffs = polyfit(alphas(alphas >= -5 & alphas <= 5), ...
                 Cl(alphas >= -5 & alphas <= 5, j), 1);

    fits(j).a0_panel = coeffs(1);            % slope [per degree]
    fits(j).aL0_panel = -coeffs(2)/coeffs(1);% zero-lift angle [deg]


    [aL0_TAT, a0_TAT] = TAT(m,p,nTAT);
    fits(j).aL0_TAT = aL0_TAT;
    fits(j).a0_TAT  = a0_TAT;
end


figure(1); hold on; grid on; box on;
colors = lines(3);

for j = 1:3
  
    plot(alphas, Cl(:,j), 'o', 'Color', colors(j,:), ...
        'MarkerSize',4, 'DisplayName',[airfoils{j} ' (Panel)']);

    % Linear-fit line
    fit_line = polyval([fits(j).a0_panel, ...
                        -fits(j).a0_panel*fits(j).aL0_panel], alphas);
    plot(alphas, fit_line, '-', 'Color', colors(j,:), 'LineWidth',1.5, ...
        'DisplayName',[airfoils{j} ' fit']);

    % Thin Airfoil Theory line
    cl_TAT = fits(j).a0_TAT*(alphas - fits(j).aL0_TAT);
    plot(alphas, cl_TAT, '--', 'Color', colors(j,:), 'LineWidth',1.2, ...
        'DisplayName',[airfoils{j} ' TAT']);
end
xlabel('\alpha (deg)'); ylabel('c_l');
title('Effect of Camber on Sectional Lift');
legend('Location','northwest');


