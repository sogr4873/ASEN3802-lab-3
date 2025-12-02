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
panels = 1:114;
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
alphas = -15:1:15;  %range i chose    
NP = 114;           
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

   
    coeffs = polyfit((alphas(alphas >= -7 & alphas <= 7)), ...
                 Cl(alphas >= -7 & alphas <= 7, j), 1);

    fits(j).a0_panel = coeffs(1);            % slope [per degree]
    fits(j).aL0_panel = (-coeffs(2)/coeffs(1));% zero-lift angle [deg]


    [aL0_TAT, a0_TAT] = TAT(m,p,nTAT);
    fits(j).aL0_TAT = aL0_TAT;
    fits(j).a0_TAT  = a0_TAT;
end


figure(1); clf; hold on; grid on;
colors = lines(3);

for j = 1:3
  
    
    plot(alphas, Cl(:,j), 'o', 'Color', colors(j,:),'MarkerSize',4, 'DisplayName',[airfoils{j} ' (Panel)']);

    % linearfit line
    fit_line = polyval([fits(j).a0_panel, -fits(j).a0_panel*fits(j).aL0_panel], alphas);
    plot(alphas, fit_line, '-', 'Color', colors(j,:), 'LineWidth',1.5, ...
        'DisplayName',[airfoils{j} ' fit']);

    % TAT line
    cl_TAT = fits(j).a0_TAT*(alphas - fits(j).aL0_TAT);
    plot(alphas, cl_TAT, '--', 'Color', colors(j,:), 'LineWidth',1.2,'DisplayName',[airfoils{j} ' TAT']);
end
xlabel('\alpha (deg)'); ylabel('c_l');
title('Effect of Camber on Sectional Lift');
legend('Location','northwest');

hold off;

%% Part 3 Task 2
%digitize from page 462 appedix 3 in book to get experimental data


cl_exp = [-0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6];
cd_exp = [0.0118 0.0107 0.0098 0.0091 0.0088 0.0090 0.0096 0.0106 ...
          0.0120 0.0145 0.0175 0.0215 0.0265];

cd_from_cl = @(cl)interp1(cl_exp, cd_exp, cl, 'linear', 'extrap');

alphas = -10:1:10;



[x12,y12] = NACAgenerator(0,0,12, 0012, 120);

for i = 1:length(alphas)
    cl_alpha_p3(i) = Vortex_Panel(x12,y12, alphas(i));
end

cd_alpha = cd_from_cl(cl_alpha_p3);

figure(5); hold on; 

plot(cl_exp,cd_exp,'ko','MarkerFaceColor','w','DisplayName','Experimental Polar');
plot(cl_alpha_p3,cd_alpha,'r-','LineWidth',2,'DisplayName','Model');
xlabel('Section Lift Coefficient c_l');
ylabel('Section Drag Coefficient c_d');
title('NACA 0012 Drag Polar');



figure(6); hold on;

plot(alphas,cd_alpha,'b-','LineWidth',2,'DisplayName','c_d(\alpha)');
xlabel('\alpha (deg)');
ylabel('c_d');
title('Profile Drag Coefficient vs Angle of Attack');






