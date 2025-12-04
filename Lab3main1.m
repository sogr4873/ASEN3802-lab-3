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
alphas = -10:10; % Define a range of angles of attack

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

% find alpha(L=0) and a_0 based on TAT
% a_0l shoul be 0 for all since all symmetric
[a_0l0006_TAT, a00006_TAT] = TAT(m0006,p0006,1000);
[a_0l0012_TAT, a00012_TAT] = TAT(m0012,p0012,1000);
[a_0l0018_TAT, a00018_TAT] = TAT(m0018,p0018,1000);

%% Task 3
alphas = -15:1:20;  %range I chose    
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


figure(3); hold on; grid on; box on;
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

%% Part 2

% %Debug Data
% 
% b = 100; a0_t = 6.3; a0_r = 6.5; c_t = 8; c_r = 10; aero_t = 0; aero_r = -2*pi/180; geo_t = 5*pi/180; geo_r = 7*pi/180; N = 5;
% %b = 100; a0_t = 2*pi; a0_r = 2*pi; c_t = 10; c_r = 10; aero_t = 0; aero_r = 0; geo_t = 5*pi/180; geo_r = 5*pi/180; N = 5;
% 
% [e, c_L, c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
% 
% % Print results
% fprintf('e     = %.4f\n', e);
% fprintf('c_L   = %.4f\n', c_L);
% fprintf('c_Di  = %.4f\n', c_Di);
% 
% % plot for delta as a function of taper ratio
% clc;clear;
% 
% AR_lib = [4, 6, 8, 10];
% taper_ratio = linspace(0.00000001, 1, 100);
% c_r = 1;
% a0_t = 2*pi; a0_r = 2*pi; aero_t = 0; aero_r = 0; geo_t = deg2rad(5); geo_r = geo_t; N = 50;
% 
% 
% for k = 1:length(AR_lib)
%     AR = AR_lib(k);
% 
%     for j = 1:length(taper_ratio)
%         c_t = taper_ratio(j)*c_r;
% 
%         b = AR*(c_r+c_t)*0.5;
%         [~,~,~,delta] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);
%         delta_vec(j) = delta;
%     end
% 
%     delta_matrix(k,:) = delta_vec;
% end
% 
% figure(4); hold on;
% xlim([0,1]);
% ylim([0,0.2]);
% for k = 1:length(AR_lib)
%     plot(taper_ratio, delta_matrix(k,:), 'LineWidth', 2)
% end
% 
% xlabel('Taper Ratio \lambda');
% ylabel('\delta (Induced Drag Factor)');
% legend('AR = '+string(AR_lib), 'Location','best');
% grid on;
% title('\delta vs \lambda with varying AR');

%% Part 3

%% task 1

b_p3 = 36;
cr_p3 = 5 + (4/12);
ct_p3 = 3 + (7/12);
N=5;
% root is 2412 tip is 0012

% first 2 rows of fit struct already have AoA for 0012 and 2412 from task 3
for i = 1:length(alphas)
    [e_p3(i),c_L_p3(i),c_Di_p3(i), delta_out_p3(i)] = PLLT(b_p3, fits(1).a0_panel, fits(2).a0_panel, ct_p3, cr_p3, fits(1).aL0_panel * (pi / 180), fits(2).aL0_panel * (pi / 180), alphas(i) * (pi /180), alphas(i) + (2*pi / 180), N);
end

figure(4);
plot(alphas, c_L_p3);
title('Coefficient of Lift vs. Angle of Attack for the Cessna 180');
xlabel('\alpha (deg)')
ylabel('C_L');


%% Task 2
%digitize from page 462 appedix 3 in book to get experimental data
cl_exp = readmatrix("NACA0012_cl_alpha_TESTALL.csv"); %cl vs alpha experimental data
cd_exp = readmatrix("NACA0012_cd_cl.csv"); %cd vs cl

cd_from_cl = @(cl)interp1(cd_exp(:,1), cd_exp(:,2), cl, 'linear', 'extrap'); % continuous cd vs cl func
cl_from_alpha = @(a)interp1(cl_exp(:,1), cl_exp(:,2), a, 'linear', 'extrap'); % continuous cl vs alpha func
cd_from_alpha = @(a)  cd_from_cl( cl_from_alpha(a) );
cd_from_alpha_exp = cd_from_alpha(cl_exp(:,1));
cd_from_alpha_t3 = cd_from_alpha(alphas);

alphas = -15:1:20;

% figure(5); 
% hold on; 
% plot(cd_exp(:,1),cd_exp(:,2),'DisplayName','Experimental Polar');
% plot(c_L_p3,cd_cl,'r-','LineWidth',2,'DisplayName','Model');
% xlabel('Section Lift Coefficient c_l');
% ylabel('Section Drag Coefficient c_d');
% title('NACA 0012 Drag Polar');
% legend
% hold off;

figure(6); 
hold on;
plot(cl_exp(:,1), cd_from_alpha_exp,'k','DisplayName','Model');
fplot(cd_from_alpha, 'bo', 'DisplayName', 'Experimental Data');
xlabel('\alpha (deg)');
ylabel('C_D_0');
title('Profile Drag Coefficient vs Angle of Attack');
hold off;

%% task 3
y = linspace(-b_p3/2,b_p3/2); %span
c_y = cr_p3 - (((cr_p3 - ct_p3)/(b_p3/2)).*(abs(y))); %chord func spanwise

S = (1/2) * (cr_p3 + ct_p3) * b_p3;
for i = 1:length(alphas)
    CD_0(i) = trapz(y,c_y.*cd_from_alpha_t3(i)) / S;
end

CD_p3 = CD_0 + c_Di_p3;

figure(7);
hold on;
plot(alphas, CD_p3,LineWidth=2);
plot(alphas, c_Di_p3,LineWidth=2);
plot(alphas, CD_0,LineWidth=2);
hold off;
title('Drag Coefficients vs. AoA');
xlabel('\alpha (deg)');
ylabel('C_D');
legend('Total (C_D)', 'Induced (C_D_i)', 'Profile (C_D_0)');

%% task 2 bonus: c_d is varying spanwise

N = 50;
y_bt2 = linspace(-b_p3/2,b_p3/2,N); %span
c_y_bt2 = cr_p3 - (((cr_p3 - ct_p3)/(b_p3/2)).*(abs(y_bt2))); %chord func spanwise
S = (1/2) * (cr_p3 + ct_p3) * b_p3;

mtot_ec = ((0 - 2) / (b_p3 / 2)) .* y_bt2 + 2;
ptot_ec = ((0 - 4) / (b_p3 / 2)) .* y_bt2 + 4;
ttot_ec = ((12 - 12) / (b_p3 / 2)) .* y_bt2 + 12; % creating linear distribution of airfoils

for i = 1:length(y_bt2)
    [xec,yec] = NACAgenerator(mtot_ec(i), ptot_ec(i), ttot_ec(i), 200, 120); % approx sectional airfoil shape along span
    for j = 1:length(alphas)
        [Cl_y(i,j)] = Vortex_Panel(xec,yec,alphas(j)); % find sectional cl for sectional airfoil
        cd_y(i,j) = cd_from_cl(Cl_y(i,j)); % find sectional cd from sectional cl
    end
end

CD_0_new = trapz(y_bt2, c_y_bt2' .* cd_y, 1) ./ S; % integrate along columns of cd_y to get one value of wing profile drag for each AoA

figure(8); 
hold on; 
plot(y_bt2,cd_y(:,16),'r-','LineWidth',2);
xlabel('y');
ylabel('c_d(y) ');
title('Sectional Drag Coefficent varying with Span');
hold off;

figure(9); 
hold on;
plot(alphas,CD_0_new,'b-','LineWidth',2,'DisplayName','c_d(\alpha)');
xlabel('\alpha (deg)');
ylabel('CD_0');
title('Profile Drag vs Angle of Attack');
hold off;

figure(6);
hold on;
plot(alphas,CD_0_new,'r-','LineWidth',2,'DisplayName','C_D_0(\alpha) Varying Span Distribution');
legend
hold off;

%% task 4
% thrust must equal profile drag, lift must equal weight, airspeed effects
% cl, cd
% altitude is at 10,000 ft, weight is 2500 lbs
S_t4 = (1/2) * (cr_p3 + ct_p3) * b_p3;

weight = 2500;

rho_10k = 17.56*10^-4; %slugs/ft^3
vel = zeros(1,length(c_L_p3));
thrust = zeros(1,length(c_L_p3));
for i = 1:length(c_L_p3)
    vel(i) = sqrt((2*2500)/(rho_10k*c_L_p3(i)*S_t4));
    thrust(i) = CD_p3(i)*(0.5)*((vel(i))^2)*rho_10k*S_t4;
end

vel_knots = vel.*(1/1.688);
vel_knots = vel_knots(17:26);
thrust = thrust(17:26);

figure(15);
grid on;
plot(vel_knots,thrust, 'Color','red',LineWidth=1.2);
xlabel('Air Speed [knots]');
ylabel('Thrust Required [lbs]');
title('Thrust Required for Steady Level Flight vs Airspeed');

