clc
clear
close all;

%% Task 1
% NACA 0018
m0018 = 0;
p0018 = 0;
t0018 = 18;
[x0018,y0018] = NACAgenerator(m0018, p0018, t0018, 0018);

% NACA 2418
m2418 = 2;
p2418 = 4;
t2418 = 18;
[x2418,y2418] = NACAgenerator(m2418, p2418, t2418, 2418);



