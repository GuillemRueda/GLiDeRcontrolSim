% Copyright © Guillem Rueda Oller, 2022
% Author: Guillem Rueda Oller (student number 5006538)
% September 2022. Thrusters Configurations
close all;
clear all; %#ok<CLALL>
clc;

Tday = 86164.0905;
mu = 3.986044418e14;
n = 2*pi/Tday;

R1 = 3; % [m]
R2 = 3; % [m]
V1 = 20e3; % [V]
V2 = -20e3; % [V]
m1 = 500; % [kg]
m2 = 1000; % [kg]

t_step = 0.1; % [s]
max_it = 20;
beta = 1e-20;

L0 = 25; % [m]
L = 20; % [m]
d_min = 15; % [m]

T = 30; % [s]
N_p = 20;
kappa_max = 10; % [s]
alpha = 10;

F_max_single = 15e-3; % [N]

n_t = 3; % [h]
N_t = n_t*3600/T;

%% Scenario 1, Configuration 1
F_max = [F_max_single; F_max_single; F_max_single; ...
         F_max_single; F_max_single; F_max_single];
tic;
[d1, v1, dv1, h1, hh1] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 1, Configuration 2
F_max = [0; F_max_single; 0; ...
         0; F_max_single; 0];
tic;
[d2, v2, dv2, h2, hh2] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 1, Configuration 3
F_max = [0; F_max_single; 0; ...
         0;            0; 0];
tic;
[d3, v3, dv3, h3, hh3] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 1, Configuration 4
F_max = [0; 1e-3; 0; ...
         0;    0; 0]; % [N]
tic;
[d4, v4, dv4, h4, hh4, i] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;

%% Plots

t = 0:(T/3600):n_t;

%% Figure 4.5

figure;
hold on;
plot(t, d1);
plot(t, d2);
plot(t, d3);
plot(t(1:(i-1)), d4(1:(i-1)), '--');
plot(t(i), d4(i), 'x', 'color', 'black','LineWidth', 2);
grid;
ylim0 = ylim();
ylim([0 30]);
xlim([0 3]);
xlabel('Time [h]')
ylabel('Distance [m]');
legend('Configuration 1: 3D thrusting (15 mN)', ...
    'Configuration 2: $\pm y$ thrusting (15 mN)', ...
    'Configuration 3: $+y$ thrusting (15 mN)', ...
    'Configuration 4: $+y$ thrusting (1 mN)', ...
    '\it{Config. 4 stops}', 'Location', 'south', ...
    'Interpreter', 'latex');

%% Figure 4.6
figure;
subplot(2, 1, 1);
hold on;
plot(t(1:(end-1)), v1);
plot(t(1:(end-1)), v2);
plot(t(1:(end-1)), v3);
grid;
xlim([0 n_t]);
xlabel('Time [h]')
ylabel('Cumulative $\left\|\Delta{V}\right\|$ [m/s]', ...
    'Interpreter', 'latex');
lgd = legend('Config. 1', 'Config. 2', 'Config. 3', ...
    'Location', 'north', 'Orientation', 'horizontal');
lgd.NumColumns = 3;

subplot(2, 1, 2);
hold on;
plot(t(1:(end-1)), sum(dv1)/(T/3600));
plot(t(1:(end-1)), sum(dv2)/(T/3600));
plot(t(1:(end-1)), sum(dv3)/(T/3600));
grid;
xlim([0 n_t]);
xlabel('Time [h]')
ylabel('$\dot{\left\|\Delta{V}\right\|}$ [(m/s)/h]', ...
    'Interpreter', 'latex');
lgd = legend('Config. 1', 'Config. 2', 'Config. 3', ...
    'Location', 'north', 'Orientation', 'horizontal');
lgd.NumColumns = 3;