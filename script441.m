% Copyright © Guillem Rueda Oller, 2022
% Author: Guillem Rueda Oller (student number 5006538)
% September 2022. Multiple Scenarios
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

t_step = 0.001; % [s]
max_it = 20;
beta = 1e-20;

T = 30; % [s]
N_p = 20;
kappa_max = 10; % [s]
alpha = 10;

F_max_single = 15e-3; % [N]
F_max = [F_max_single; F_max_single; F_max_single; ...
         F_max_single; F_max_single; F_max_single];

n_t = 24; % [h]
N_t = n_t*3600/T;

%% Scenario 1
L0 = 25; % [m]
L = 20; % [m]
d_min = 15; % [m]
tic;
[d1, v1, dv1, h1, hh1] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 2
L0 = 20; % [m]
L = 20; % [m]
d_min = 15; % [m]
tic;
[d2, v2, dv2, h2, hh2] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 3
L0 = 15; % [m]
L = 20; % [m]
d_min = 15; % [m]
tic;
[d3, v3, dv3, h3, hh3] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 4
L0 = 40; % [m]
L = 35; % [m]
d_min = 15; % [m]
tic;
[d4, v4, dv4, h4, hh4] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 5
L0 = 35; % [m]
L = 35; % [m]
d_min = 15; % [m]
tic;
[d5, v5, dv5, h5, hh5] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 6
L0 = 30; % [m]
L = 35; % [m]
d_min = 15; % [m]
tic;
[d6, v6, dv6, h6, hh6] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;

%% Plots

t = 0:(T/3600):n_t;

%% Figure 4.1
figure;
hold on;
plot(t, d1);
plot(t, d2, '--');
plot(t, d3);
plot(t, d4);
plot(t, d5, '--');
plot(t, d6);
grid;
ylim0 = ylim();
ylim([0 ylim0(2)]);
xlim([0 3]);
xlabel('Time [h]')
ylabel('Distance [m]');
lgd = legend('Scenario 1', 'Scenario 2', 'Scenario 3', ...
    'Scenario 4', 'Scenario 5', 'Scenario 6', 'Location', 'south', ...
    'Orientation', 'horizontal');
lgd.NumColumns = 3;

%% Figure 4.2
figure;
subplot(2,2,1);
hold on;
plot(t, h2, '--', 'color', [0.850, 0.325, 0.098]);
plot(t, h5, '--', 'color', [0.466, 0.674, 0.188]);
grid;
xlim([0 n_t]);
ylim1 = ylim();
ylim([0 ylim1(2)]);
xticks(0:6:n_t);
xlabel('Time [h]')
ylabel('$\Delta{r}$ [m]', 'Interpreter', 'latex');

subplot(2,2,3);
hold on;
plot(t(1:(end-1)), (h2(2:end)-h2(1:(end-1)))/(T/3600), '--', ...
    'color', [0.850, 0.325, 0.098]);
plot(t(1:(end-1)), (h5(2:end)-h5(1:(end-1)))/(T/3600), '--', ...
    'color', [0.466, 0.674, 0.188]);
grid;
xlim([0 n_t]);
ylim2 = ylim();
ylim([0 ylim2(2)]);
xticks(0:6:n_t);
xlabel('Time [h]')
ylabel('$\dot{r}$ [m/h]', 'Interpreter', 'latex');

subplot(2,2,2);
hold on;
plot(t, hh2, '--', 'color', [0.850, 0.325, 0.098]);
plot(t, hh5, '--', 'color', [0.466, 0.674, 0.188]);
grid;
xlim([0 n_t]);
ylim([0 ylim1(2)]);
xticks(0:6:n_t);
ylim2 = ylim();
ylim([0 ylim2(2)]);
xlabel('Time [h]')
ylabel('$\Delta{a}$ [m]', 'Interpreter', 'latex');

subplot(2,2,4);
hold on;
plot(t(1:(end-1)), (hh2(2:end)-hh2(1:(end-1)))/(T/3600), '--', ...
    'color', [0.850, 0.325, 0.098]);
plot(t(1:(end-1)), (hh5(2:end)-hh5(1:(end-1)))/(T/3600), '--', ...
    'color', [0.466, 0.674, 0.188]);
grid;
xlim([0 n_t]);
ylim([0 ylim2(2)]);
xticks(0:6:n_t);
ylim2 = ylim();
ylim([0 ylim2(2)]);
xlabel('Time [h]')
ylabel('$\dot{a}$ [m/h]', 'Interpreter', 'latex');

%% Figure 4.3
figure;
subplot(2, 1, 1);
hold on;
plot(t(1:(end-1)), v1);
plot(t(1:(end-1)), v2, '--');
plot(t(1:(end-1)), v3);
plot(t(1:(end-1)), v4);
plot(t(1:(end-1)), v5, '--');
plot(t(1:(end-1)), v6);
grid;
xlim([0 n_t]);
xticks(0:2:n_t);
xlabel('Time [h]')
ylabel('Cumulative $\left\|\Delta{V}\right\|$ [m/s]', ...
    'Interpreter', 'latex');
lgd = legend('Scenario 1', 'Scenario 2', 'Scenario 3', ...
    'Scenario 4', 'Scenario 5', 'Scenario 6', 'Location', 'north', ...
    'Orientation', 'horizontal');
lgd.NumColumns = 3;

subplot(2, 1, 2);
hold on;
plot(t(1:(end-1)), sum(dv1)/(T/3600));
plot(t(1:(end-1)), sum(dv2)/(T/3600), '--');
plot(t(1:(end-1)), sum(dv3)/(T/3600));
plot(t(1:(end-1)), sum(dv4)/(T/3600));
plot(t(1:(end-1)), sum(dv5)/(T/3600), '--');
plot(t(1:(end-1)), sum(dv6)/(T/3600));
grid;
xlim([0 n_t]);
xticks(0:2:n_t);
xlabel('Time [h]')
ylabel('$\dot{\left\|\Delta{V}\right\|}$ [(m/s)/h]', ...
    'Interpreter', 'latex');
lgd = legend('Scenario 1', 'Scenario 2', 'Scenario 3', ...
    'Scenario 4', 'Scenario 5', 'Scenario 6', 'Location', 'north', ...
    'Orientation', 'horizontal');
lgd.NumColumns = 3;