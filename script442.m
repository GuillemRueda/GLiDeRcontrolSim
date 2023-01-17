% Copyright © Guillem Rueda Oller, 2022
% Author: Guillem Rueda Oller (student number 5006538)
% September 2022. Tuning Controller Parameters
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

F_max_single = 15e-3; % [N]
F_max = [F_max_single; F_max_single; F_max_single; ...
         F_max_single; F_max_single; F_max_single];

n_t = 3; % [h]

%% Scenario 1, Tunning 1
T = 30; % [s]
N_p = 20;
kappa_max = 10; % [s]
alpha = 10;
N_t = n_t*3600/T;
tic;
[d1, v1, dv1, h1, hh1] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 1, Tunning 2
T = 60; % [s]
N_p = 10;
kappa_max = 20; % [s]
alpha = 10;
N_t = n_t*3600/T;
tic;
[d2, v2, dv2, h2, hh2] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 1, Tunning 3
T = 90; % [s]
N_p = 4;
kappa_max = 30; % [s]
alpha = 10;
N_t = n_t*3600/T;
tic;
[d3, v3, dv3, h3, hh3] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 1, Tunning 4
T = 30; % [s]
N_p = 20;
kappa_max = 12; % [s]
alpha = 10;
N_t = n_t*3600/T;
tic;
[d4, v4, dv4, h4, hh4] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 1, Tunning 5
T = 30; % [s]
N_p = 20;
kappa_max = 8; % [s]
alpha = 10;
N_t = n_t*3600/T;
tic;
[d5, v5, dv5, h5, hh5] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 1, Tunning 6
T = 30; % [s]
N_p = 20;
kappa_max = 10; % [s]
alpha = 1;
N_t = n_t*3600/T;
tic;
[d6, v6, dv6, h6, hh6] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;
%% Scenario 1, Tunning 7
T = 30; % [s]
N_p = 20;
kappa_max = 10; % [s]
alpha = 1000;
N_t = n_t*3600/T;
tic;
[d7, v7, dv7, h7, hh7] = execute(T, N_p, kappa_max, alpha, d_min, ...
    L, L0, F_max, N_t, R1, R2, V1, V2, m1, m2, t_step, max_it, ...
    beta, mu, n);
toc;

%% Plot. Figure 4.4

t = 0:(T/3600):n_t;
t2 = 0:(2*T/3600):n_t;
t3 = 0:(3*T/3600):n_t;
figure;

subplot(3, 1, 1);
hold on;
plot(t, d1);
plot(t2, d2);
plot(t3, d3);
grid;
ylim0 = ylim();
ylim([0 ylim0(2)]);
xlim([0 3]);
xlabel('Time [h]')
ylabel('Distance [m]');
legend('Tunning 1: $T=30$ s, $N_p=20$, $\kappa_{max}=10$ s', ...
    'Tunning 2: $T=60$ s, $N_p=10$, $\kappa_{max}=20$ s', ...
    'Tunning 3: $T=90$ s, $N_p=7$, $\phantom{0}\kappa_{max}=30$ s', ...
    'Location', 'southeast', 'Interpreter', 'latex');

subplot(3, 1, 2);
hold on;
plot(t, d1);
plot(t, d4, 'color', [0.4940 0.1840 0.5560]);
plot(t, d5, 'color', [0.4660 0.6740 0.1880]);
grid;
ylim0 = ylim();
ylim([0 ylim0(2)]);
xlim([0 3]);
xlabel('Time [h]')
ylabel('Distance [m]');
legend('Tunning 1: $\kappa_{max}=10$ s', ...
    'Tunning 4: $\kappa_{max}=12$ s', ...
    'Tunning 5: $\kappa_{max}=8$ s', ...
    'Location', 'southeast', 'Interpreter', 'latex');

subplot(3, 1, 3);
hold on;
plot(t, d1);
plot(t, d6, 'color', [0.3010 0.7450 0.9330]);
plot(t, d7, 'color', [0.6350 0.0780 0.1840]);
grid;
ylim0 = ylim();
ylim([0 ylim0(2)]);
xlim([0 3]);
xlabel('Time [h]')
ylabel('Distance [m]');
legend('Tunning 1: $\alpha=10$', ...
    'Tunning 6: $\alpha=1$', ...
    'Tunning 7: $\alpha=1000$', ...
    'Location', 'southeast', 'Interpreter', 'latex');