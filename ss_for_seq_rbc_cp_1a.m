% ss_for_seq_rbc_cp_1a

% Miguel Alvarado

% Para hallar el estado estacionario del sistema de ecuaciones definido
% en el m.file "sist_for_ss_rbc_cp_1a.m" 

clear;
clc;

global p_beta p_delta p_alfa p_sigma p_v p_rho ...
       c_ss k_ss z_ss ...
       g h p_theta p_pi;
 
% Definimos el vector de parámetros P

P = [p_beta; p_delta; p_alfa; p_sigma; p_v];

%p_rho   = (1/beta)-1;

% Hallando el estado estacionario

opt = optimset('Display','iter');

% Tomaremos un valor inicial (guess) para todas las variables del modelo

X0=[0.5;
    0.5;
    0.5];

[solstar] = fsolve(@(X) sist_for_ss_rbc_cp_1a(X,P),X0,opt);

% La solución hallada será, considerando la construcción de X en la función
% sist_for_ss_rbc_cp_1, se tiene:

c_ss = solstar(1,1);
k_ss = solstar(2,1);
z_ss = solstar(3,1);

display (c_ss);
display (k_ss);
display (z_ss);

