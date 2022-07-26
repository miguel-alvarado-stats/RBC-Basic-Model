% rbc_cp_1a

% Miguel Alvarado

clear global;
clear all;
clc;
% D:\matlab_my_files\malvarado_programs\ejercicio1_a;
format compact;
format short;

global p_beta p_delta p_alfa p_sigma p_v p_rho ...
       c_ss k_ss z_ss ...
       g h p_theta p_pi;

% parametros

p_beta  = 0.9825;
p_delta = 0.025;
p_alfa  = 1/3;
p_sigma = 2;
p_v     = 1;
p_rho   = 0.9;

% Calculo Estado EStacionario (ver archivo ss_for_seq_rbc_cp_1.m)

run ss_for_seq_rbc_cp_1a.m

g = ((1-p_alfa)/(p_v+p_alfa));
h = ((p_v+1)/(p_v+p_alfa));

% En el modelo log-linealizado (ver archivo pdf, sistema 2), denotamos los 
% siguientes parámetros

p_theta = ((1-p_alfa)^g)*(z_ss^h)*(c_ss^(-p_sigma*g))*(k_ss^(p_alfa*h)); % theta
p_pi = p_alfa*p_beta*((1-p_alfa)^g)*(z_ss^h)*(c_ss^(-p_sigma*g))*(k_ss^(-p_v*g)); % pi

a1 = (p_sigma+(p_pi*p_sigma*g));
a2 = (p_pi*p_v*g);
a3 = -(p_pi*h);

b1 = ((p_theta*p_sigma*g)+c_ss)/k_ss;
b2 = -((p_theta*p_alfa*h)+((1-p_delta)*k_ss))/k_ss;
b3 = -(p_theta*h)/k_ss;

% modelo matricial;

A1 = [p_sigma 0 0;
      b1 b2 b3;
      0 0 p_rho];
      
A2 = [a1 a2 a3;
      0 -1 0;
      0 0 1];
  
A3 = [0 a1 a2 a3;
      0 0 0 0;
      -1 0 0 0];
  
% formal estructural  
 
A1_inv = inv(A1);

A = A1_inv*A2;

B = A1_inv*A3;

[Q,F] = jordan(A); %F matriz de autovalores de A, Q matriz de autovectores de A

[U,T] = schur(A);  %T matriz de autovalores de A, U matriz de autovectores de A

% Verificamos las condiciones de Blanchard & Kahn (estabilidad del sistema)

variables = 3;
vcontrol = 1; %variablles libres {c_(t)}

root_g_1=0; % contador
root_l_1=0; % contador

for i=1:variables
if abs(F(i,i))>=1
    root_g_1=root_g_1 + 1;
else
    root_l_1=root_l_1 + 1;
end
end;

if root_l_1==vcontrol
    display('SISTEMA ESTABLE'), F
else
    display('SISTEMA INESTABLE'), break
end

% Funciones de Politica (FP) & IRF

QQ=Q^-1;

q=zeros(root_l_1,variables); % vector donde recogeré los coeficientes para armar la FP-IRF

j=1; % contador

for i=1:variables
    if abs(F(i,i))<1
        q(j,:)=QQ(i,:);
        j=j+1;
    end
end

% vectores donde recogeré la trayectoria de las variables del modelo

z_tray = zeros(50,1);
c_tray = zeros(50,1);
k_tray = zeros(51,1);

% Dinámica del shock de 1% en el AR(1) del modelo

z_tray(1,1) = 1; % Shock en t = 1
k_tray(1,1) = 0; % Cero al ser variable predeterminada (estado) en t = 1.
c_tray(1,1) = -(q(1,2)/q(1,1))*k_tray(1,1) - (q(1,3)/q(1,1))*z_tray(1,1); % Ecuación (20) en t = 1

for i=2:50
    z_tray(i,1) = p_rho*z_tray(i-1,1); % Ecuación (22)
    k_tray(i,1) = -(b1*c_tray(i-1,1)) -(b2*k_tray(i-1,1)) -(b3*z_tray(i-1,1)); % Del sistema log-linealizado: % Ecuación (21)
    c_tray(i,1) = -(q(1,2)/q(1,1))*k_tray(i,1) - (q(1,3)/q(1,1))*z_tray(i,1); % Ecuación (20)
end

% calculamos un periodo adicional para el stock de capital.

k_tray(51,1) = -(b1*c_tray(50,1)) -(b2*k_tray(50,1)) -(b3*z_tray(50,1));

% calculamos el resto de las trayectorias:

% vectores donde recogeré la trayectoria de las variables restantes del modelo

lambda_tray = zeros(50,1);
n_tray = zeros(50,1);
i_tray = zeros(50,1);
y_tray = zeros(50,1);
wp_tray = zeros(50,1);
Rp_tray = zeros(50,1);


for i=1:50
    lambda_tray(i,1) = (-1)*p_sigma*c_tray(i,1); % Ecuacion (23)
    n_tray(i,1) = (-p_sigma/(p_v+p_alfa))*c_tray(i,1) + (1/(p_v+p_alfa))*z_tray(i,1) + (p_alfa/(p_v+p_alfa))*k_tray(i,1); % Ecuacion (25)
    i_tray(i,1) = (1/p_delta)*(k_tray(i+1,1)-k_tray(i,1)*(1-p_delta)); % Ecuacion (27)
    y_tray(i,1) = (c_ss/(c_ss + p_delta*k_ss))*c_tray(i,1) + (p_delta*k_ss/(c_ss + p_delta*k_ss))*i_tray(i,1); % Ecuacion (29)
    wp_tray(i,1) = z_tray(i,1) + p_alfa*k_tray(i,1) - p_alfa*n_tray(i,1); % Ecuacion (31)
    Rp_tray(i,1) = z_tray(i,1) + (p_alfa - 1)*k_tray(i,1) + (1 - p_alfa)*n_tray(i,1); % Ecuacion (33)
end

% Ajusto el lag que realiza dynare en el stock de capital

k_tray = [k_tray(2:51,1)];

% Variables para los gráficos.

time = [0:1:49];
l_cero = zeros(50,1);

% Grafico de las IRF's

% Trayectoria Tecnología: log(z)
subplot(3,3,1), plot(time,z_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(z_tray)+0.05]), title('Tecnología: log(z)')

% Trayectoria Consumo: c
subplot(3,3,2), plot(time,c_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(c_tray)+0.05]), title('Consumo: c')

% Trayectoria Stock de Capital: k
subplot(3,3,3), plot(time,k_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(k_tray)+0.05]), title('Stock de Capital: k')

% Trayectoria Trabajo: n
subplot(3,3,4), plot(time,n_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(n_tray)-0.05 max(n_tray)+0.05]), title('Trabajo: n')

% Trayectoria Inversión: i
subplot(3,3,5), plot(time,i_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(i_tray)-0.05 max(i_tray)+0.05]), title('Inversión: i')

% Trayectoria Producto: y
subplot(3,3,6), plot(time,y_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(y_tray)-0.05 max(y_tray)+0.05]), title('Producto: y')

% Trayectoria Salarios: wp
subplot(3,3,7), plot(time,wp_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 -0.05 max(wp_tray)+0.05]), title('Salarios: wp')

% Trayectoria Costo del Capital: Rp
subplot(3,3,8), plot(time,Rp_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(Rp_tray)-0.05 max(Rp_tray)+0.05]), title('Costo del Capital: Rp')

% Trayectoria Lambda: Lambda
subplot(3,3,9), plot(time,lambda_tray,'b-',time,l_cero,'r-','Linewidth',1.0), axis([0 49 min(lambda_tray)-0.05 0.05]), title('Lambda: Lambda')