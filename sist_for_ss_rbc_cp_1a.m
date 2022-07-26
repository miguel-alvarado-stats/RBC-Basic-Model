% sist_for_ss_rbc_cp_1a

% Miguel Alvarado

% Escribimos una funcion (sist_for_ss_rbc_cp_1a), para nuestro sistema no 
% lineal deterministico F = 0, que determina el Estado Estacionario (EE) 
% de las variables de interes.
% El sistema viene de las CPO, tomamos las variables en estado estacionario
% esto es: x_(t) = x_(t+1) = x , for all t.
% Esta función acepta por input un vector X, donde X es un vector de variables 
% desconocidas (los EE) que luego son asignadas a una 
% respectiva variable del sistema (no importa el orden) y que genere (la 
% funcion) un vector F(X)como output.

% Además, añadimos como argumento de la función, un otro input P, donde P 
% es un vector de parámetros que se definen abajo de manera general
% (o bien, se definen abajo sus valores numericos).

function F=sist_for_ss_rbc_cp_1a(X, P)

c      = X(1);
k      = X(2);
z      = X(3);

p_beta  = P(1);
p_delta = P(2);
p_alfa  = P(3);
p_sigma = P(4);
p_v     = P(5);

global g h;

g = ((1-p_alfa)/(p_v+p_alfa));
h = ((p_v+1)/(p_v+p_alfa));

% (Ver archivo pdf, sistema 1)

F = [p_alfa*p_beta*((1-p_alfa)^g)*(z^h)*(c^(-p_sigma*g))*(k^(-p_v*g)) + (p_beta*(1-p_delta)) - 1;
     ((1-p_alfa)^g)*(z^k)*(c^(-p_sigma*g))*(k^(p_alfa*h)) - (p_delta*k) - c;
     z - 1];
end
    
    