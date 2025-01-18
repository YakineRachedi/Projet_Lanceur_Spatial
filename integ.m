function dydt = integ(t,RVM)
% RVM : vecteur de taille 5 : R de taille 2
%                             V de taille 2
%                             M de taille 1

global j;
global theta;
global norm_T;
global q;
% Extraction
R = RVM(1:2);
V = RVM(3:4);
M = RVM(5);
norm_R = norm(R,2);
norm_V = norm(V,2);

% Constante
mu = 3.986e14;
cx = 0.1;
Rt = 6378137;
H = 7000;
rho = 1.225*exp(-(norm_R-Rt)/H);

% W
W = -mu*M*R/(norm_R^3);
% D
D = -cx*rho*norm_V*V;
% T
gamma = asin(R'*V/(norm_R*norm_V));
er = R/norm_R;
eh = [-R(2);R(1)]/norm_R;
u = eh*cos(gamma+theta)+er*sin(gamma+theta);
T = norm_T*u;

dydt = [V;(T+W+D)/M;-q];
end