% % Cas test 1 :
% fc = @func1;
% u_bound=[];
% l_bound=[];
% x0 = [1;1];

% % Cas test 2 :
% fc = @func2;
% u_bound = [0;3;2;0;-1];
% l_bound = [-2;1;0;-3;-3];
% x0=[-1;2;1;-2;-2];
% fc([-1.2366; 2.4616;1.1911;-0.2144;-1.6165])

% % Cas test 3 :
% fc = @func3;
% u_bound=[];
% l_bound=[];
% x0=[100000;50000;10000];
% size(fc([145349; 31215; 7933]))

% Cas Ariane
Vc = sqrt(3.986e14/6628137);
Vp = Vc+0.2*Vc;
fc = @(x) pb_etagement(x,Vp);
u_bound=[1000000;1000000;1000000];
l_bound=[1000;1000;1000];
x0=[100000;50000;10000];

% Algo utilisé
bfgs = 1;
% Paramètres à modifier
h = 1e-6*x0;
rho = 1e-1;
tol = 1e-2;

me = SQP(fc, u_bound, l_bound, x0, h, rho, tol, bfgs);

