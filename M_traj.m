% Partie 2
me = [17575; 8374; 3663];
mu = 1000;
k = [0.1; 0.15; 0.2];
ms = me.*k;
alpha = [15, 10, 10];
Ve = [2600, 3000, 4400];
% Initial : t=0
t = 0;
Rt = 6378137;
cx = 0.1;
R = [Rt; 0];

% Guess initial de thetas
thetas = pi/180*[5;1.2;10;5];

V = 100*[cos(thetas(1));sin(thetas(1))];

M = mu+sum(ms)+sum(me);

RVM = [R;V;M];

global j;
global theta;
global norm_T;
global q;

% Pour l'affichage de la trajectoire
Rx = [];  % Stockage des positions en x
Ry = [];  % Stockage des positions en y
norm_R = []; % Stockage des normes de R
norm_V = [];     % Normes des vitesses
masses = [];     % Masse en fonction du temps
stock_t = [];
fprintf("Masse totale = %f\n", M);
etage = [0;0;0];
Mf = [0;0;0];
for j=1:3
    % On est à l'étage j
    theta = thetas(j+1);
    norm_T = alpha(j)*RVM(5);
    q = norm_T/Ve(j);
    tej = t + me(j)/q;
    % On calcule le nouveau RVM par intégration ode45
    [t_vals,RVM_] = ode45(@integ,[t tej],RVM);
    
    % Calcul des normes locales
    norm_R_local = vecnorm(RVM_(:, 1:2), 2, 2);
    norm_V_local = vecnorm(RVM_(:, 3:4), 2, 2); % Norme de la vitesse (vx, vy)
    masses_local = RVM_(:, 5);                 % Masse au fil du temps
    norm_R = [norm_R; norm_R_local];
    norm_V = [norm_V; norm_V_local];
    masses = [masses; masses_local];

    etage(j) = size(masses,1);

    % Maj des masses
    RVM = RVM_(end,:)';

    % Mise à jour de la masse après consommation de carburant
    fprintf("M avant séparation de l'étage = %f\n", RVM(5));
    RVM(5) = RVM(5) - ms(j);  % Retirer la masse de l'étage utilisé
    fprintf("M après séparation de l'étage = %f\n", RVM(5));
    Mf(j) = RVM(5);

    % Mise à jour des positions et temps
    stock_t = [stock_t; t_vals];
    Rx = [Rx; RVM_(:, 1)];
    Ry = [Ry; RVM_(:, 2)];

    % Mise à jour du temps
    t = tej;
end

% On veut une vitesse finale de 7.75e3

Rc = 6378137+250000;
fprintf("norm V : %f\n",norm(RVM(3:4),2));
fprintf("norm R : %f\n",norm(RVM(1:2),2)-Rt);
fprintf("V : %f %f\n",RVM(3:4));
fprintf("R : %f %f\n",RVM(1:2));
fprintf("Contraite altitude : %f\n",norm(RVM(1:2),2)-Rc);
fprintf("Contrainte ortho : %f\n",RVM(1:2)'*RVM(3:4));

% Position
subplot(2, 2, 1);
%plot(stock_t, norm_R-Rt, 'b-', 'LineWidth', 2);
plot(stock_t(1:etage(1)), norm_R(1:etage(1))-Rt, 'b-', 'LineWidth', 2);
hold on
plot(stock_t(etage(1)+1:etage(2)), norm_R(etage(1)+1:etage(2))-Rt, 'g-', 'LineWidth', 2);
hold on
plot(stock_t(etage(2)+1:end), norm_R(etage(2)+1:end)-Rt, 'r-', 'LineWidth', 2);

yline(250000, 'r--', 'Objectif');
xlabel('Temps (s)');
ylabel('Norme de R (m)');
title("Évolution de la position |R| en fonction du temps");
grid on;

% Vitesse
subplot(2, 2, 2);
%plot(stock_t, norm_V, 'r-', 'LineWidth', 2);
plot(stock_t(1:etage(1)), norm_V(1:etage(1)), 'b-', 'LineWidth', 2);
hold on
plot(stock_t(etage(1)+1:etage(2)), norm_V(etage(1)+1:etage(2)), 'g-', 'LineWidth', 2);
hold on
plot(stock_t(etage(2)+1:end), norm_V(etage(2)+1:end), 'r-', 'LineWidth', 2);

yline(sqrt(3.986e14/Rc), 'r--', 'Objectif');
xlabel('Temps (s)');
ylabel('Norme de V (m/s)');
title("Évolution de la vitesse |V| en fonction du temps");
grid on;

% Masse
subplot(2, 2, 3);
% plot(stock_t, masses, 'g-', 'LineWidth', 2);
plot(stock_t(1:etage(1)+1), masses(1:etage(1)+1), 'b-', 'LineWidth', 2);
hold on
plot([stock_t(etage(1)+1), stock_t(etage(1)+1)],[masses(etage(1)), Mf(1)],'b-', 'LineWidth', 2);
hold on
plot(stock_t(etage(1)+1:etage(2)+1), masses(etage(1)+1:etage(2)+1), 'g-', 'LineWidth', 2);
hold on
plot([stock_t(etage(2)+1), stock_t(etage(2)+1)],[masses(etage(2)), Mf(2)],'g-', 'LineWidth', 2);
hold on
plot(stock_t(etage(2)+1:end), masses(etage(2)+1:end), 'r-', 'LineWidth', 2);
hold on
plot([stock_t(etage(3)), stock_t(etage(3))],[masses(etage(3)), Mf(3)],'r-', 'LineWidth', 2);
xlabel('Temps (s)');
ylabel('Masse (kg)');
title("Évolution de la masse en fonction du temps");
grid on;

%Terre
subplot(2,2,4);
phi = linspace(-pi/2, pi/2, 500);
x_terre = Rt * cos(phi);
y_terre = Rt * sin(phi);

x_orbite = Rc * cos(phi);
y_orbite = Rc * sin(phi);

hold on;
%plot(Rx, Ry, 'k-', 'LineWidth',2);
plot(x_terre, y_terre, 'k-', 'LineWidth',2);
plot(x_orbite, y_orbite, 'k-.', 'LineWidth',2);
plot(Rx(1:etage(1)), Ry(1:etage(1)), 'b-', 'LineWidth',2);
hold on
plot(Rx(etage(1)+1:etage(2)), Ry(etage(1)+1:etage(2)), 'g-', 'LineWidth',2);
hold on
plot(Rx(etage(2)+1:end), Ry(etage(2)+1:end), 'r-', 'LineWidth',2);
xlabel('X');
ylabel('Y');
axis equal;

R = RVM(1:2);
V = RVM(3:4);
norm_R = norm(R,2);
norm_V = norm(V,2);

fprintf("Angle : %f\n", 180*asin(R'*V/(norm_R*norm_V)/pi));

RVM_f = simulateur(me,thetas);
RVM;
