% Pb final

% Initialisation
Rt = 6378137;
Hc = 250000;
Rc = Rt + Hc;
Vc = sqrt(3.986e14/Rc);
Vp = Vc;
Delta_V = 0.2*Vc;
thetas = [0.34761; 0.155789; 0.503194; 0.003426];

global theta0 
global theta1

for i=1:3

    Vp = Vp+Delta_V;

    % Création de la fonction d'étagement à Vp fixé
    fc = @(x) pb_etagement(x,Vp);

    % Conditions iniiales et limites
    u_bound=[];
    l_bound=[];
    x0=[100000;50000;10000];
    % Algo utilisé
    bfgs = 1;
    % Paramètres à modifier
    h = 1e-6*x0;
    rho = 1e-1;
    tol = 1e-2;

    % Calcul de me
    me = SQP(fc, u_bound, l_bound, x0, h, rho, tol, bfgs);

    % Conditions initiales de thetas
    thetas_ini = thetas(3:4);
    h2 = 10e-6*thetas_ini;
    rho2 = 1e-2;
    tol2 = 1e-2;

    % Là c'est une fonction de thetas seulement 
    fc2 = @(thetas23) pb_angle(thetas23,me);
    
    thetas = [0;0;0;0];
    ref = 0;
    
    % On fait un balayage
    for theta0=0:5
        for theta1=0:0.2:2
            theta0 = theta0*pi/180
            theta1 = theta1*pi/180
            thetas_res = SQP(fc2, [], [], thetas_ini, h2, rho2, tol2, bfgs);
            if (fc2(thetas_res)<ref)
                ref = fc2(thetas_res);
                thetas = [theta0;theta1;thetas_res];
            end
        end
    end

    Vr = -fc2(thetas_res);
    Delta_V = Vc - Vr;
    fprintf("Delta_V = %f\n", Delta_V)

end 

simulateur_print(me, thetas)