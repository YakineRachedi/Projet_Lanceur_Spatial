function [RVM]=simulateur(me, thetas23)

    global theta0;
    global theta1;
    mu = 1000;
    k = [0.1; 0.15; 0.2];
    ms = me.*k;
    alpha = [15, 10, 10];
    Ve = [2600, 3000, 4400];
    thetas = [theta0; theta1; thetas23];
    % Initial : t=0
    ti = 0;
    Rt = 6378137;
    R = [Rt; 0];
    
    V = 100*[cos(thetas(1));sin(thetas(1))];
    
    M = mu+sum(ms)+sum(me);
    
    RVM = [R;V;M];
    global j;
    global theta;
    global norm_T;
    global q;
    
    for j=1:3
        theta = thetas(j+1);
        norm_T = alpha(j)*RVM(5);
        q = norm_T/Ve(j);
        tf = ti + me(j)/q;
        % On calcule le nouveau RVM par intégration ode45
        [t_vals,RVM_] = ode45(@integ,[ti tf],RVM);
        
        
        % Maj des masses
        RVM = RVM_(end,:)';
        % Mise à jour de la masse après consommation de carburant
        RVM(5) = RVM(5) - ms(j);  % Retirer la masse de l'étage utilisé
   
        % Mise à jour du temps
        ti = tf;
    end
end