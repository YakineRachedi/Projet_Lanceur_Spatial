function [xk]=SQP(fc, u_bound, l_bound, x0, h, rho, tol, bfgs)

    global nfonc;
    nfonc = 0;

    % Initialisation
    lambdk=0;
    xk = x0;
    xk_1 = x0+ones(size(xk));
    gf_1 = 1;
    gc_1 = 0;
    Hk=0;
    indicateur = 0;
    echec = 0;
    
    k = 0;
    while k < 100000 & norm(xk-xk_1)>tol & rho<10^15
        % Calcul des fonctions et des gradients
        if (k==0 | indicateur~=0)
            [gf, gc] = gradient(@probleme,fc,xk,h);
            [fk, ck] = probleme(fc,xk);
        end
        % Calcul du Hessien
        Hk = hessien(lambdk, xk, xk_1 ,gf, gc, gf_1, gc_1, Hk, bfgs, indicateur);
        % Plus besoin de la remettre a l'identite
        indicateur = 1;
        % Modif du Hessien
        Hk = modif_H(Hk);
        % Résolution du pb quad
        [lambdk, dk] = pb_quad(Hk, gf, gc, ck);
        
        % Vérification du dk trouvé
        [xk_new, indicateur] = globalisation(@probleme,fc, xk, dk, rho, gf, fk, ck, indicateur, u_bound, l_bound);
        
        if indicateur == 0
            echec = echec+1;
            if echec >= 2
                rho = 10*rho;
            end
        else 
            echec = 0;
        end
        
        % Affichage
        xk_print=sprintf('%f ', xk);
        ck_print= sprintf('%f ', ck);
        lambdk_print= sprintf('%f ', lambdk);
        fprintf("Iter = %d | nfonc = %d | xk = (%s) | rho = %d | echec = %d | f(xk) = %f | c(xk) = %s | lambdk = %s | norm(L(xk,lambdk)) = %f\n", k, nfonc, xk_print, rho, echec, fk, ck_print, lambdk_print, norm(gf+gc*lambdk))
        
        if (indicateur~=0)
            % Stockage des valeurs pour l'itération suivante
            xk_1 = xk;
            xk = xk_new;
            gf_1 = gf;
            gc_1 = gc;
        end
        k = k+1;
    end
    xk_print=sprintf('%f ', xk);
    fprintf("xk_final = (%s)\n", xk_print)
end
