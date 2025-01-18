function [xk_new, indicateur] = globalisation(probleme, fc, xk, dk, rho, gf, fk, ck, indicateur, u_bound, l_bound)
% Entrée
% fc : fonction qui donne en sortie f(x) et c(x)
% xk : point courant
% dk : déplacement à l'itération k
% rho : coefficient de pénalisation
% gf : gradient de f au point xk
% fk : valeur de f au point xk
% ck : valeur de c au point xk
% Sortie
% xk_new : point pour l'itération suivante
    
    % Cas sans conditions au bord
    if (isempty(u_bound) && isempty(l_bound))
        s = 1; 
        it = 0;
        xk_new = xk+dk;
        % Fixe dans la boucle
        Fk = fk+rho*norm(ck,1);
        dFd = gf'*dk-rho*norm(ck,1);
        % Initialisation
        [fk_s, ck_s]=probleme(fc,xk_new);
        Fk_s = fk_s+rho*norm(ck_s,1);
    
        while (Fk_s > Fk+0.1*s*dFd & it<10)
            s = s/2;
            xk_new = xk+s*dk;
            [fk_s, ck_s]=probleme(fc,xk_new);
            Fk_s = fk_s+rho*norm(ck_s,1);
            it = it+1;
        end
        if it==10
            indicateur = 0;
            xk_new = xk;
        end
    % Cas avec conditions au bord
    else
        s = 1; 
        it = 0;
        xk_new = max(min(xk+dk, u_bound), l_bound);
        % Fixe dans la boucle
        Fk = fk+rho*norm(ck,1);
        dFd = gf'*dk-rho*norm(ck,1);
        % Initialisation
        [fk_s, ck_s]=probleme(fc,xk_new);
        Fk_s = fk_s+rho*norm(ck_s,1);
    
        while (Fk_s > Fk+0.1*s*dFd & it<10)
            s = s/2;
            xk_new = max(min(xk+s*dk, u_bound), l_bound);
            [fk_s, ck_s]=probleme(fc,xk_new);
            Fk_s = fk_s+rho*norm(ck_s,1);
            it = it+1;
        end
        if it==10
            indicateur = 0;
            xk_new = xk;
        end
    end
end