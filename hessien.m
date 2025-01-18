function [H]=hessien(lambdk, xk, xk_1, gf, gc, gf_1, gc_1, Hk_1, bfgs, indicateur)
% Entrée
% lambdk : multiplicateur de lagrange
% xk : point courant
% xk_1 : point précédent
% gf : gradient de f au point xk
% gf_1 : gradient de f au point xk_1
% gc : gradient de c au point xk
% gc_1 : gradient de c au point xk_1
% Hk_1 : hessien au point précédent
% bfgs : 1 si BFGS et 0 si SR1
% indicateur : 0 si on doit mettre la hessienne a l'identite
% Sortie
% H : hessien du quasi Newton
    if indicateur == 1
        dk_1 = xk-xk_1;
        yk_1 = (gf-gf_1)+(gc-gc_1)*lambdk;
        if (bfgs==1)
            cond1 = yk_1'*dk_1;
            if cond1>0
                H = Hk_1+(yk_1*yk_1')/cond1 - (Hk_1*dk_1*(dk_1)'*Hk_1)/(dk_1'*Hk_1*dk_1);
            else 
                H = Hk_1;
            end
        elseif (bfgs==0)
            % Stocker (yk_1-Hk_1*dk_1) ?
            cond2 = dk_1'*(yk_1-Hk_1*dk_1);
            if cond2 ~= 0
                H = Hk_1+(yk_1-Hk_1*dk_1)*(yk_1-Hk_1*dk_1)'/cond2;
            else 
                H = Hk_1;
            end
        else
            error("Erreur : Algo dispo -> BFGS, SR1 ")
        end
    elseif indicateur == 0
        H = eye(length(xk)); % Initialisation avec une matrice identité
    else
        error("Invalid indicator");
    end
end