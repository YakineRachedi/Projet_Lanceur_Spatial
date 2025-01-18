function [lambdk, dk] = pb_quad(Hk, gf, gc, ck)
% Entrée
% Hk : hessien modifié qui est définie positive
% gc : gradient de f au point courant
% gf : gradient de c au point courant
% ck : valeur de c au point courant
% Sortie
% lambdk : multiplicateur de lagrange à l'itération k
% dk : déplacement à l'itération k
    Q_1 = inv(Hk);
    lambdk = -inv(gc'*Q_1*gc)*(gc'*Q_1*gf-ck);
    dk = -Q_1*(gc*lambdk+gf);
end