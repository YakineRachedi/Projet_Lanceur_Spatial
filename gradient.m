function [gf, gc]=gradient(probleme, func, x0, h)
% Entrée
% x0 : point initial -> vecteur colonne de taille n
% func : fonction qui renvoie f:Rn->R et c:Rn dans Rm
% h : increment -> vecteur colonne de taille n
% Sortie
% gf : gradient de f -> vecteur colonne de taille n
% gc : Jacobienne des contraintes c -> matrice nxm
    
    n = length(x0);
    [f_x0, c_x0] = probleme(func,x0);
    m = length(c_x0);
    gf = zeros(n,1);
    gc = zeros(n,m);

    for i=1:n
        xc = x0;
        xc(i)= xc(i)+h(i);
        [f_xc,c_xc] = probleme(func,xc);
        gf(i)=(f_xc-f_x0)/h(i);
        gc(i,:)=(c_xc-c_x0)/h(i);
    end
end 