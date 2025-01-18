function [Hk] = modif_H(Hk)
% Entrée
% Hk : hessien du quasi Newton obtenu, il n'est pas forcément définie
% positif
% Sortie
% H : hessien du quasi newton modifié
    
    
    tau = min(eig(Hk));
    Hk = Hk+(abs(tau)+0.0001)*eye(size(Hk,1));
    
    % tau = min(eig(Hk));
    % if tau < 0.1
    %     Hk = Hk+(abs(tau)+0.01)*eye(size(Hk,1));
    % end
   
end