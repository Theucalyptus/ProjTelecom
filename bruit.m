function [HB] = bruit(H, Ns, M, eb_n0)
%% Ajoute un bruit gaussien à un signal H
% Entrées : H le signal à bruiter
%           Ns facteur de sur-échantillonage
%           M Ordre de la modulation
%           eb_n0 Rapport Eb/N0
% Sorties : HB le signal bruité
    
    Px = mean(abs(H).^2); % calcul de la puissance du signal sans bruit
    
    sigma_carre = Px*Ns / (2*log2(M)*eb_n0);
    sigma_n = sqrt(sigma_carre);
    bruit = sigma_n*randn(1, length(H));
    HB = H + bruit;

end

