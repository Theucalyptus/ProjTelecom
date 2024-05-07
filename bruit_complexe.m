function [HB] = bruit_complexe(H_bdb, Ns, M, eb_n0)
%% Ajoute un bruit gaussien à un signal H
% Entrées : H le signal à bruiter
%           Ns facteur de sur-échantillonage
%           M Ordre de la modulation
%           eb_n0 Rapport Eb/N0
% Sorties : HB le signal bruité
    
    Pxe = mean(abs(H_bdb).^2); % calcul de la puissance du signal sans bruit
    sigma_carre = Pxe*Ns / (2*log2(M)*eb_n0);
    sigma_n = sqrt(sigma_carre);
    nI = sigma_n*randn(1, length(H_bdb));
    nQ = sigma_n*randn(1, length(H_bdb));    
    bruit = nI + 1i*nQ;
    HB = H_bdb + bruit;

end

