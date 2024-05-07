function [BitsDecode] = dvbs_demod(H_b, Rb, Fe, Fp)
% Réalise un modulateur DVBS
%   Entrées : H_b : signal bruite
%             Rb: débit binaire voulue
%             Fe : fréquence d'échantillonage
%             Fp : fréquence porteuse
%   


    %% Constantes
    Te = 1 / Fe; % Temps d'échantillonage
    Tb=1/Rb; % Temps binaire
    Ts=Tb/log2(M); % Temps symbole
    Ns=round(Ts/Te); % Facteur de sur-échantillonage

    



end

