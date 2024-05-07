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

    %% Démodulation
    temps = linspace(0, Te, length(H));
    I = H_b.*2*cos(2*pi*Fp*temps);
    Q = H_b.*2*cos(2*pi*Fp*temps);

    L=10;
    ROLL_OFF=0.35;
    B = rcosdesign(ROLL_OFF, L, Ns);
    I_filtre = filter(B, 1, I);
    Q_filtre = filter(B, 1, Q);
    Hr = I_filtre + 1i*Q_filtre;
    Hr = Hr(L/2:end);

    %% Décision symboles
    seuilR = 0; % seuil sur la partie réelle
    seuilI = 0; % seuil sur la partie imaginaire
    NbSym = length(H_b) / Ns;
    N0 = Ns; % instant d'échantillonage
    Hr_ech = Hr(N0 + [0:NbSym-1]*Ns);
    DecAk = real(Hr_ech) > seuilR;
    DecBk = imag(Hr_ech) > seuilI;


    %% Dé-mapping
    BitsDecode = ;



end

