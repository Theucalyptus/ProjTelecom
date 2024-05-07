function [mapping, h_bdb, h_p, Ns, M] = dvbs_mod(Bits, Rb,Fe, Fp)
% Réalise un modulateur DVBS
%   Entrées : Bits: information binaire à transmettre
%             Rb: débit binaire voulue
%             Fe : fréquence d'échantillonage
%             Fp : fréquence porteuse
%   Sortie :    mapping = [I Q] le mapping de la QPSK
%               H_e le signal émis en sortie du modulateur

    %% CONSTANTES
    ROLL_OFF=0.35;
    L=10;
    M=4; % ordre de la modulation

    % Initialisation des constantes du programme
    Te = 1 / Fe; % Temps d'échantillonage
    Tb=1/Rb; % Temps binaire
    Ts=Tb/log2(M); % Temps symbole
    Ns=round(Ts/Te); % Facteur de sur-échantillonage
   
    %% MAPPING
    % mapping de gray
    % 0 0 -> I=0, Q=0
    % 0 1 -> I=0, Q=1
    % 1 1 -> I=1, Q=1
    % 1 0 -> I=1, Q=0
    Ak = Bits(1:2:end);
    Bk = Bits(2:2:end);
    Dk = Ak + 1i*Bk;
    mapping = Dk;

    %% MODULATEUR BANDE DE BASE
    B = rcosdesign(ROLL_OFF, L, Ns);
    u = zeros(1, Ns);
    u(1) = 1;
    k = kron(Dk', u);
    length(k)
    h_bdb=filter(B, 1, k); % signal bande de base
    h_bdb=h_bdb(L/2:end); % suppression des valeurs nulles à cause du retard du filtre

    %% PASSAGE SUR PORTEUSE
    temps = linspace(0, (length(Dk)-1)*Ts, length(h_bdb));  
    h=h_bdb.*exp(2i*pi*Fp*temps); % signal transporté sur porteuse

end

