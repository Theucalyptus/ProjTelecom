Rb=3000; % débit binaire
Fp = 3e3; % fréquence porteuse
Fe = 24e3; % fréquence d'échantillonnage

NBBITS=2000;
Bits=randi([0 1], NBBITS, 1);

EbN0_db=0:1:16; % rapport signal à bruit en Db
EbN0=10.^(EbN0_db/10); % en rapport


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
% 0 0 -> I=-1, Q=-1
% 0 1 -> I=-1, Q=1
% 1 1 -> I=1, Q=1
% 1 0 -> I=1, Q=-1
Ak = 2*Bits(1:2:end) - 1;
Bk = 2*Bits(2:2:end) - 1;
Dk = Ak + 1i*Bk;
mapping = Dk;

%% MODULATEUR BANDE DE BASE
B = rcosdesign(ROLL_OFF, L, Ns, 'sqrt');
u = zeros(1, Ns);
u(1) = 1;
k = kron(Dk', u);
h_bdb=filter(B, 1, k); % signal bande de base
h_bdb=h_bdb(L/2+1:end); % suppression des valeurs nulles à cause du retard du filtre
h_bdb=[h_bdb, zeros(1, L/2)]; % ajout de zero à la fin
temps = linspace(0, (length(Dk)-1)*Ts, length(h_bdb));  


figure 
hold on
plot(temps, real(h_bdb))
plot(temps, imag(h_bdb))
xlabel("Temps (s)")
ylabel("Signal")
legend("En phase", "En quadrature")
title("Signal transmis en bande de base")


%% PASSAGE SUR PORTEUSE
h_p=h_bdb.*exp(2i*pi*Fp*temps); % signal transporté sur porteuse
h_p=real(h_p);

%% Tracé signal sur porteuse temporel
figure 
plot(h_p)
xlabel("Temps (s)")
ylabel("Signal")
title("Signal transmis sur porteuse")

%% Tracé des DSP
figure
dsp = pwelch(h_bdb, [],[],[],Fe,'twosided');
ech_freq=linspace(-Fe/2, Fe/2, length(dsp));
semilogy(ech_freq, fftshift(dsp));
hold on
dsp = pwelch(h_p, [],[],[],Fe,'twosided');
ech_freq=linspace(Fp-Fe/2, Fp+Fe/2, length(dsp));
semilogy(ech_freq, fftshift(dsp));
legend('Bande de base', 'Transportée sur porteuse');
title('DSP')

TEB = zeros(length(EbN0), 1); % vecteur des TEB
for j=1:length(EbN0)
    ebn0 = EbN0(j);
    h_bruite = bruit(h_p, Ns, M, ebn0);
   

    %% Démodulation
    length(temps)
    length(h_bruite)
    I = 2*h_bruite.*cos(2*pi*Fp*temps);
    Q = 2*h_bruite.*sin(2*pi*Fp*temps);

    L=10;
    ROLL_OFF=0.35;
    B = rcosdesign(ROLL_OFF, L, Ns);
    I_filtre = filter(B, 1, I);
    I_filtre = [I_filtre(L/2+1:end), zeros(1, L/2)];
    Q_filtre = filter(B, 1, Q);
    Q_filtre = [Q_filtre(L/2+1:end), zeros(1, L/2)];
    Hr = I_filtre - 1i*Q_filtre;

    %% Décision symboles
    seuilR = 0; % seuil sur la partie réelle
    seuilI = 0; % seuil sur la partie imaginaire
    NbSym = length(Dk);
    N0 = Ns; % instant d'échantillonage
    Hr_ech = Hr(N0 + [0:NbSym-1]*Ns);
    DecAk = real(Hr_ech) > seuilR;
    DecBk = imag(Hr_ech) < seuilI;


    BitsDecodes = zeros(NBBITS, 1);
    %% Dé-mapping
    for i=1:length(DecAk)
       if DecAk(i)
           BitsDecodes(2*i-1) = 1;
       end
       if DecBk(i)
           BitsDecodes(2*i) = 1;
       end
    end

    TEB(j) = sum(BitsDecodes ~= Bits) / NBBITS;

end

%% Trace du TEB
figure 
hold on
semilogy(EbN0_db, TEB)
legend("QPSK")
title("TEB en fonction du bruit")
xlabel("Eb/N0 (db)")
ylabel("TEB")
yscale('log')
xticks(EbN0_db)