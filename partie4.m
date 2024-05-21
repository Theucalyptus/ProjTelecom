Rb=3000; % débit binaire
Fe = 6e3;

NBBITS=50000;
Bits=randi([0 1], NBBITS, 1);

EbN0_db=0:1:6; % rapport signal à bruit en Db
EbN0=10.^(EbN0_db/10); % en rapport


%% CONSTANTES FILTES
ROLL_OFF=0.35;
L=10;
M=4; % ordre de la modulation

% Initialisation des constantes du programme
Te = 1 / Fe; % Temps d'échantillonage
Tb=1/Rb; % Temps binaire
Ts=Tb*log2(M); % Temps symbole
Ns=round(Ts/Te); % Facteur de sur-échantillonage

%% MAPPING
% mapping de gray QPSK
% 0 0 -> I=-1, Q=-1
% 0 1 -> I=-1, Q=1
% 1 1 -> I=1, Q=1
% 1 0 -> I=1, Q=-1
Ak_qpsk = 2*Bits(1:2:end) - 1;
Bk_qpsk = 2*Bits(2:2:end) - 1;
Dk_qpsk = Ak_qpsk + 1i*Bk_qpsk;

ind1 = 1:2:NBBITS;
ind2 = 2:2:NBBITS;
SymAsk = 2*Bits(ind1) + Bits(ind2); % symboles 4-ASK
Dk_4ask = pammod(SymAsk, M, 0, "gray");

figure
hold on
plot(Ak_qpsk, Bk_qpsk, "b*")
plot(Dk_4ask, "ro");
title("Constellation en sortie du mapping sans bruit")
xlabel("Ak")
ylabel("Bk")
legend("QPSK", "4-ASK");

%% MODULATEUR BANDE DE BASE
B = rcosdesign(ROLL_OFF, L, Ns, 'sqrt');
u = zeros(1, Ns);
u(1) = 1;
k_qpsk = kron(Dk_qpsk', u);
k_qpsk = [k_qpsk, zeros(1, L/2*Ns)];
h_qpsk=filter(B, 1, k_qpsk); % signal bande de base
h_qpsk=h_qpsk(L/2*Ns+1:end); % suppression des valeurs nulles à cause du retard du filtre
temps = linspace(0, (length(Dk_qpsk)-1)*Ts, length(h_qpsk));  


k_ask = kron(Dk_4ask', u);
k_ask = [k_ask, zeros(1, L/2*Ns)];
h_ask=filter(B, 1, k_ask); % signal bande de base
h_ask=h_ask(L/2*Ns+1:end); % suppression des valeurs nulles à cause du retard du filtre



%% Tracé des DSP
figure
hold on
dsp = pwelch(h_qpsk, [],[],[],Fe,'twosided');
ech_freq=linspace(-Fe/2, Fe/2, length(dsp));
semilogy(ech_freq, fftshift(dsp));
dsp = pwelch(h_ask, [],[],[],Fe,'twosided');
ech_freq=linspace(-Fe/2, Fe/2, length(dsp));
semilogy(ech_freq, fftshift(dsp));
legend('QPSK', 'ASK');
title('DSP')

TEB_qpsk = zeros(length(EbN0), 1); % vecteur des TEB
TEB_qpsk_th = zeros(length(EbN0), 1); % vecteur des TEB
TEB_ask = zeros(length(EbN0), 1); % vecteur des TEB
TEB_ask_th = zeros(length(EbN0), 1); % vecteur des TEB

for j=1:length(EbN0)
    ebn0 = EbN0(j);
    TES_ask_th = 2 * (1 - 1/M) * qfunc(sqrt(6*ebn0*log2(M) / (M*M - 1)));
    TEB_ask_th(j) = TES_ask_th / log2(M);
    TEB_qpsk_th(j) = qfunc(sqrt(4*EbN0(j))*sin(pi/M)); % Es=2*Eb et TEB = TES/log2(M)



    h_bruite_qpsk = bruit_complexe(h_qpsk, Ns, M, ebn0);
    h_bruite_ask = bruit_complexe(h_ask, Ns, M, ebn0);

    %% DEMODULATION QPSK BANDE DE BASE
    h_bruite_qpsk = [h_bruite_qpsk, zeros(1, L/2*Ns)];
    B = rcosdesign(ROLL_OFF, L, Ns, "sqrt");
    Hr = filter(B, 1, h_bruite_qpsk);
    Hr = Hr(L/2*Ns+1:end);

    %% Décision symboles
    seuilR = 0; % seuil sur la partie réelle
    seuilI = 0; % seuil sur la partie imaginaire
    Hr_ech = Hr(1:Ns:end);
    % figure
    % hold on
    % plot(real(Hr_ech), imag(Hr_ech), "r*")
    % title(strcat("Constellation pour Eb/N0=", strcat(num2str(EbN0_db(j)), "db")))
    % xlabel("Ak_qpsk")
    % ylabel("Bk_qpsk")
    DecAk_qpsk = real(Hr_ech) > seuilR;
    DecBk_qpsk = imag(Hr_ech) < seuilI;


    BitsDecodes = zeros(NBBITS, 1);
    %% Dé-mapping
    for i=1:length(DecAk_qpsk)
       if DecAk_qpsk(i)
           BitsDecodes(2*i-1) = 1;
       end
       if DecBk_qpsk(i)
           BitsDecodes(2*i) = 1;
       end
    end
    TEB_qpsk(j) = sum(BitsDecodes ~= Bits) / NBBITS;


    %% DEMODULATION 4-ASK BANDE DE BASE
    h_bruite_ask = [h_bruite_ask, zeros(1, L/2*Ns)];
    B = rcosdesign(ROLL_OFF, L, Ns, "sqrt");
    Hr = filter(B, 1, h_bruite_ask);
    Hr = Hr(L/2*Ns+1:end);

    %% Décision symboles
    NbSym = length(Dk_4ask);
    Hr_ech = Hr(1:Ns:end);
    figure
    hold on
    plot(real(Hr_ech), zeros(length(Hr_ech), 1), "b*")
    title(strcat("Constellation pour Eb/N0=", strcat(num2str(EbN0_db(j)), "db")))
    xlabel("Ak")

    DecSymAsk = pamdemod(Hr_ech, M, 0, "gray");
    BitsDecodes_ask = zeros(2*length(DecSymAsk), 1);
    for i=1:length(DecSymAsk)
        sym = DecSymAsk(i);
        if sym == 1 || sym == 3
            BitsDecodes_ask(2*i) = 1;
        end

        if sym >= 2
            BitsDecodes_ask(2*i-1) = 1;
        end
    end
    TEB_ask(j) = sum(BitsDecodes_ask ~= Bits) / NBBITS;

end

%% Trace du TEB
figure 
hold on
semilogy(EbN0_db, TEB_qpsk, "m-*")
semilogy(EbN0_db, TEB_ask, "r-*")
semilogy(EbN0_db, TEB_ask_th, "b--o")
semilogy(EbN0_db, TEB_qpsk_th, "g--o")
legend("QPSK exp", "4-ASK exp", "4-ASK Théorique", "QPSK Théorique")
title("TEB en fonction du bruit")
xlabel("Eb/N0 (db)")
ylabel("TEB")
yscale('log')
grid("on")
xticks(EbN0_db)