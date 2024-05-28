Rb=3000; % débit binaire
Fe = 6e3;

NBBITS=300000;
Bits=randi([0 1], NBBITS, 1);

EbN0_db=0:1:6; % rapport signal à bruit en Db
EbN0=10.^(EbN0_db/10); % en rapport


%% CONSTANTES FILTRES DVBS
ROLL_OFF=0.35;
L=10;
M=4; % ordre de la modulation


%% CONSTANTES FILTRES DVBS2
ROLL_OFF_s2=0.2;
L=10;
M_s2=8;

% Initialisation des constantes du programme
Te = 1 / Fe; % Temps d'échantillonage
Tb=1/Rb; % Temps binaire

Ts=Tb*log2(M); % Temps symbole QPSK
Ns=round(Ts/Te); % Facteur de sur-échantillonage

Ts_s2=Tb*log2(M_s2); % Temps symbole 8-PSK
Ns_s2=round(Ts_s2/Te);


%% MAPPING
% mapping de gray QPSK
% 0 0 -> I=-1, Q=-1
% 0 1 -> I=-1, Q=1
% 1 1 -> I=1, Q=1
% 1 0 -> I=1, Q=-1
Ak_qpsk = 2*Bits(1:2:end) - 1;
Bk_qpsk = 2*Bits(2:2:end) - 1;
Dk_qpsk = Ak_qpsk + 1i*Bk_qpsk;

% mapping de gray 8-PSK

Dk_s2 = zeros(NBBITS/log2(M_s2), 1);
for i=1:length(Dk_s2)
    j=1 + log2(M_s2)*(i-1);
    SymBits = Bits(j:j+2)';
    if SymBits == [0, 0, 0]
        Dk_s2(i) = exp(1i * pi/8 * 1);
    elseif SymBits == [1, 0, 0]
        Dk_s2(i) = exp(1i * pi/8 * 3);
    elseif SymBits == [1, 1, 0]
        Dk_s2(i) = exp(1i * pi/8 * 5);
    elseif SymBits == [0 1 0]
        Dk_s2(i) = exp(1i * pi/8 * 7);
    elseif SymBits == [0 1 1]
        Dk_s2(i) = exp(1i * pi/8 * 9);
    elseif SymBits == [1 1 1]
        Dk_s2(i) = exp(1i * pi/8 * 11);
    elseif SymBits == [1 0 1]
        Dk_s2(i) = exp(1i * pi/8 * 13);
    elseif SymBits == [0 0 1]
        Dk_s2(i) = exp(1i * pi/8 * 15);
    else
        sprintf("ERREUR");
    end
end

figure
hold on
plot(Ak_qpsk, Bk_qpsk, "c*")
plot(Dk_s2, "rx");
title("Constellation en sortie du mapping sans bruit")
xlabel("Ak")
ylabel("Bk")
legend("DVB-S", "DVB-S2");

%% MODULATEUR BANDE DE BASE
% dvbs
B = rcosdesign(ROLL_OFF, L, Ns, 'sqrt');
u = zeros(1, Ns);
u(1) = 1;
k_qpsk = kron(Dk_qpsk', u);
k_qpsk = [k_qpsk, zeros(1, Ns*L/2)];
h_qpsk=filter(B, 1, k_qpsk); % signal bande de base
h_qpsk=h_qpsk(L/2*Ns+1:end); % suppression des valeurs nulles à cause du retard du filtre

%dvbs2
B_s2 = rcosdesign(ROLL_OFF_s2, L, Ns_s2, 'sqrt');
u_s2 = zeros(1, Ns_s2);
u_s2(1) = 1;
k_s2 = kron(Dk_s2', u_s2);
k_s2 = [k_s2, zeros(1, Ns_s2*L/2)];
h_s2=filter(B_s2, 1, k_s2); % signal bande de base
h_s2=h_s2(Ns_s2*L/2+1:end); % suppression des valeurs nulles à cause du retard du filtre



%% Tracé des DSP
figure
hold on
dsp = pwelch(h_qpsk, [],[],[],Fe,'twosided');
ech_freq=linspace(-Fe/2, Fe/2, length(dsp));
semilogy(ech_freq, fftshift(dsp));
dsp = pwelch(h_s2, [],[],[],Fe,'twosided');
ech_freq=linspace(-Fe/2, Fe/2, length(dsp));
semilogy(ech_freq, fftshift(dsp));
legend('QPSK (DVB-S)', '8-PSK (DVB-S2)');
title('DSP')

TEB_qpsk = zeros(length(EbN0), 1); % vecteur des TEB
TEB_qpsk_th = zeros(length(EbN0), 1); % vecteur des TEB
TEB_8psk = zeros(length(EbN0), 1); % vecteur des TEB
TEB_8psk_th = zeros(length(EbN0), 1); % vecteur des TEB

for j=1:length(EbN0)
    ebn0 = EbN0(j);
    TEB_8psk_th(j) = 2 * qfunc(sqrt(2*3*EbN0(j))*sin(pi/M_s2)) / log2(M_s2); % Es=3*Eb et TEB=TES / log2(M);
    TEB_qpsk_th(j) = qfunc(sqrt(4*EbN0(j))*sin(pi/M)); % Es=2*Eb et TEB = TES/log2(M);



    h_bruite_qpsk = bruit_complexe(h_qpsk, Ns, M, ebn0);
    h_bruite_8psk = bruit_complexe(h_s2, Ns_s2, M_s2, ebn0);
    %h_bruite_8psk = h_s2;

    %% DEMODULATION QPSK BANDE DE BASE
    h_bruite_qpsk = [h_bruite_qpsk, zeros(1, L/2*Ns)];
    B = rcosdesign(ROLL_OFF, L, Ns, "sqrt");
    Hr = filter(B, 1, h_bruite_qpsk);
    Hr = Hr(L/2*Ns+1:end);

    %% Décision symboles
    seuilR = 0; % seuil sur la partie réelle
    seuilI = 0; % seuil sur la partie imaginaire
    Hr_ech = Hr(1:Ns:end);
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


    %% DEMODULATION 8-PSK DVB-S2
    h_bruite_8psk = [h_bruite_8psk, zeros(1, L/2*Ns_s2)];
    B = rcosdesign(ROLL_OFF_s2, L, Ns_s2, "sqrt");
    Hr = filter(B, 1, h_bruite_8psk);
    Hr = Hr(L/2*Ns_s2+1:end);

    %% Décision symboles
    Hr_ech = Hr(1:Ns_s2:end);
    figure
    hold on
    plot(real(Hr_ech), imag(Hr_ech), "g*")
    title(strcat("Constellation pour Eb/N0=", strcat(num2str(EbN0_db(j)), "db")))
    xlabel("Ak")
    ylabel("Bk")

    Dec1 = real(Hr_ech) > 0; % partie réelle positive
    Dec2 = imag(Hr_ech) < 0; % partie imaginaire positive
    Dec3 = abs(real(Hr_ech)) > abs(imag(Hr_ech)); % partie réelle > partie imaginaire
    BitsDecodes_s2 = zeros(log2(M_s2)*length(Hr_ech), 1);
    for i=1:length(Hr_ech)
        w=1 + (i-1)*log2(M_s2);
        if Dec1(i) && Dec2(i) && Dec3(i)
            BitsDecodes_s2(w:w+2) = [0 0 0];
        elseif Dec1(i) && Dec2(i) && ~Dec3(i)
            BitsDecodes_s2(w:w+2) = [1 0 0];
        elseif ~Dec1(i) && Dec2(i) && Dec3(i)
            BitsDecodes_s2(w:w+2) = [0 1 0];
        elseif ~Dec1(i) && Dec2(i) && ~Dec3(i)
            BitsDecodes_s2(w:w+2) = [1 1 0];
        elseif Dec1(i) && ~Dec2(i) && Dec3(i)
            BitsDecodes_s2(w:w+2) = [0 0 1];
        elseif Dec1(i) && ~Dec2(i) && ~Dec3(i)
            BitsDecodes_s2(w:w+2) = [1 0 1];
        elseif ~Dec1(i) && ~Dec2(i) && Dec3(i)
            BitsDecodes_s2(w:w+2) = [0 1 1];
        else
            BitsDecodes_s2(w:w+2) = [1 1 1];
        end
    end

    TEB_8psk(j) = sum(BitsDecodes_s2 ~= Bits) / NBBITS;

end

%% Trace du TEB
figure 
hold on
semilogy(EbN0_db, TEB_qpsk_th, "g--*")
semilogy(EbN0_db, TEB_qpsk, "m-*")
semilogy(EbN0_db, TEB_8psk_th, "b--x")
semilogy(EbN0_db, TEB_8psk, "r-x")
legend("QPSK théorique", "QPSK exp", "8-PSK théorique", "8-PSK exp")
title("TEB en fonction du bruit")
xlabel("Eb/N0 (db)")
ylabel("TEB")
yscale('log')
grid("on")
xticks(EbN0_db)