RB=3000; % débit binaire
FP = 3e3; % fréquence porteuse
FE = 24e3; % fréquence d'échantillonnage

NBBITS=2000;
BITS=randi([0 1], NBBITS, 1);

EbN0_db=0:1:8; % rapport signal à bruit en Db
EbN0=10.^(EbN0_db/10); % en rapport


[Mapping, h_bdb, h, Ns, M] = dvbs_mod(BITS, RB, FE, FP); % Modulation

figure
dsp = pwelch(h_bdb, [],[],[],FE,'twosided');
ech_freq=linspace(-FE/2, FE/2, length(dsp));
semilogy(ech_freq, fftshift(dsp));
hold on
dsp = pwelch(h, [],[],[],FE,'twosided');
ech_freq=linspace(-FE/2, FE/2, length(dsp));
semilogy(ech_freq, fftshift(dsp));
legend('Bande de base', 'Transportée sur porteuse');
title('DSP')

TEB = zeros(length(EnN0), 1); % vecteur des TEB
for ebn0=EbN0
    h_bruite = bruit(h, Ns, M, EbN0(4));
    %BitsDecode =  
end


