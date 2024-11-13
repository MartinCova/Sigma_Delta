%%======================================================================%%
%  ETII 4 - ESE - FEM - CPE Lyon                                         %
%  TP Sigma Delta                                                        % 
%  Etude d'un Modulateur SD d'ordre 2 CIFB type "Booser Wooley"          %
%  Rémy Cellier - remy.cellier@cpe.fr                                    %
%%======================================================================%%

close all;
clear variables;
clc;

%%Initialisation
%===============


% Paramètres du modulateur
OSR         = 64;           % Over Sampling Ratio
fv          = 1e3;          % Frequence du signal d'entree
N           = 2^20;         % Nombre echantillon pour la simu temporelle Simulink
BA          = 20e3;         % Bande audio pour plot fft
Av          = -6;           % Input signal amplitude in dBFs
FS          = 1;            % Pleine échelle (~Vdd)

% Paramètres de la simulation
fs=2*BA*OSR;                            % Frequence d'echantillonnage d'entree
ts=1/fs;                                % Periode d'echantillonnage d'entree
Nin=2^(ceil(log(N/OSR)/log(2))-1);      % Creation du nb de point pour fftin
Tstop=(N-1)*ts;                         % Temps de fin de la simulation (cf parametre du model simulink)
  
% Paramètres d'entrées
Amp=10^(Av/20)*FS;                      % Input sine wave amplitude
fin=double(1/ts/N*2*round(fv*ts*N/2));  % adjust fin for windowing
win=2*pi*fin;                           % fréquence (rad/s) du signal d'entree

% Comportements théoriques CIFB ordre 2
z=tf('z',ts);
NTF=(1-z^-1)^2;
STF=z^-2;



A1 = 10000;
A2 = A1;
GBW = 50e6;
tau = A1/GBW;
C = 100e-12;
a1 = 1/3;
a2 = 3;
R11 = 1/(a1*fs*C);
R22 = 1/(a2*fs*C);
SNR_liste = [];



%% Simulation
%============
open('Model_SD_CIFB_k2_TC_param.mdl');
sim('Model_SD_CIFB_k2_TC_param.mdl');

%% Exploitation des résultats
%============================

% Génération des axes
time = ts * (0:1:N-1);
frequency = (1/Tstop) * (0:1:N-1);

% Calcul des FFT
FFT_in = (fft(in.*blackman(N),N))/(N/2/2); 
FFT_out = (fft(out_SD.*blackman(N),N))/(N/2/2); 
FFT_out_decim = (fft(out_decim.*blackman(N),N))/(N/2/2); 
% /!\ Le niveau de bruit est sur-évalué de 3dB (100% d'erreur) par cette
% méthode car bruit non moyenné sur de multiple FFTs

% Extraction du gabatit de la NTF
[num_NTF,den_NTF]=tfdata(NTF,'v');
[G_NTF,w_NTF]=freqz(num_NTF,den_NTF,N);
f_NTF=w_NTF/(2*pi)*fs;    %Dénormalisation de vecteur pulsation w_NTF borné entre 0 et pi

%% A Faire

% paramètres sigma delta schema bloc



%============================
% Calcul SNR (=> ENOB) et THD dans bande utile sur signal out
% Utiliser Toolbox R. Schreier

%SNR_theo = 6.02*1 + 1.76 -5.17 + 30*log10(OSR); % ordre 1
SNR_theo = 6.02*1 + 1.76 - 12.89 + 50*log10(OSR); % ordre 2

%ENOB_theo = floor(1.5*SNR_theo/9); % ordre 1
ENOB_theo = floor(2.5*SNR_theo/15); % ordre 2

idx_fv = floor(length(frequency)*fv/fs+1);
idx_BA = floor(length(frequency)*BA/(2*fs)+1);

SNR = calculateSNR(FFT_out(1:idx_BA),idx_fv, 3);

%ENOB = floor(1.5*SNR/9); % ordre 1
ENOB = floor(2.5*SNR/15); % ordre 2

harmo_out = FFT_out.*0;
k = 1;

while k*fv < BA
    idx = floor(length(frequency)*k*fv/fs+1);
    harmo_out(idx-3:idx+3) = FFT_out(idx-3:idx+3);
    k = k+1;
end

THD = calculateSNR(harmo_out(1:idx_BA),idx_fv, 3);

%%

% Figures
figure(1);
hold on;
grid on;
axis([0 1/ts -250 0]);
plot(frequency,dbv(FFT_in),'red','LineWidth',3);
plot(frequency,dbv(FFT_out),'green');
plot(f_NTF,dbv(G_NTF)-max(dbv(G_NTF)),'black')
plot([BA BA],[-250 0],'--');
plot([1/(2*ts) 1/(2*ts)],[-250 0],'LineWidth',3);
xlabel('Frequency (Hz)');
ylabel('Modulator output spectrum (dB)');
set(gca,'Xscale', 'log');
legend('Entrée \Sigma\Delta','Sortie \Sigma\Delta','Gabarit NTF','Limite Bande Utile','Limite DFT','Location','NorthWest');
hold off;

%Remarque : 26dB de diff entre NTF et OUT car 6n+10log(OSR)=26 
%(6n => quantification, 10 log(OSR) => Noise Shapper)

figure(2); % Afficher SNR, ENOB et THD
hold on;
grid on;
axis([0 1/(ts*OSR*2) -250 0]);
plot(frequency,dbv(FFT_in),'red');
plot(frequency,dbv(FFT_out_decim),'green');
xlabel('Frequency (Hz)');
ylabel('Decimator output spectrum (dB)');
set(gca,'Xscale', 'log');

hold off;

figure(3);
hold on;
grid on;
T = 1:3000;
plot(time(T), in(T));
plot(time(T), out_decim(T))
hold off;