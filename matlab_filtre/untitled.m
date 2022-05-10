clc;
clear all;

frekans_pass=300; 
frekans_stop=400; 
Zayiflama_pass=3; % db
Zayiflama_stop=60; %

f_sample=400*10^3; 

N_frekans_sayisi=512;

n=ceil(Zayiflama_stop/(20*log10(frekans_stop/frekans_pass)));

[As_nominal,Bs_nominal]=butter(n,1,'s');  % nominal alçak geçiren filtre

[H_nominal,w_nominal]=freqs(As_nominal,Bs_nominal);

[As,Bs]=butter(n,2*pi*frekans_pass,'s') ;  %kesim frekansı frekans_pass olan filtre

[H,w]=freqs(As,Bs,N_frekans_sayisi);


% çizdirme:
subplot(3,1,1);
plot(w_nominal,abs(H_nominal));
grid on;
xlabel("Frequency")
ylabel("H(e^jw)")
title("Frequency Response of Normalized Filter")
% semilogy(w_nominal,abs(H_nominal))


subplot(3,1,2);
plot(w/(2*pi),abs(H));
grid on;
xlabel("Frequency")
ylabel("H(e^jw)")
title("Frequency Response for 300 Hz Cut-Off Frequency")


subplot(3,1,3);
plot(w,angle(H));
grid on;
xlabel("Frequency")
ylabel("Phase")
title("Phase Response")




