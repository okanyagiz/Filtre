clc;
clear;

syms z s 

fpass = 300;
fstop = 400;

Fs = 4*10^3;
Ts = 1/Fs;
wpyeni = 2*Fs*tan(pi*fpass/Fs); %% Ön sarma etkisi
wsyeni = 2*Fs*tan(pi*fstop/Fs); %% Ön sarma etkisi

frekans_pass = wpyeni/(2*pi);
frekans_stop = wsyeni/(2*pi);
Zayiflama_pass=3; % db
Zayiflama_stop=60; %

derece=ceil(Zayiflama_stop/(20*log10(frekans_stop/frekans_pass))); %% Filtre derecesinin hesaplanması

for K = 1:(derece-1)/2

    katsayi(K) = cos((2*K+derece-1)*pi/(2*derece));
    polinom(K) = (s^2 -2*s*katsayi(K) + 1);
end

%% Polinomların Değerleri 

r0 = (s+1);     
r1 = (s^2 + 0.1364*s +1);   %% Bunları bu şekilde yazmamım sebebi başka bir kodda bu değerleri kullanmak istememden dolayı.
r2 = (s^2 + 0.4069*s +1);
r3 = (s^2 + 0.6698*s +1);
r4 = (s^2 + 0.9201*s +1);
r5 = (s^2 + 1.1534*s +1);
r6 = (s^2 + 1.3651*s +1);
r7 = (s^2 + 1.5514*s +1);
r8 = (s^2 + 1.7088*s +1);
r9 = (s^2 + 1.8344*s +1);
r10 = (s^2 + 1.9258*s +1);
r11 = (s^2 + 1.9814*s +1);


%% wc frekansında Transfer Fonksiyonu Polinomları

s0 = (s/wpyeni+1);     
s1 = ((s/wpyeni)^2 + 0.1364*s/wpyeni +1);   
s2 = ((s/wpyeni)^2 + 0.4069*s/wpyeni +1);
s3 = ((s/wpyeni)^2 + 0.6698*s/wpyeni +1);
s4 = ((s/wpyeni)^2 + 0.9201*s/wpyeni +1);
s5 = ((s/wpyeni)^2 + 1.1534*s/wpyeni +1);
s6 = ((s/wpyeni)^2 + 1.3651*s/wpyeni +1);
s7 = ((s/wpyeni)^2 + 1.5514*s/wpyeni +1);
s8 = ((s/wpyeni)^2 + 1.7088*s/wpyeni +1);
s9 = ((s/wpyeni)^2 + 1.8344*s/wpyeni +1);
s10 = ((s/wpyeni)^2 + 1.9258*s/wpyeni +1);
s11 = ((s/wpyeni)^2 + 1.9814*s/wpyeni +1);


%% S Domeninde Transfer Fonksiyonu

B(s) = r0*r1*r2*r3*r4*r5*r6*r7*r8*r9*r10*r11;
H_s(s) = 1/B(s);
B1(s) = s0*s1*s2*s3*s4*s5*s6*s7*s8*s9*s10*s11;
c = coeffs(B);
d = coeffs(B1);

[LP_Pay,LP_Payda]=lp2lp(1,double(c),wpyeni);    
[H_jw_surekli_zaman_LP, w_surekli_zaman_LP]=freqs(LP_Pay,LP_Payda,512);         % Ön sarma frekanslarındaki Analog Filtre
subplot(2,1,1)
grid on;
plot(w_surekli_zaman_LP,abs(H_jw_surekli_zaman_LP))
xlabel('w (rad/sn)')
title('Analog Filtre Matlab TOOL')


[Pay_H_z_LP_bilinear,Payda_H_z_LP_bilinear]=bilinear(LP_Pay,LP_Payda,Fs);
[H_jw_ayrik_zaman_LP, w_ayrik_zaman_LP]=freqz(Pay_H_z_LP_bilinear,Payda_H_z_LP_bilinear,512);
subplot(2,1,2)
grid on;
plot(Fs.*w_ayrik_zaman_LP,abs(H_jw_ayrik_zaman_LP));
xlabel('W')
title('Dijital Filtre MATLAB')



%% Z domenine Geçiş // s = (2/T)*(1-z^-1)/(1+z^+1)

z_p = 4.1653*(z-1); %% s'den z'ye bilineer dönüşüm
coe1 = (z+1);       %% wc = wpyeni kesim frekansındaki sistemi elde etmek için dönüşümün aldığı hal => s = (2/T)*(1-z^-1)/((1+z^+1)*wpyeni)
coe2 = (z+1)^2;

z0 = (z+1)^23;
z1 = (z_p) + 1*coe1;
z2 = (z_p^2 + 0.1364*z_p*coe1 + 1*coe2);
z3 = (z_p^2 + 0.4069*z_p*coe1 + 1*coe2);
z4 = (z_p^2 + 0.6698*z_p*coe1 + 1*coe2);
z5 = (z_p^2 + 0.9201*z_p*coe1 + 1*coe2);
z6 = (z_p^2 + 1.1534*z_p*coe1 + 1*coe2);
z7 = (z_p^2 + 1.3651*z_p*coe1 + 1*coe2);
z8 = (z_p^2 + 1.5514*z_p*coe1 + 1*coe2);
z9 = (z_p^2 + 1.7088*z_p*coe1 + 1*coe2);
z10 = (z_p^2 + 1.8344*z_p*coe1 + 1*coe2);
z11 = (z_p^2 + 1.9258*z_p*coe1 + 1*coe2);
z12 = (z_p^2 + 1.9814*z_p*coe1 + 1*coe2);

H_Pay(z) = z0;
H_Payda(z) = z1*z2*z3*z4*z5*z6*z7*z8*z9*z10*z11*z12;

H_z(z) = H_Pay/H_Payda;

c1 = eval(coeffs(H_Pay(z)));
c2 = eval(coeffs(H_Payda(z)));

% c1 = c1/c2(1);
% c2 = c2/c2(1);

[H_jw_ayrik_zaman_LP_bilinear, w_ayrik_zaman_LP_bilinear]=freqz(c1,c2,512);
figure;
grid on;
plot(Fs.*w_ayrik_zaman_LP_bilinear,abs(H_jw_ayrik_zaman_LP_bilinear));
xlabel('W')
title('Transfer Fonksiyonu Aracılığıyla Elde Edilen Dijital Filtre')


% y(n) = -c2(2)*y(n-1) -c2(3)*y(n-2) -c2(4)*y(n-3) -c2(5)*y(n-4) -c2(6)*y(n-5) -c2(7)*y(n-6) -c2(8)*y(n-7) -c2(9)*y(n-8) -c2(10)*y(n-9) ...
%     -c2(11)*y(n-10) -c2(12)*y(n-11) -c2(13)*y(n-12) -c2(14)*y(n-13) -c2(15)*y(n-14) -c2(16)*y(n-15) -c2(17)*y(n-16) -c2(18)*y(n-17) -c2(19)*y(n-18) ...
%     -c2(20)*y(n-19) -c2(21)*y(n-20) -c2(22)*y(n-21) -c2(23)*y(n-22) -c2(24)*y(n-23) + c1(1)*x(n) +c1(2)*x(n-1) +c1(3)*x(n-2) +c1(4)*x(n-3) +c1(5)*x(n-4)...
%     +c1(6)*x(n-5) +c1(7)*x(n-6) +c1(8)*x(n-7) +c1(9)*x(n-8) +c1(10)*x(n-9) +c1(11)*x(n-10) +c1(12)*x(n-11) +c1(13)*x(n-12) +c1(14)*x(n-13) +c1(15)*x(n-14)...
%     +c1(16)*x(n-15) +c1(17)*x(n-16) +c1(18)*x(n-17) +c1(19)*x(n-18) ...
%     +c1(20)*x(n-19) +c1(21)*x(n-20) +c1(22)*x(n-21) +c1(23)*x(n-22) +c1(24)*x(n-23);


n = linspace(0,200,200);
input = sin(2*pi*n/Fs);
y = filter(c1,c2,input);
figure;
plot(n,y)

[R,p,C] = residuez(c1,c2);
for i = 1:derece
    temp(i) = R(i)/(1-p(i)*z^-1);
end

H_p_z(z) = temp(1) + temp(2) + temp(3) + temp(4) + temp(5) + temp(6) + temp(7)+ temp(8)+ temp(9) + temp(10) + temp(11) + temp(12) + temp(13) + temp(14)+ temp(15) ...
    + temp(16) + temp(17) + temp(18) + temp(19) + temp(20) + temp(21) + temp(22) + temp(23) +C;
