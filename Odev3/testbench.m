function testbench
close all

numSamples = 161;  % Number of input vectors
Fs = 4000;      % Sampling frequency in Hz
sinFreq = 500; % Input sine frequency in Hz
sinFreq1 = 400;
sinFreq2 = 300;

% Create input data
data = 5 * sin( 2 * pi * (1:numSamples) / (Fs/sinFreq));
temp1 = 5 * sin( 2 * pi * (1:numSamples) / (Fs/sinFreq1));
temp2 = 5 * sin( 2 * pi * (1:numSamples) / (Fs/sinFreq2));
rng('default');
noise = 2*(rand(1,numSamples)-0.5);

indata = data + noise;
indata1 = temp1 + noise;
indata2 = temp2 + noise;

outdata = zeros(1, numSamples);
outdata1 = zeros(1,numSamples);
outdata2 = zeros(1,numSamples);

% Apply filter to each input sample
for n = 1:161
  % Call to design
  outdata(n) = mlhdlc_fir(indata(n));
end
for n= 1:161
      outdata1(n) = mlhdlc_fir(indata1(n));
end
for n= 1:161
      outdata2(n) = mlhdlc_fir(indata2(n));
end      

figure('Name', [mfilename, '_io_plot']);
subplot(2,1,1); plot(data);
axis([1 numSamples -6 6]);
title(['Input = ',num2str(sinFreq),' Hz']);
subplot(2,1,2); plot(noise);
axis([1 numSamples -6 6]);
title('Noise');

% Plot input and output of filter
figure;
subplot(3,2,1); plot(indata);
axis([1 numSamples -6 6]);
title('500 Hz Input');
subplot(3,2,2); plot(outdata);
axis([1 numSamples -6 6]);
title('Filtered Output');

subplot(3,2,3); plot(indata1);
axis([1 numSamples -6 6]);
title('400 Hz Input');
subplot(3,2,4); plot(outdata1);
axis([1 numSamples -6 6]);
title('Filtered Output');

subplot(3,2,5); plot(indata2);
axis([1 numSamples -6 6]);
title('300 Hz Input');
subplot(3,2,6); plot(outdata2);
axis([1 numSamples -6 6]);
title('Filtered Output');

