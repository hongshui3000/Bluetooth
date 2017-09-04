
clc
clear 
close all

% Sampling frequency 
Fs = 8e3;


load Audio8KHzOrignal
load nearspeech
load farspeech
num = 256000;
% nearSpeech = 1/2*Audio8KHzOrignal(1:num).';
nearSpeech = v(1:num).';
SNR = 0;
micNoise = mean(abs(nearSpeech))/10^(SNR/20)*randn(1,num);
nearSpeech = nearSpeech+micNoise;

nearSpeechFix = round(nearSpeech*32768);


fid= fopen('nearSpeechDec.txt','w');
for ii = 1:num 
    fprintf(fid,'%d\n',nearSpeechFix(ii));
end
fclose(fid);
    

figure(1)
subplot(2,1,1)
plot(nearSpeech,'b');grid on
subplot(2,1,2)
plot(nearSpeechFix,'k');grid on






