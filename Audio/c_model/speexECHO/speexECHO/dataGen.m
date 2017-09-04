
clc
clear 
close all

% Sampling frequency 
Fs = 8e3;


load Audio8KHzOrignal
load nearspeech
load farspeech
num = 256000;
farSpeech = Audio8KHzOrignal(1:num).'+0.0001*randn(1,num);
% farSpeech = x(1:num).'+0.0001*randn(1,num);
% farSpeech = randn(1,num);
nearSpeech = v(1:num).'+0.0001*randn(1,num);
% nearSpeech = zeros(1,num)+0.0001*randn(1,num);
% nearSpeech = [ zeros(1,num/2) randn(1,num/2)];

% Length of CIR between speaker and mic
ECHO_LEN = round((128e-3)*Fs);
% Addtional echo attenuation in dB
ADAT = -10;
% Additional echo delay in sample points
ADDE = round((0e-3)*Fs);
% Generate CIR between speaker and mic
echoCIR1 = randn(1,ECHO_LEN);
echoCIR2 = randn(1,ECHO_LEN);
% normalize CIR
echoCIR1 = echoCIR1/sqrt(sum(echoCIR1.^2));
echoCIR2 = echoCIR2/sqrt(sum(echoCIR2.^2));
% Apply additional attenuation and delay into CIR
echoCIR1 = [ zeros(1,ADDE) 10^(ADAT/20)*echoCIR1 ];
echoCIR2 = [ zeros(1,ADDE) 10^(ADAT/20)*echoCIR2 ];
% load CIR
% echoCIR = CIR(1:10:end);
% echoCIR = randn(1,96)/8;
% echoCIR = 1;



% farSpeechEcho = [fftfilt(echoCIR1,farSpeech(1:num/2)) fftfilt(echoCIR2,farSpeech(num/2+1:end))];
farSpeechEcho = fftfilt(echoCIR1,farSpeech(1:num));
SNR = 10;
micNoise = mean(abs(farSpeechEcho))/10^(SNR/20)*randn(1,num);
micSignal = nearSpeech+farSpeechEcho+micNoise;



fid= fopen('farSpeech.txt','w');
for ii = 1:num 
    fprintf(fid,'%f\n',farSpeech(ii)*32767);
end
fclose(fid);

fid= fopen('farSpeechEcho.txt','w');
for ii = 1:num 
    fprintf(fid,'%f\n',farSpeechEcho(ii)*32767);
end
fclose(fid);

fid= fopen('nearSpeech.txt','w');
for ii = 1:num 
    fprintf(fid,'%f\n',nearSpeech(ii)*32767);
end
fclose(fid);
    
fid= fopen('micSignal.txt','w');
for ii = 1:num 
    fprintf(fid,'%f\n',micSignal(ii)*32767);
end
fclose(fid);



figure(1)
subplot(2,1,1)
plot(farSpeech,'b');grid on
subplot(2,1,2)
plot(micSignal,'k');grid on





