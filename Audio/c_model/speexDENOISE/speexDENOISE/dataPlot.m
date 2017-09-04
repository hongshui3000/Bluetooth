clc
clear
close all


fid = fopen('nearSpeechDec.txt','r');
nearSpeech = fscanf(fid,'%d');
fclose(fid);

fid = fopen('nearSpeechDenoiseDec.txt','r');
nearSpeechDenoise = fscanf(fid,'%d');
fclose(fid);

fid = fopen('VAD.txt','r');
VAD = fscanf(fid,'%d');
fclose(fid);

figure(2)
subplot(2,1,1)
plot(nearSpeech,'b');grid on
title('before Noise Suppression')
subplot(2,1,2)
plot(nearSpeechDenoise,'k');grid on
title('after Noise Suppression')

figure(3)
plot(nearSpeech/32768,'b');grid on
hold on
plot(kron(VAD,ones(160,1)),'k');grid on
axis([-inf inf -2 2])

soundsc(nearSpeechDenoise,8e3);

% fid = fopen('speechTime.txt','r');
% speechTime = fscanf(fid,'%f');
% fclose(fid);
% 
% fid = fopen('speechFreq.txt','r');
% speechFreq = fscanf(fid,'%f');
% fclose(fid);
% 
% fid = fopen('Window.txt','r');
% window = fscanf(fid,'%f');
% fclose(fid);
% 
% figure(4)
% subplot(2,1,1)
% plot(speechTime,'b');grid on
% subplot(2,1,2)
% plot(speechFreq,'k');grid on
% 
% figure(5)
% plot(window/max(window),'b');grid on
% 
% 
