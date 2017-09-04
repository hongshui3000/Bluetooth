clc
clear
close all

fid = fopen('farSpeechEcho.txt','r');
farSpeechEcho = fscanf(fid,'%f');
fclose(fid);

fid = fopen('nearSpeech.txt','r');
nearSpeech = fscanf(fid,'%f');
fclose(fid);

fid = fopen('micSignal.txt','r');
micSignal = fscanf(fid,'%f');
fclose(fid);

fid = fopen('nearSpeech_aec.txt','r');
nearSpeech_aec = fscanf(fid,'%f');
fclose(fid);

fid = fopen('nearSpeech_ns.txt','r');
nearSpeech_ns = fscanf(fid,'%f');
fclose(fid);

fid = fopen('farSpeechEchoEsti.txt','r');
farSpeechEchoEsti = fscanf(fid,'%f');
fclose(fid);



soundsc(nearSpeech_aec,8e3);

smoothFir = ones(1,1024);
ERLE = filter(smoothFir,1,(farSpeechEchoEsti-farSpeechEcho).^2)./filter(smoothFir,1,farSpeechEcho(1:end).^2);
ERLE = -10*log10(ERLE);

figure(1)
subplot(3,1,1)
plot(farSpeechEcho,'b');grid on
subplot(3,1,2)
plot(nearSpeech,'k');grid on
subplot(3,1,3)
plot(micSignal,'r');grid on

figure(2)
subplot(2,1,1)
plot(farSpeechEcho,'b');hold on
plot(farSpeechEchoEsti,'k');hold on
plot(farSpeechEchoEsti-farSpeechEcho,'r');grid on
subplot(2,1,2)
plot(ERLE,'b');grid on

figure(3)
subplot(3,1,1)
plot(nearSpeech,'b');grid on
subplot(3,1,2)
plot(nearSpeech_aec,'k');grid on
subplot(3,1,3)
plot(nearSpeech_ns,'r');grid on







