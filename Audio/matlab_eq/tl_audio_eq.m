% Description
%   This is a program for testing cascaded bi-quad filters
%	In this example, three SoS filters are cascaded on an exisiting frequency response
%	to tune the response. However, it is possible to add more stages as needed.
%
%
% Usage
%   <Example usage of the function.>
%
% Change Log
%
%   2/18/16
%   - Initial version
%
% =============================================================================
% <Filename>.m
%   Author: HJ
%
% COPYRIGHT 2016 Telink.
% All rights reserved.  Telink proprietary and confidential.
%
% =============================================================================*/

clear variables;
close all;

fs = 16e3;  %sample frequency
fvec = linspace(0, fs/2, 1000);		%Vector of frequency points for evaluating filter performance

%----------------------------
% 1st stage filter
%----------------------------
if(1)
[b1, a1]= tl_biquad('Lowshelf', 3.9e3, fs, 6, 'S', 1);
mag1 = freqz(b1, a1, fvec, fs);
gd1 = grpdelay(b1, a1, fvec, fs)';

b1i = round(4096*b1);
a1i = round(-1024*a1);
mag1i = freqz(b1i, a1i, fvec, fs);

figure(1);
subplot(2, 1, 1);
semilogx(fvec, 20*log10(abs(mag1)));
grid on;
xlabel('Freq (Hz)'); ylabel('dB');
title('1st Freq response');
subplot(2, 1, 2);
semilogx(fvec, gd1/fs*1e6);
ylabel('us');
xlabel('Freq (Hz)');
grid on;
title('1st Group delay'); 
else
    mag1 = ones(size(fvec));
    gd1 = 0*mag1;
end

%-------------------------------
%	2nd stage filter
%-------------------------------
if(1)
[b2, a2]= tl_biquad('Peak', 1.2e3, fs, 4, 'Q', 2);
mag2 = freqz(b2, a2, fvec, fs);
gd2 = grpdelay(b2, a2, fvec, fs)';

b2i = round(4096*b2);
a2i = round(-1024*a2);
mag2i = freqz(b2i, a2i, fvec, fs);

figure(2);
subplot(2, 1, 1);
semilogx(fvec, 20*log10(abs(mag2)));
grid on;
xlabel('Freq (Hz)'); ylabel('dB');
title('2nd Freq response');
subplot(2, 1, 2);
semilogx(fvec, gd2/fs*1e6);
ylabel('us');
xlabel('Freq (Hz)');
grid on;
title('2nd Group delay');
else
    mag2 = ones(size(fvec));
    gd2 = 0*mag1;
end

%---------------------------
%	3nd stage filter
%---------------------------
if(1)
[b3, a3]= tl_biquad('Peak', 3e3, fs, 4, 'Q', 1);
mag3 = freqz(b3, a3, fvec, fs);
gd3 = grpdelay(b3, a3, fvec, fs)';

b3i = round(4096*b3);
a3i = round(-1024*a3);
mag3i = freqz(b3i, a3i, fvec, fs);

figure(8);
subplot(2, 1, 1);
semilogx(fvec, 20*log10(abs(mag3)));
grid on;
xlabel('Freq (Hz)'); ylabel('dB');
title('3rd Freq response');
subplot(2, 1, 2);
semilogx(fvec, gd3/fs*1e6);
ylabel('us');
xlabel('Freq (Hz)');
grid on;
title('3rd Group delay');
else
    mag3 = ones(size(fvec));
    gd3 = 0*mag3;
end

%-----------------------------------------------
%	Filter response with all stagess cascaded
%-----------------------------------------------
mag = mag1.*mag2.*mag3;
gd = gd1 + gd2 + gd3;
magi = mag1i.*mag2i.*mag3i;

figure(3);
subplot(2, 1, 1);
semilogx(fvec, 20*log10(abs(mag)), fvec, 20*log10(abs(magi))-36);
grid on;
legend('float', 'int');
xlabel('Freq (Hz)'); ylabel('dB');
title('Freq response after all filter');
subplot(2, 1, 2);
semilogx(fvec, gd/fs*1e6);
ylabel('us');
xlabel('Freq (Hz)');
grid on;
title('Group delay after all filter');

%-----------------------------------------------------
% Further process the filter response data to 
% constrain filter data to interested range
%-----------------------------------------------------
idx = (fvec > 150) & (fvec < 7300); 
sfvec = fvec(idx);
dbeq = 20*log10(mag(idx));

%------------------------------------------------------
%	load UEI test data and interpolate the data to the 
%	evaluation frequency points so that the designed 
%	filter can be added to check tuned response
%------------------------------------------------------
load uei.mat;
fuei = freq150;
suei = Spec_data150;
%interpolate to desired points
dbuei = spline(fuei, suei, sfvec);
figure(4);
semilogx(fuei, suei, sfvec, dbuei);
grid on;
xlabel('Freq (Hz)'); ylabel('dB');
legend('uei orig', 'interp');
title('UEI measurements');

%--------------------------------------------
%
% add equalization response to test data
%
%--------------------------------------------
dball = dbeq + dbuei;

figure(5);
semilogx(sfvec, dball);
grid on;
xlabel('Freq (Hz)'); ylabel('dB');
title('Freq response after all processing');

%--------------------------------------------
%
% Export the filter coef for program use
%
%--------------------------------------------
filter1 = [b1i, a1i(2:end)];
filter2 = [b2i, a2i(2:end)];
filter3 = [b3i, a3i(2:end)];

%export coef in ascii file
fid = fopen('filter-tap.txt', 'w');
fprintf(fid, '%d %d %d %d %d\n', filter1, filter2, filter3);
fclose(fid);

%auto-generate Telink tcdb command for writing filter parameters into flash
fid = fopen('filter-script.txt', 'w');
fprintf(fid, '.\\tcdb.exe wf 71000 -s 4k -e\n');
coef = filter1;
fprintf(fid, '.\\tcdb.exe wf 71000 %s %s %s %s %s\n', ...
    dec2hex(typecast(int32(coef(1)), 'uint32'), 8), dec2hex(typecast(int32(coef(2)), 'uint32'), 8), ...
    dec2hex(typecast(int32(coef(3)), 'uint32'), 8), dec2hex(typecast(int32(coef(4)), 'uint32'), 8), ...
    dec2hex(typecast(int32(coef(5)), 'uint32'), 8) );

coef = filter2;
fprintf(fid, '.\\tcdb.exe wf 71014 %s %s %s %s %s\n', ...
    dec2hex(typecast(int32(coef(1)), 'uint32'), 8), dec2hex(typecast(int32(coef(2)), 'uint32'), 8), ...
    dec2hex(typecast(int32(coef(3)), 'uint32'), 8), dec2hex(typecast(int32(coef(4)), 'uint32'), 8), ...
    dec2hex(typecast(int32(coef(5)), 'uint32'), 8) );

coef = filter3;
fprintf(fid, '.\\tcdb.exe wf 71028 %s %s %s %s %s\n', ...
    dec2hex(typecast(int32(coef(1)), 'uint32'), 8), dec2hex(typecast(int32(coef(2)), 'uint32'), 8), ...
    dec2hex(typecast(int32(coef(3)), 'uint32'), 8), dec2hex(typecast(int32(coef(4)), 'uint32'), 8), ...
    dec2hex(typecast(int32(coef(5)), 'uint32'), 8) );

fclose(fid);
