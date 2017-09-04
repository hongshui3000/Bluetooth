% Description
%   This is a program for designing biquad digital filters with various characteristics.
%
%
% Usage
%   <Example usage of the function.>
%
% Change Log
%
%   2/18/16
%   - convert to a generic function
%
%   5/3/12
%   - Initial Version created based on Bristow-Johnson's cookbook
%
%
% =============================================================================
% <Filename>.m
%   Author: HJ
%
% COPYRIGHT 2012 - 2016 Telink.
% All rights reserved.  Telink proprietary and confidential.
%
% =============================================================================*/

%----------------------------------------------------------------------------
%
%         WARNING!   WARNING!	WARNING!
%
%	Please do not modify this function directly as it is the basic
%	building blocks for filter tuning. Just need to change parameters
%	in the main program for real tuning
%
%-----------------------------------------------------------------------------


function [bcoef acoef] = tl_biquad(fType, f0, fs, gaindB, option, opval)

%fType = 'Notch';  %'LPF', 'HPF', 'BPF', 'APF', 'Peak', 'Notch', 'Lowshelf', 'Highshelf'
                  %'Notch1'
%f0 = 4000;      %Hz, the significant frequency.
                %It indicates the center freq, corner freq, or shelf
                %midpoint freq depending on which type of filter is being
                %designed
%fs = 16e3;      %Sampling frequency

%gaindB = 1.5;    %used in peaking and shelving filters

%option = 'BW';   %'Q', 'BW', 'S', determines which of the following three 
                % parameters is used in filter design
%Q = 0.707;      %The filter Q value. In peaking filter, A*Q is the actual Q value
%bw = 2;      %The bandwidth in octaves
                %BPF: between -3dB points
                %Notch: between -3dB points
                %Peak: between midpoints (gaindB/2)
%slope = 1;      %shelf slope parameter, The shelf slope in dB/octave remains proportional
                % to "slope" for fixed value of f0/fs and gaindB.

                
%Intermediate values used in parameter calculation
A = 10^(gaindB/40);
omega = 2*pi*f0/fs;
sn = sin(omega);
cs = cos(omega);

%
%    alpha = sin(w0)/(2*Q)                                       (case: Q)
%          = sin(w0)*sinh( ln(2)/2 * BW * w0/sin(w0) )           (case: BW)
%          = sin(w0)/2 * sqrt( (A + 1/A)*(1/S - 1) + 2 )         (case: S)
%
%        FYI: The relationship between bandwidth and Q is
%             1/Q = 2*sinh(ln(2)/2*BW*w0/sin(w0))     (digital filter w BLT)
%        or   1/Q = 2*sinh(ln(2)/2*BW)             (analog filter prototype)
%
%        The relationship between shelf slope and Q is
%             1/Q = sqrt((A + 1/A)*(1/S - 1) + 2)
%

if(strcmp(option, 'Q'))
    Q = opval;
    alpha = sn / (2*Q);
elseif(strcmp(option, 'BW'))
    bw = opval;
    alpha = sn*sinh( log(2)/2*bw*omega/sn );
elseif(strcmp(option, 'S'))
    slope = opval;
    alpha = sn/2*sqrt( (A + 1/A)*(1/slope - 1) + 2);
else
    disp('Unrecognized option');
    return;
end

%Generate parameters for different filters
switch(fType)
    case 'LPF'
        b0 = (1 - cs)/2;
        b1 = 1 - cs;
        b2 = (1 - cs)/2;
        a0 = 1 + alpha;
        a1 = -2 * cs;
        a2 = 1 - alpha;
        
    case 'HPF'
        b0 = (1 + cs) /2;
        b1 = -(1 + cs);
        b2 = (1 + cs) /2;
        a0 = 1 + alpha;
        a1 = -2 * cs;
        a2 = 1 - alpha;
        
    case 'BPF'
        b0 = alpha;
        b1 = 0;
        b2 = -alpha;
        a0 = 1 + alpha;
        a1 = -2 * cs;
        a2 = 1 - alpha;
        
    case 'Notch1'
        %single pole classic notch, tune the a1 value for filter width
        b0 = 1;
        b1 = -1*exp(1j*omega);
        b2 = 0;
        a0 = 1;
        a1 = -0.9*exp(1j*omega);
        a2 = 0;
    case 'Notch'
        %Biquad Notch
        b0 = 1;
        b1 = -2 * cs;
        b2 = 1;
        a0 = 1 + alpha;
        a1 = -2 * cs;
        a2 = 1 - alpha;
        
    case 'APF'  %All pass filter
        b0 = 1 - alpha;
        b1 = -2 * cs;
        b2 = 1 + alpha;
        a0 = 1 + alpha;
        a1 = -2 * cs;
        a2 = 1 - alpha;
        
    case 'Peak'
        b0 = 1 + (alpha * A);
        b1 = -2 * cs;
        b2 = 1 - (alpha * A);
        a0 = 1 + (alpha /A);
        a1 = -2 * cs;
        a2 = 1 - (alpha /A);
        
    case 'Lowshelf'
        b0 = A * ((A + 1) - (A - 1) * cs + 2*sqrt(A)*alpha);
        b1 = 2 * A * ((A - 1) - (A + 1) * cs);
        b2 = A * ((A + 1) - (A - 1) * cs - 2*sqrt(A)*alpha);
        a0 = (A + 1) + (A - 1) * cs + 2*sqrt(A)*alpha;
        a1 = -2 * ((A - 1) + (A + 1) * cs);
        a2 = (A + 1) + (A - 1) * cs - 2*sqrt(A)*alpha;
        
    case 'Highshelf'
        b0 = A * ((A + 1) + (A - 1) * cs + 2*sqrt(A)*alpha);
        b1 = -2 * A * ((A - 1) + (A + 1) * cs);
        b2 = A * ((A + 1) + (A - 1) * cs - 2*sqrt(A)*alpha);
        a0 = (A + 1) - (A - 1) * cs + 2*sqrt(A)*alpha;
        a1 = 2 * ((A - 1) - (A + 1) * cs);
        a2 = (A + 1) - (A - 1) * cs - 2*sqrt(A)*alpha;
end

bcoef = [b0 b1 b2]/a0;
acoef = [a0 a1 a2]/a0;

if(0)
	%Debug loop for showing the filter response and delay information
    b = bcoef; a = acoef;
    %plot the response of the filter
    fvec = linspace(0, fs/2, 1000); 
    mag = freqz(b, a, fvec, fs);
    gd = grpdelay(b, a, fvec, fs);

    fvtool(b, a);
    figure(1);
    subplot(2, 1, 1);
    plot(fvec, 20*log10(abs(mag)));
    grid on;
    xlabel('Freq (Hz)'); ylabel('dB');
    title('Freq response');
    subplot(2, 1, 2);
    plot(fvec, gd/fs*1e6);
    ylabel('us');
    xlabel('Freq (Hz)');
    grid on;
    title('Group delay');

    figure(2);
    subplot(2, 1, 1);
    semilogx(fvec, 20*log10(abs(mag)));
    grid on;
    xlabel('Freq (Hz)'); ylabel('dB');
    title('Freq response');
    subplot(2, 1, 2);
    semilogx(fvec, gd/fs*1e6);
    ylabel('us');
    xlabel('Freq (Hz)');
    grid on;
    title('Group delay');
end