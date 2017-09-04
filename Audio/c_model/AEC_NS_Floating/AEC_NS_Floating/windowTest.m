clc
clear
close all

N = 128;

x = 0.5-0.5*cos(2*pi*(0:127)/N);

y = zeros(1,N);
for ii = 0:N-1
    tmp = 4*ii/N;
    inv  = 0;
    if(tmp<1)
        tmp1 = tmp;
    elseif(tmp<2)
        tmp1 = 2-tmp;
        inv  =1;
    elseif(tmp<3)
        tmp1 = tmp-2;
        inv = 1;
    else
        tmp1 = 4-tmp;
    end
    
    tmp2 = 1.271903*tmp1;
    tmp3 = (0.5-0.5*cos(0.5*pi*tmp2))^2;
    
    if(inv)
        tmp3 = 1-tmp3;
    end
    
    y(ii+1) = sqrt(tmp3);
end

figure(1)
plot(x,'b');hold on
plot(y,'k');grid on
legend('Hanning','speex window')
            