x = sin(2*pi*(0:10000)*1e6/8e6);
y = (sin(2*pi*(0:10000)*10e3/8e6));
z = sign(y);
z(z<0) = 0.01;

xin = x.*z;
xout = filter([1 -1], [1 -(2^11-1)/2^11], xin);

figure(1); 
idx = 0:10000;
plot(idx, x, idx, y, idx, z); grid on;

figure(2); 
plot(idx, xin, idx, xout);
grid on;

