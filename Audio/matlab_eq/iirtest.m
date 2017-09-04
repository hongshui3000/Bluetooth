cu = zeros(length(y)+1, 3);
for k = 1:length(y)

        a = s2a; b=s2b; g = s2g; 
        si = y(k); 
        cu(k, 1) = (g*si - a(2)*cu(k, 2) - a(3)*cu(k, 3))/a(1);
        yout(k) = b(1)*cu(k, 1) + b(2)*cu(k, 2) + b(3)*cu(k, 3);
        cu(k+1, 3) = cu(k, 2); cu(k+1, 2) = cu(k, 1);
end