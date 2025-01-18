function [y]=c(x3,vp)
    ve = [2600; 3000; 4400];
    k = [0.1; 0.15; 0.2];
    gamma = k./(1+k);
    const = ve(3)*(1-gamma(3)*x3);

    y = ve(1)*log((ve(1)-const)/(ve(1)*gamma(1))) + ve(2)*log((ve(2)-const)/(ve(2)*gamma(2))) + ve(3)*log(x3)-vp;
end

function [y]=c_prime(x3)
    ve = [2600; 3000; 4400];
    k = [0.1; 0.15; 0.2];
    gamma = k./(1+k);
    const = ve(3)*(1-gamma(3)*x3);

    y = ve(1)*ve(3)*gamma(3)/(ve(1)-const) + ve(2)*ve(3)*gamma(3)/(ve(2)-const) + ve(3)/x3;
end

ve = [2600; 3000; 4400];
k = [0.1; 0.15; 0.2];
gamma = k./(1+k);
mu = 1000;

Vc = sqrt(3.986e14/6628137);
Vp = Vc+0.2*Vc;
x3 = 4;
for i=1:200
    x3 = x3-c(x3,Vp)/c_prime(x3);
end

const = ve(3)*(1-gamma(3)*x3);
x1 = -(const/ve(1)-1)/gamma(1)
x2 = -(const/ve(2)-1)/gamma(2)
x3

me3 = mu*(x3-1)/(1+k(3)-k(3)*x3)
me2 = (mu+me3*(1+k(3)))*(x2-1)/(1+k(2)-k(2)*x2)
me1 = (mu+me3*(1+k(3))+me2*(1+k(2)))*(x1-1)/(1+k(1)-k(1)*x1)

[me1,me2,me3]