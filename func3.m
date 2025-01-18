function [f, c]=func3(x)
    mu = 1700;
    D_V = 11527;
    v = [2647.2, 2922.4, 4344.3];
    k = [0.1101, 0.1532, 0.2154];
    f = mu;
    c = D_V;
    for i=3:-1:1
        f = f + (1+k(i))*x(i); 
        c = c - v(i)*log(f/(f-x(i)));
    end
end