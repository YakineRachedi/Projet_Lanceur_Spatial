function [f, c]=probleme(fc,x)
    global nfonc;
    nfonc = nfonc+1;
    [f,c] = fc(x);
end