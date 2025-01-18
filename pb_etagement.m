function [f, c]=pb_etagement(me,Vp)
    % me : [me1, me2, me3]

    mu = 1000;
    v = [2600; 3000; 4400];
    k = [0.1; 0.15; 0.2];
    ms = me.*k;
    M0 = mu+sum(ms)+sum(me);
    c = Vp;
    f = M0;
    Mi = mu;
    for i=3:-1:1
        Mf = Mi+ms(i);
        Mi = Mf+me(i);
        c = c - v(i)*log(Mi/Mf);
    end
end