function [f, c]=pb_angle(thetas23, me)
    % thetas23 : [theta2; theta3] 

    global theta0
    global theta1
    RVM_f = simulateur(me, thetas23);
    f = -norm(RVM_f(3:4),2);
    Rc = 6378137+250000;
    mu = 3.986e14;
    c = [(norm(RVM_f(1:2),2)-Rc)/Rc; (RVM_f(1:2)'*RVM_f(3:4))/(Rc*sqrt(mu/Rc))];
end