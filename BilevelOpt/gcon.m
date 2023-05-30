function [c,ceq] = gcon(x)
    format long;
    global g;
% transformed input for fmincon
    c = g(x);
    ceq = [];
end