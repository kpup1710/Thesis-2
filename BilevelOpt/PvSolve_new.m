function [ x,t] = PvSolve_new(f,A,b,g,v,hatd)
format long;
options = optimset('Display','off','Algorithm','active-set');
% Pv Solver to find w,t in direction v-y0: w = y0 + t(v-y0)
% We consider the feasible set: Ax <= b && g(x) <= 0
% % where g: non-linear constraint function
% f: X --> Y strictly quasiconvex function
% y0: start point
% v-y0 desired direction
% w \in bound(f(X)) intersect with v-y0 direction

% number of dimensions of space X
n = size(A, 2);
% number of dimensions of space Y
p = length(v);
% \hat{d} ray direction
if isempty(hatd)
   hatd = ones(p, 1); 
end
%d = [0.2;1;1];
%d = [1;1.3];

% xt = (x, t)
function fmint = fmint(xt)
    fmint = xt(n + 1);
end

xt0 = [zeros(n-1,1);1;1];

function [c,ceq] = nonlconPv(xt)
    x = xt(1:n);
    t = xt(n+1);
    
    % Condition: f(x) - t(v-y0) - y0 <= 0
    % Meaning: f(x) is an inner point of Y
    
    c = [];
    for i = 1 : p
       %c = [c, f{i}(x) - t*(v(i)-y0(i))-y0(i)];
       c = [c, f{i}(x) - v(i)-t*hatd(i)];
    end
    %c = [c, cellfun(@(c) c(x), f) - t*(v-y0)-y0];
    
    % Condition: g(x) <= 0   
    c = [c, g(x)];
    % for Problem 8
    % c = [c, g{1}(x), g{2}(x)];
    
    % Condition: t >= 1
    %c = [c, 1 - t];
    c = [c, t];
    
    % Condition: Ax <= b
    % Will be added by constraint: [A'; zeros(1, size(A, 1))]',b - a, b in
    % fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
    % b: column vector
    ceq = [];
end

[xt_sol, t_sol, exitflag] = fmincon(@fmint, xt0, [A'; zeros(1, size(A, 1))]',b,[],[],[],[], @nonlconPv, options);

if(exitflag == -2)
     disp('No solution found. (When solve Pv)');
     return;
end

x = xt_sol;
x(n+1) = [];

t = t_sol;
end