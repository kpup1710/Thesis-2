function [x_sol, y_sol, exitflag] = phi_monotonic(n,p,g,A,b,h,f,y,options)
    format long;
    function [c, ceq] = nonlcon(x)
        c = [g(x)];
        % for Problem 8
        % c = [g{1}(x), g{2}(x)];
        
        for i = 1 : p
            c = [c, f{i}(x) - y(i)];
        end
        ceq = [];
    end
    % case 9
    [x_sol, y_sol, exitflag] = fmincon(h, [ones(1, n-1) 1], A, b', [], [], [], [], @nonlcon, options);
    if(exitflag == -2)
        disp('No solution found. (When solve min Phi(y*))');
        return;
    end
end