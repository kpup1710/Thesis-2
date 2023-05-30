function [V_loc, alpha_loc, x_sol_loc, y_sol_loc] = ParTask_new(h,f,A,b,g,VArray,hatd,p,m,M)
    format longG;
    alpha_loc = Inf;
    V_loc = [];
    
    for i = 1 : size(VArray,1)
        v = VArray(i,:)';
        [x,t] = PvSolve_new(f,A,b,g,v,hatd);
        w = v + t*hatd;
        y = cellfun(@(c) c(x),f)';
        DrawLineBetweenTwoPoints(v, w);
        plot(y(1),y(2),'go','markersize',3);
        
        alpha_loc = h(x);
        x_sol_loc = x;
        y_sol_loc = y;
        
        V_loc = [V_loc;GetNewProperElements(v,w,p,m,M)];
    end
end