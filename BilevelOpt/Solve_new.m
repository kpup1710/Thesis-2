function [x_sol, h_sol] = Solve_new(n,p,g,A,b,paddingbd,f,options,h,MAX_NUM_LOOP,PARALLEL,THREAD_LOAD,hatd,epsilon)
    %%% STEP 0: INITIALIZATION %%%
    %% Initialize
    format longG;
    x_sol = [];
    h_sol = [];
   
    m = zeros(p, 1);
    %% Find m, M
    function [c, ceq] = gcon(x)
        c = [g(x)];
        % for Problem 8
        % c = [g{1}(x), g{2}(x)];
        ceq = [];
    end

    xInfs = [];
    for i = 1 : p
        % case 9
        [xInf, ~, exitflag] = fmincon(f{i}, [ones(1, n-1) 1], A, b', [], [], [], [], @gcon, options);
        % [xInf, ~, exitflag] = fmincon(f{i}, [zeros(1, n-1) 1], A, b', [], [], [], [], @gcon, options);
        m(i) = f{i}(xInf);
        xInfs = [xInfs; xInf];
    end
    fprintf('m = %s\n', sprintf(' %5.4f ', m));
    
    %% Specify M by maximum value for every min solution
    yInfs = [];
    for i = 1 : p
        yInfs = [yInfs; cellfun(@(c) c(xInfs(i,:)), f)];
    end
%     disp(yInfs);
    M = zeros(p, 1);
    for i = 1 : p
        M(i) = max(yInfs(:,i));
    end
    fprintf('M = %s\n', sprintf(' %5.4f ', M));
%     M = m+paddingbd;
%     M = M + 5;
    fprintf('Box [m, M] = [%s%s]\n',sprintf('%4.2f,',m),sprintf('%4.2f,',M));
    
    %% Specify upper bound
    ub = [];
%     for i = 1 : p
%         ub = [ub; h(xInfs(i,:))];
%     end
    alpha = Inf; % upper bound
    fprintf('Initial upper bound: alpha = min(h(xInf(i))) = %f\n', alpha);
    beta = -Inf; % lower bound
    fprintf('Initial lower bound: beta = h() = %f\n', beta);
    
    %% Specify lower bound
    V = M';
    [x_phi, phiVal, exitflag] = phi_monotonic(n,p,g,A,b,h,f,M,options);
    fprintf('\nFind first lower bound: x* = [%5.3f, %5.3f], phi(M) = %f\n', x_phi, phiVal);
%     disp(x_phi);
%     disp(A*x_phi'-b);
%     disp(g(x_phi));
%     disp(h(x_phi));
%     disp(exitflag);
    
    %% Prepare for parallel
    if PARALLEL
        pool = gcp();
        % manual copies
        P_h = h;
        P_f = f;
        P_A = A;
        P_b = b;
        P_g = g;
        P_hatd = hatd;
        P_p = p;
        P_m = m;
        P_M = M;
        All_alpha_loc = zeros(pool.NumWorkers, 1);
        All_x_sol_loc = zeros(pool.NumWorkers, n);
        All_y_sol_loc = zeros(pool.NumWorkers, p);
        VArray = cell(pool.NumWorkers,1);
    end
    
    %% Main loop of algorithm
    for loop = 1 : MAX_NUM_LOOP
        fprintf('\nLoop #%d\n', loop);
        
        %% Find an efficient solution
        if PARALLEL
            All_V_loc = cell(pool.NumWorkers, 1);
            numWorkingThreads = min(pool.NumWorkers,size(V,1));
            totalLoad = min(size(V,1),THREAD_LOAD*pool.NumWorkers);
            for i = 1:numWorkingThreads
                VArray{i} = V(i:pool.NumWorkers:totalLoad,:);
            end
            parfor i = 1:numWorkingThreads
                [V_loc, alpha_loc, x_sol_loc, y_sol_loc] = ParTask_new(P_h,P_f,P_A,P_b,P_g,VArray{i},P_hatd,P_p,P_m,P_M);
                All_V_loc{i} = V_loc;
                All_alpha_loc(i) = alpha_loc;
                All_x_sol_loc(i,:) = x_sol_loc;
                All_y_sol_loc(i,:) = y_sol_loc;
            end
            All_V_loc = vertcat(All_V_loc{:});
        else
            numWorkingThreads = 1;
            totalLoad = 1;
            [All_V_loc, All_alpha_loc, All_x_sol_loc, All_y_sol_loc] = ParTask_new(h,f,A,b,g,V(1,:),hatd,p,m,M);
        end
        
        disp(All_x_sol_loc);
        disp(V(1,:));
        %% Update upper bound
        V(1:totalLoad,:) = [];
        phiVal(1:totalLoad) = [];
        [min_alpha_loc, i] = min(All_alpha_loc(1:numWorkingThreads));
        if IsSmallerThan(min_alpha_loc, alpha)
            alpha = min_alpha_loc;
            fprintf('UPDATE upper bound: alpha = %f\n', alpha);
%             disp(All_x_sol_loc);
%             disp(All_y_sol_loc);
              disp(All_alpha_loc);
            if PARALLEL
                x_sol = All_x_sol_loc(i,:);
                h_sol = All_alpha_loc(i,:);
            else
                x_sol = All_x_sol_loc;
                h_sol = All_alpha_loc;
            end
        end
        gap = alpha-beta;
        fprintf('Gap: %f\n',gap);
        if IsSmallerThan(gap, epsilon)
            fprintf('\n\nALGO TERMINATES.\nalpha = %f; beta = %f; error = %f < %f\n', alpha, beta, alpha-beta, epsilon);
            return;
        end
        
        % for each newly found point
        for i = 1 : size(All_V_loc, 1)
            new_vk = All_V_loc(i,:)';
            [x_phi, newPhiVal] = phi_monotonic(n,p,g,A,b,h,f,new_vk,options);
            if isempty(V) || IsSmallerThan(newPhiVal,phiVal(1)) % in case of first place left insert
                phiVal = [newPhiVal;phiVal];
                V = [new_vk';V];
            else
                for j = size(V,1):-1:1 % newly found points usually give bigger phi
                    if IsSmallerThan(phiVal(j),newPhiVal)
                        phiVal = [phiVal(1:j);newPhiVal;phiVal(j+1:end)];
                        V = [V(1:j,:);new_vk';V(j+1:end,:)];
                        break;
                    end
                end
            end
        end
        
        %% Update lower bound
        if phiVal(1) > beta 
        % if IsSmallerThan(beta, phiVal(1))
            beta = phiVal(1);
            fprintf('UPDATE lower bound: beta = %f\n', beta);
        end
    end
end