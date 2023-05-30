%%%%%% MAIN PROGRAM %%%%%%
clear all;
format longG;
% Example problems
PROBLEM_SET = 7;
PARALLEL = 0;
% 1. Benson2005 - 2D biobjective nonconvex
% 2. 3D convex problem
PROBLEMS; % Procedure to get problem data

% Some variables
global options paddingbd epsilon hatd A b g f h n p 
% phi = @(y) PHI_PRODUCT(y); % PHI_PRODUCT / PHI_SUM
MAX_NUM_LOOP = 1000;
THREAD_LOAD = 10; % number of P(v) problems solved by a thread

num_reps = 1;
timer = zeros(num_reps,1);
for i = 1 : num_reps
tic
% case 9
[x_sol, y_sol] = Solve_new(n,p,g,A,b,paddingbd,f,options,h,MAX_NUM_LOOP,PARALLEL,THREAD_LOAD,hatd,epsilon);
% [x_sol, y_sol] = Solve_new(n,p,g,g,A,b,paddingbd,f,options,h,MAX_NUM_LOOP,PARALLEL,THREAD_LOAD,hatd,epsilon);
timer(i) = toc;
end
timer(timer>mean(timer)) = [];

disp(x_sol);
fprintf('\nx = [%f, %f]', x_sol);
fprintf('\nh(x) = %f\n', y_sol);
fprintf('[PARALLEL=%d] Average time (%d reps): %4.4f\n', PARALLEL, num_reps, mean(timer));
axis equal;