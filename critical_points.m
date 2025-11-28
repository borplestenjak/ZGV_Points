function [lambda, mu] = critical_points(A,B,C,opts)

% [lambda, mu] = critical_points(A,B,C,opts) finds critical points (lambda,mu) 
% for the eigenvalue problem (A + lambda*B + mu*C)x = 0, such that 
% a) 2D points: there exist nonzero vectors x and y such that 
%    (A + lambda*B + mu*C)*x=0, y'*(A + lambda*B + mu*C)=0, y'*B*x=0
% b) ZGV points: 2D points where mu is a simple eigenvalue 
%
% Keep matrices of the problem small, i.e., not larger than 40 x 40.
%
% Input:
%   - A,B,C : square matrices n x n
%   - opts : options
%
% Options in opts:
%   - goal ('2D'): set to 'ZGV' to compute only ZGV points
%   - sc_steps (0): steps of staircase algorithm for an optional initial
%     reduction of the delta matrices (when it works it might speed up the computation)
%   - membtol (1e-6): treshold for preventing to output multiple instances of critical points
%   - multtol (1e-4): treshold for identifying a possible multiple eigenvalue
%   - YBXtol (1e-6): treshold for the condition |y'*B*x| < YBXtol*norm(B) for the 2D point
%   - YCXtol (1e-6): treshold for the condition |y'*C*x| > YCXtol*norm(C) for the ZGV point
%   - GStol  (1e-6): treshold for the condition sigma(n-1) > GStol*sigma(1) for geometric simplicity
%   - GN_refine (0): set to 1 to apply Gauss-Newton refinement of computed points
%   - other options for singgep or staircase_step_cr_np
%   - inviter (1): use inverse iteration for eigenvectors or slow svd (0)
%
% Output:
%   - lambda, mu: coordinates of critical points

% This is Algorithm 1 in B. Plestenjak: On properties and numerical 
% computation of critical points of eigencurves of bivariate matrix pencils

% Bor Plestenjak 2024, revised in 2025

if nargin<4, opts=[]; end

if isfield(opts,'fp_type') && is_numeric_type_supported(opts.fp_type)  
    class_t = opts.fp_type;   
else
    class_t = superiorfloat(A,B,C);
end

if isfield(opts,'goal'),      goal = opts.goal;         else,  goal = '2D';    end
if isfield(opts,'sc_steps'),  sc_steps = opts.sc_steps; else,  sc_steps = 0;   end
if isfield(opts,'membtol'),   membtol = opts.membtol;   else,  membtol = 1e-6; end
if isfield(opts,'YBXtol'),    YBXtol = opts.YBXtol;     else,  YBXtol = 1e-6;  end
if isfield(opts,'YCXtol'),    YCXtol = opts.YCXtol;     else,  YCXtol = 1e-6;  end
if isfield(opts,'GStol'),     GStol = opts.GStol;       else,  GStol = 1e-6;   end
if isfield(opts,'GN_refine'), GN_refine = opts.GN_refine; else,  GN_refine   = 0;   end
if isfield(opts,'multtol'),   multtol = opts.multtol;   else,  multtol = 1e-4; end
if ~isfield(opts,'sep'),      opts.sep = 1e2*sqrt(eps(class_t)); end
if ~isfield(opts,'sepinf'),   opts.sepinf = 1e-1*eps(class_t); end

% Make sure all inputs are of the same numeric type.
if ~isa(A,class_t), A = numeric_t(A,class_t); end
if ~isa(B,class_t), B = numeric_t(B,class_t); end
if ~isa(C,class_t), C = numeric_t(C,class_t); end

lambda = numeric_t([],class_t); 
mu = numeric_t([],class_t); 
n = size(A,2);
Z = zeros(n,class_t);

% linearization for the second equation
A2 = [A Z;  B A];
B2 = [B Z;  Z B];
C2 = [C Z;  Z C];

[Delta0,Delta1,Delta2] = twopar_delta(A,-B,-C,A2,-B2,-C2);

% optional initial staircase reduction of singular pencils
if sc_steps>0
    Delta = {Delta0,Delta1,Delta2};
    for k=1:sc_steps
        Delta = staircase_step_cr_np(Delta,opts);
    end
    Delta0 = Delta{1};  
    Delta1 = Delta{2};  
    Delta2 = Delta{3};
end

% We solve a singular two-parameter eigenvalue problem by first applying
% singgep to the pencil (Delta1,Delta0) and then by inserting computed
% lambda's in the GEP and checking the 2D criteria for the computed mu's. 
lambda_pos = singgep(Delta1,Delta0,opts);
m = length(lambda_pos);
normB = norm(B);
normC = norm(C);
solution = [];
for k = 1:m
    M = A + lambda_pos(k)*B;
    [X,D,Y] = eig(M,-C);
    mu_pos = diag(D);
    for j = 1:length(mu_pos)
        if isfinite(mu_pos(j))
            x1 = X(:,j);
            y1 = Y(:,j);
            lambda_cand = lambda_pos(k);
            mu_cand = mu_pos(j);
            condB = abs(y1'*B*x1)/(norm(y1)*norm(x1));
            condC = abs(y1'*C*x1)/(norm(y1)*norm(x1));
            if ~is_in_set(solution,[lambda_cand mu_cand],membtol)
                add_solution = 0;
                if strcmp(goal,'ZGV')
                    % if we are computing ZGV points, we first check if |y'*B*x|=0 and |y'*C*x| is nonzero
                    if (condB<YBXtol*normB) && (condC>YCXtol*normC)
                        % if this is true, we additionally check if mu is geometrically simple
                        s = svd(A + lambda_cand*B + mu_cand*C);
                        if s(n-1) > GStol*s(1)
                            add_solution = 1;
                        end
                    end
                end
                if strcmp(goal,'2D')
                    % if we are computing 2D points, we check that either |y'*B*x|=0 or 
                    % geometric multiplicity is at least two
                    if condB < YBXtol*normB
                        add_solution = 1;
                    else
                        % we estimate multiplicity from the number of close eigenvalues
                        pos = abs(mu_pos-mu_cand)<multtol*(1+abs(mu_cand));
                        alg_mult = sum(pos);
                        if alg_mult>1
                            s = svd(A + lambda_cand*B + mu_cand*C);
                            if s(n-1) <= GStol*s(1)
                                add_solution = 1;
                            end
                        end
                    end
                end
                if add_solution
                    if GN_refine                      
                        % optional refinement with Gauss-Newton method
                        [lambda_candC,mu_candC,~,~,flag,err,step] = ZGV_GaussNewton(A,B,C,lambda_cand,mu_cand,x1,y1);
                        if err(end)<err(1)
                            lambda_cand = lambda_candC;
                            mu_cand = mu_candC;
                            if is_in_set(solution,[lambda_cand mu_cand],membtol)
                                % additional check of the refined point
                                add_solution = 0;
                            end
                        end
                    end
                    if add_solution
                        solution = [solution; [lambda_cand mu_cand]];
                    end
                end
            end
        end
    end
end

if ~isempty(solution)
    lambda = solution(:,1);
    mu = solution(:,2);
else
    lambda = [];
    mu =  [];
end
