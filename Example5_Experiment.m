% Example 7.3 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils

% Performs M runs on a random triple of n x n real matrices and reports
% median computational times and median maximal error of computed
% eigenvalues. Correct eigenvalues are computed in higher precision using
% Advanpix MCT

% Bor Plestenjak 2024

% Change the following two parameters
n = 5; % size of matrices (keep it below 40)
M = 10; % number of runs

ns = n*(n-1);

t0 = [];
t1 = [];
t2 = [];
t3 = [];
err1 = [];
err2 = [];
err3 = [];
find0 = [];
find1 = [];
find2 = [];
find3 = [];
clear scurr scurr1 scurr2 scurr3

for k = 1:M

    % we construct random matrices for the example
    rng(k)
    A = randn(n);
    B = randn(n);
    C = randn(n);
    
    % For "exact" eigenvalues we use MFRD + Gauss-Newton and refinement in
    % higher precision
    opts0 = [];
    opts0.mprefine = 1;
    opts0.show = 0;
    tic
    scurr(k) = rng;
    [lambda,mu] = critical_points_MFRD(A,B,C,opts0);
    t0(k) = toc;
    sol0 = double([lambda mu]);
    find0(k) = length(sol0);
    if length(sol0)~=ns
        disp('Incorrect number of solution MFRD + MP')
    end

    % first method is Algorithm 1
    opts1 = [];
    opts1.sc_steps=0;
    opts1.show = 0;
    tic
    scurr1(k) = rng;
    [lambda1,mu1] = critical_points(A,B,C,opts1);
    t1(k) = toc;
    sol1 = [lambda1 mu1];
    find1(k) = length(sol1);
    if length(sol1)<ns
        p = matchrows(sol1,sol0);
        err1(k) = max(vecnorm(sol1-sol0(p,:),2,2)./vecnorm(sol0(p,:),2,2));
    elseif length(sol1)==ns
        p = matchrows(sol0,sol1);
        err1(k) = max(vecnorm(sol1(p,:)-sol0,2,2)./vecnorm(sol0,2,2));
    else
        p = matchrows(sol0,sol1);
        err1(k) = max(vecnorm(sol1(p,:)-sol0,2,2)./vecnorm(sol0,2,2));
    end

    % second option is Algorithm 2
    opts2 = [];
    tic
    scurr2(k) = rng;
    [lambda2,mu2] = critical_points_project(A,B,C,opts2);
    t2(k) = toc;
    sol2 = [lambda2 mu2];
    find2(k) = length(sol2);
    if length(sol2)<ns
        p = matchrows(sol2,sol0);
        err2(k) = max(vecnorm(sol2-sol0(p,:),2,2)./vecnorm(sol0(p,:),2,2));
    elseif length(sol1)==ns
        p = matchrows(sol0,sol2);
        err2(k) = max(vecnorm(sol2(p,:)-sol0,2,2)./vecnorm(sol0,2,2));
    else
        p = matchrows(sol0,sol2);
        err2(k) = max(vecnorm(sol2(p,:)-sol0,2,2)./vecnorm(sol0,2,2));
    end

    % third option is MFRD + Gauss-Newton
    optsMFRD = [];
    tic
    scurr3(k) = rng;
    [lambda3,mu3] = critical_points_MFRD(A,B,C,optsMFRD);
    t3(k) = toc;
    sol3 = [lambda3 mu3];
    find3(k) = length(sol3);
    if length(sol3)<ns
        p = matchrows(sol3,sol0);
        err3(k) = max(vecnorm(sol3-sol0(p,:),2,2)./vecnorm(sol0(p,:),2,2));
    elseif length(sol3)==ns
        p = matchrows(sol0,sol3);
        err3(k) = max(vecnorm(sol3(p,:)-sol0,2,2)./vecnorm(sol0,2,2));
    else
        p = matchrows(sol0,sol3);
        err3(k) = max(vecnorm(sol3(p,:)-sol0,2,2)./vecnorm(sol0,2,2));
    end

end

find1
find2
find3
t0
t1
err1
t2
err2
t3
err3
median_times = [median(t0) median(t1) median(t2) median(t3)]
median_error = [median(err1) median(err2) median(err3)]