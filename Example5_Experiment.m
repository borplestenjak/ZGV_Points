% Example 7.3 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils

% Performs M runs on a random triple of n x n real matrices and reports
% median computational times and median maximal error of computed
% eigenvalues. Correct eigenvalues are computed in higher precision using
% Advanpix MCT

% Bor Plestenjak 2024, revised in 2025

% -----------------------------------------------------------------------
% next setting increases chances of reproducible results 
% comment it if speed is more important then reproducibility
maxNumCompThreads(1);  
rng(1,'twister')
% -----------------------------------------------------------------------

mp.Digits(34);
generate_sample = 1;
test_algorithm = [1 1 1 1 1];

for n = 5:5:45

fprintf("-------------------------------------\n")
fprintf("n=%d\n",n)
fprintf("-------------------------------------\n")

M = 10; % number of runs

ns = n*(n-1);

t0 = [];
t1 = [];
t2 = [];
t1r = [];
t2r = [];
t3 = [];
err1 = [];
err2 = [];
err1r = [];
err2r = [];
err3 = [];
find0 = [];
find1 = [];
find2 = [];
find1r = [];
find2r = [];
find3 = [];

for k = 1:M
    fprintf(" k=%d\n",k)
    % we construct random matrices for the example
    if generate_sample==1
        rng(k)
        A = randn(n);
        B = randn(n);
        C = randn(n);
        
        % For "exact" eigenvalues we use MFRD + Gauss-Newton and refinement in
        % higher precision
        opts0 = [];
        opts0.mprefine = 1;
        opts0.show = 0;
        opts0.delta = 2e-3;
        opts0.opts_mp.tol = 1e6*eps('mp');
        opts0.opts_mp.maxsteps = 5;
        opts0.maxsteps = 10;
        
        tic
        rng(k);
        [lambda,mu] = critical_points_MFRD(A,B,C,opts0);
        t0(k) = toc;
        sol0 = double([lambda mu]);
        find0(k) = length(sol0);
        fprintf('MFRD + MP found %d solutions in %8.4f secconds\n',find0(k),t0(k))
        if length(sol0)~=ns
            disp('Incorrect number of solution MFRD + MP')
        end
    end

    if test_algorithm(1)==1
        % first method is Algorithm 1 without refinement
        opts1 = [];
        opts1.sc_steps=0;
        opts1.show = 0;
        opts1.GN_refine = 0;
        tic
        rng(k);
        [lambda1,mu1] = critical_points(A,B,C,opts1);
        t1(k) = toc;
        sol1 = [lambda1 mu1];
        find1(k) = length(sol1);
        fprintf('Alg1  found %d solutions in %8.4f secconds\n',find1(k),t1(k))
        if length(sol1)<length(sol0)
            p = matchrows(sol1,sol0);
            err1(k) = max(vecnorm(sol1-sol0(p,:),2,2)./(1+vecnorm(sol0(p,:),2,2)));
            tmp = ones(length(sol0),1);
            tmp(p) = 0;
            miss = find(tmp==1);
            missing_pairs_1 = sol0(miss,:)
        elseif length(sol1)==length(sol0)
            p = matchrows(sol0,sol1);
            err1(k) = max(vecnorm(sol1(p,:)-sol0,2,2)./(1+vecnorm(sol0,2,2)));
        else
            p = matchrows(sol0,sol1);
            err1(k) = max(vecnorm(sol1(p,:)-sol0,2,2)./(1+vecnorm(sol0,2,2)));
            tmp = ones(length(sol1),1);
            tmp(p) = 0;
            miss = find(tmp==1);
            extra_pairs_1 = sol1(miss,:)
        end
    end

    if test_algorithm(2)==1
        % Algorithm 1 with refinement
        opts1r = [];
        opts1r.sc_steps=0;
        opts1r.show = 0;
        opts1r.GN_refine = 1;
        tic
        rng(k);
        [lambda1r,mu1r] = critical_points(A,B,C,opts1r);
        t1r(k) = toc;
        sol1r = [lambda1r mu1r];
        find1r(k) = length(sol1r);
        fprintf('Alg1R found %d solutions in %8.4f secconds\n',find1r(k),t1r(k))
        if length(sol1r)<length(sol0)
            p = matchrows(sol1r,sol0);
            err1r(k) = max(vecnorm(sol1r-sol0(p,:),2,2)./(1+vecnorm(sol0(p,:),2,2)));
            tmp = ones(length(sol0),1);
            tmp(p) = 0;
            miss = find(tmp==1);
            missing_pairs_1r = sol0(miss,:)
        elseif length(sol1r)==length(sol0)
            p = matchrows(sol0,sol1r);
            err1r(k) = max(vecnorm(sol1r(p,:)-sol0,2,2)./(1+vecnorm(sol0,2,2)));
        else
            p = matchrows(sol0,sol1r);
            err1r(k) = max(vecnorm(sol1r(p,:)-sol0,2,2)./(1+vecnorm(sol0,2,2)));
            tmp = ones(length(sol1r),1);
            tmp(p) = 0;
            miss = find(tmp==1);
            extra_pairs_1r = sol1r(miss,:)
        end
    end

    if test_algorithm(3)==1
        % Algorithm 3 without GN refinement
        opts2 = [];
        opts2.show = 0;
        opts2.cplx_rnd = 1;
        opts2.refine = 1;               % option for multipareig
        opts2.GN_refine = 0;
        opts2.sep = 1e-5;
        opts2.YBXtol = 5e-5;
        tic
        rng(k);
        [lambda2,mu2] = critical_points_project(A,B,C,opts2);
        t2(k) = toc;
        sol2 = [lambda2 mu2];
        find2(k) = length(sol2);
        fprintf('Alg3  found %d solutions in %8.4f secconds\n',find2(k),t2(k))
        if length(sol2)<length(sol0)
            p = matchrows(sol2,sol0);
            err2(k) = max(vecnorm(sol2-sol0(p,:),2,2)./(1+vecnorm(sol0(p,:),2,2)));
            tmp = ones(length(sol0),1);
            tmp(p) = 0;
            miss = find(tmp==1);
            missing_pairs_2 = sol0(miss,:)
        elseif length(sol2)==length(sol0)
            p = matchrows(sol0,sol2);
            err2(k) = max(vecnorm(sol2(p,:)-sol0,2,2)./(1+vecnorm(sol0,2,2)));
        else
            p = matchrows(sol0,sol2);
            err2(k) = max(vecnorm(sol2(p,:)-sol0,2,2)./(1+vecnorm(sol0,2,2)));
            tmp = ones(length(sol2),1);
            tmp(p) = 0;
            miss = find(tmp==1);
            extra_pairs_2 = sol2(miss,:)
        end
    end

    if test_algorithm(4)==1
        % Algorithm 3 with refinement
        opts2r = [];
        opts2r.show = 0;
        opts2r.cplx_rnd = 1;
        opts2r.refine = 1;
        opts2r.GN_refine = 1;
        opts2r.sep = 1e-5;
        opts2r.YBXtol = 5e-5;
        tic
        rng(k);
        [lambda2r,mu2r] = critical_points_project(A,B,C,opts2r);
        t2r(k) = toc;
        sol2r = [lambda2r mu2r];
        find2r(k) = length(sol2r);
        fprintf('Alg3R found %d solutions in %8.4f secconds\n',find2r(k),t2r(k))
        if length(sol2r)<length(sol0)
            p = matchrows(sol2r,sol0);
            err2r(k) = max(vecnorm(sol2r-sol0(p,:),2,2)./(1+vecnorm(sol0(p,:),2,2)));
            tmp = ones(length(sol0),1);
            tmp(p) = 0;
            miss = find(tmp==1);
            missing_pairs_2r = sol0(miss,:)
        elseif length(sol2r)==length(sol0)
            p = matchrows(sol0,sol2r);
            err2r(k) = max(vecnorm(sol2r(p,:)-sol0,2,2)./(1+vecnorm(sol0,2,2)));
        else
            p = matchrows(sol0,sol2r);
            err2r(k) = max(vecnorm(sol2r(p,:)-sol0,2,2)./(1+vecnorm(sol0,2,2)));
            tmp = ones(length(sol2r),1);
            tmp(p) = 0;
            miss = find(tmp==1);
            extra_pairs_2r = sol2r(miss,:)
        end
    end

    if test_algorithm(5)==1
        % third option is MFRD + Gauss-Newton
        optsMFRD = [];
        optsMFRD.delta = 2e-3;
        optsMFRD.maxsteps = 10;
        tic
        rng(k);
        [lambda3,mu3] = critical_points_MFRD(A,B,C,optsMFRD);
        t3(k) = toc;
        sol3 = [lambda3 mu3];
        find3(k) = length(sol3);
        fprintf('MFRD  found %d solutions in %8.4f secconds\n',find3(k),t3(k))
        if length(sol3)<length(sol0)
            p = matchrows(sol3,sol0);
            err3(k) = max(vecnorm(sol3-sol0(p,:),2,2)./(1+vecnorm(sol0(p,:),2,2)));
            tmp = ones(length(sol0),1);
            tmp(p) = 0;
            miss = find(tmp==1);
            missing_pairs_3 = sol0(miss,:)
        elseif length(sol3)==length(sol0)
            p = matchrows(sol0,sol3);
            err3(k) = max(vecnorm(sol3(p,:)-sol0,2,2)./(1+vecnorm(sol0,2,2)));
        else
            p = matchrows(sol0,sol3);
            err3(k) = max(vecnorm(sol3(p,:)-sol0,2,2)./(1+vecnorm(sol0,2,2)));
            tmp = ones(length(sol3),1);
            tmp(p) = 0;
            miss = find(tmp==1);
            extra_pairs_3 = sol3(miss,:)
        end
    end

end

t0
if test_algorithm(1)==1
    find1
    err1
    t1
end
if test_algorithm(2)==1
    find1r
    err1r
    t1r
end
if test_algorithm(3)==1
    find2
    err2
    t2
end
if test_algorithm(4)==1
    find2r
    err2r
    t2r
end
if test_algorithm(5)==1
    find3
    err3
    t3
end

median_times = [median(t0) median(t1) median(t1r) median(t2) median(t2r) median(t3)]
median_error = [median(err1) median(err1r) median(err2) median(err2r) median(err3)]
end
