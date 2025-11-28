% Table 2 in Example 7.5 from B. Plestenjak: On properties 
% and numerical computation of critical points of eigencurves of bivariate 
% matrix pencils

% Bor Plestenjak 2925

% -----------------------------------------------------------------------
% next setting increases chances of reproducible results 
% comment it if speed is more important then reproducibility
maxNumCompThreads(1); % could 
rng(1,'twister')
% -----------------------------------------------------------------------

n = 10; % size of matrices 

A = 5*eye(n)+1*diag(ones(n-2,1),2)+1*diag(ones(n-2,1),-2);
C = eye(n);
B = 0.5*eye(n)+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
opts = [];

N = 20; % number of runs
goal1 = 64; % target number of all 2D points
goal2 = 39; % target number of ZGV points
 
% we first compute all critical points using MFRD + Gauss_Newton
found_2d = zeros(10,1);
opts.goal = '2D';
opts.refine = 1;
for k = 1:10
    opts.delta = 10^(-k);
    tmp = 0;
    tmp01 = 0;
    fprintf(' delta = 1e-%2d,', k)
    for m = 1:N
        rng(m)
        [lambda0,mu0] = critical_points_MFRD(A,-B,-C,opts); 
        fprintf(" %3d,", length(lambda0))
        tmp = tmp + length(lambda0);
        if length(lambda0)==goal1
            tmp01 = tmp01 + 1;
        end
    end
    found_2d(k) = tmp/(N*goal1);
    fprintf(' in %d runs found %8.3e of all 2D points \n',N,found_2d(k))
end

% then we compute all ZGV points using MFRD + Gauss_Newton
found_ZGV = zeros(10,1);
opts.goal = 'ZGV';
for k = 1:10
    opts.delta = 10^(-k); %*(1+1i);
    tmp = 0;
    tmp01 = 0;
    fprintf(' delta = 1e-%2d,', k)
    for m = 1:N
        rng(m)
        [lambda0,mu0] = critical_points_MFRD(A,-B,-C,opts); 
        fprintf(" %3d,", length(lambda0))
        tmp = tmp + length(lambda0);
        if length(lambda0)==goal2
            tmp01 = tmp01 + 1;
        end
    end
    found_ZGV(k) = tmp/(N*goal2);
    fprintf(' in %d runs found %8.3e of all ZGV points \n',N,found_ZGV(k))
end

% empirical probability of finding all 2D point od type d)
only2D = (found_2d*goal1-found_ZGV*goal2)/(goal1-goal2);

fprintf('   delta  |   all 2D found  | only 2D of type (d) |  ZGV found  | \n')
for k = 1:10
    fprintf(' 1e-%2d,       %8.3e         %8.3e           %8.3e \n',k,found_2d(k),only2D(k),found_ZGV(k))
end
    
    
    