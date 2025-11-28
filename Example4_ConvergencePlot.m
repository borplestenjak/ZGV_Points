% Example 7.2 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils
%
% This example requires Advanpix MCT

% Bor Plestenjak 2024, revised in 2025

% % -----------------------------------------------------------------------
% % next setting increases chances of reproducible results 
% % comment it if speed is more important then reproducibility
maxNumCompThreads(1); % could 
rng(1,'twister')
% % -----------------------------------------------------------------------

clear all

A = [1 2 3 0; 2 0 1 0; 3 1 1 0; 0 0 0 -3];
B = [1 0 1 0; 0 1 1 0; 1 1 0 0; 0 0 0 -3];
C = [2 1 0 0; 1 3 0 0; 0 0 1 0; 0 0 0 1];
n = 4;

mp.Digits(100);
rng(1);
opts = [];
opts.delta = 1e-2;
opts.tol = 1e-90;
opts.showplot = 1;
opts.show = 0;
opts.svmult = opts.delta;

PlotSettings

[lambda, mu, cand_lambda, cand_mu] = critical_points_MFRD(mp(A),mp(B),mp(C),opts);
legend('ZGV points','','','','','','intersections','Location','SouthWest')
xlabel('iteration')
ylabel('error')