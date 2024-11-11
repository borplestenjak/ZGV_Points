% Example 7.5 (left hand side Figure 7.6) from B. Plestenjak: On properties 
% and numerical computation of critical points of eigencurves of bivariate 
% matrix pencils

% Bor Plestenjak

n = 100;

A = 5*eye(n)+diag(ones(n-2,1),2)+diag(ones(n-2,1),-2);
C = eye(n);
B = 0.5*eye(n)+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
opts = [];

disp('Computing critical points')
target = -0.15;
neigs = 100;
opts.use_eigs = 1;
opts.target = target;
opts.neigs = 2*neigs;
opts.delta = 1e-3;
rng(1)
[lambda3, mu3] = critical_points_MFRD(A,-B,-C,opts);
if length(lambda3)>neigs
    lambda3 = lambda3(1:neigs);
    mu3 = mu3(1:neigs);
end

PlotSettings

disp('Computing eigencurves')
% we use higher precision so that eigencurves look smooth
MA = mp(A);
MB = mp(B);
MC = mp(C);
n = size(A,1);
if 1==1
    pts = 2000;
    lam = linspace(mp(-0.5),mp(0.5),pts);
    H = [];
    for k = 1:pts
       H(:,k) = double(eig(MA-lam(k)*MB,MC)); 
    end
    pl = kron(ones(n,1),double(lam));
    scatter(pl(:),real(H(:)),5,abs(imag(H(:))),'filled')
    axis([-0.4 0.4 3 7])
    hold on
    plot(real(lambda3),real(mu3),'sb',[target target],[2 8],'--r','MarkerSize',5)
    hold off
    legend('eigencurves','ZGV and 2D points','target')
end

