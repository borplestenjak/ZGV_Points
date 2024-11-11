% Example 7.5 (right hand side Figure 7.6) from B. Plestenjak: On properties 
% and numerical computation of critical points of eigencurves of bivariate 
% matrix pencils

% Bor Plestenjak

n = 100;

A = 5*eye(n)+diag(ones(n-2,1),2)+diag(ones(n-2,1),-2);
C = eye(n);
B = 0.5*eye(n)+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
opts = [];

disp('Computing critical points')
delta = 1e-3;
target = [0 5];
neigs = 8;
opts.use_jd = 1;
opts.target = target;
opts.neigs = neigs + 5;
opts.delta = delta;
opts_jd = [];
opts_jd.harmonic = 1;
opts_jd.innersteps = 3;
opts_jd.showinfo = 1;
opts_jd.minsize = 10;
opts_jd.maxsize = 15;
opts_jd.target = target;

A2 = A;
B2 = (1+delta)*B;
C2 = C;

opts_jd.M1 = inv(A-target(1)*B-target(2)*C);
opts_jd.M2 = inv(A2-target(1)*B2-target(2)*C2);
opts_jd.maxsteps = 600;
opts.opts_jd = opts_jd;

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
    lam = linspace(mp(-0.25),mp(0.25),pts);
    H = [];
    for k = 1:pts
       H(:,k) = double(eig(MA-lam(k)*MB,MC)); 
    end
    pl = kron(ones(n,1),double(lam));
    scatter(pl(:),real(H(:)),5,abs(imag(H(:))),'filled')
    axis([-0.25 0.25 4.75 5.25])
    hold on
    plot(real(lambda3),real(mu3),'sb',target(1),target(2),'dr','MarkerSize',10)
    hold off
    legend('eigencurves','ZGV and 2D points','target')
end

dif = [lambda3 mu3]-target;
[~,ord] = sort(diag(dif*dif'));
closest4 = [lambda3(ord(1:4)) mu3(ord(1:4))]
