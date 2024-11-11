% Example 7.5 (right hand side Figure 7.5) from B. Plestenjak: On properties 
% and numerical computation of critical points of eigencurves of bivariate 
% matrix pencils

% Bor Plestenjak

n = 20;

A = 5*eye(n)+1*diag(ones(n-2,1),2)+1*diag(ones(n-2,1),-2);
C = eye(n);
B = 0.5*eye(n)+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
opts = [];

PlotSettings

% figure
pts = 3000;
lam = linspace(-0.5,0.5,pts);
H = [];
for k = 1:pts
   H(:,k) = eig(A-lam(k)*B,C); 
end
pl = kron(ones(n,1),lam);
scatter(pl(:),real(H(:)),10,abs(imag(H(:))),'filled')

% we first compute all critical points using MFRD + Gauss_Newton
opts.gensize = 0;
opts.delta = 1e-4*(1+1i);
opts.goal = '2D';
[lambda0,mu0] = critical_points_MFRD(A,-B,-C,opts); 
sol0 = [lambda0 mu0];

indr = find(abs(imag(lambda0))<1e-6);
lambdar = real(lambda0(indr));
mur = real(mu0(indr));

% we compute only ZGV points using MFRD + Gauss_Newton
opts.gensize = 0;
opts.delta = 1e-4*(1+1i);
opts.goal = 'ZGV';
[lambda1,mu1] = critical_points_MFRD(A,-B,-C,opts); 
sol1 = [lambda1 mu1];

indr1 = find(abs(imag(lambda1))<1e-6);
lambdar1 = real(lambda1(indr1));
mur1 = real(mu1(indr1));

hold on
    plot(lambdar,mur,'.r','MarkerSize',30)
    plot(lambdar1,mur1,'.k','MarkerSize',30)
hold off
axis([-0.25 0.25 5.6 7.2])