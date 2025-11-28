% Example 7.8 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils

% Bor Plestenjak 2024

% -----------------------------------------------------------------------
% next setting increases chances of reproducible results 
% comment it if speed is more important then reproducibility
maxNumCompThreads(1); % could 
rng(1,'twister')
% -----------------------------------------------------------------------

% We are looking for ZGV points of (lambda^2*L2+lambda*L1+L0+omega^2*M)x=0
n = 3;
L2 = [-1  0.5  0;  0.5 -2  0.5;  0  0.5  -3];
L1 = [ 1 -0.25 0; -0.25 2 -0.25; 0 -0.25 -3];
L0 = [-1  0    0;  0   -2  0;    0  0    -3];
M =  [ 2  1    0;  1    3  1;    0  1     4];

% We can linearize this as a bilinear ZGV problem 
MA = [L0 L1; zeros(n) eye(n)];
MB = [zeros(n) L2; -eye(n) zeros(n)];
MC = [M zeros(n); zeros(n) zeros(n)];
[lambda,mu] = critical_points(MA,MB,MC);

indr = find(abs(imag(lambda))<1e-6);
lambdar = real(lambda(indr));
mur = real(mu(indr));
ZGV_points_omega = [lambdar mur sqrt(mur)]

PlotSettings
pts = 1000;
kl = linspace(-2,2,pts);
H = [];
for k = 1:pts
    H(:,k) = sort(real(sqrt(eig(L0+kl(k)*L1+kl(k)^2*L2,-M)))); 
end
pl = kron(ones(n,1),kl);
scatter(pl(:),real(H(:)),10,abs(imag(H(:))),'filled')
hold on
plot(lambdar,sqrt(mur),'.r','MarkerSize',40)
hold off

axis([-1.5 1.5 0 2])
xlabel('\lambda')
ylabel('\omega')


