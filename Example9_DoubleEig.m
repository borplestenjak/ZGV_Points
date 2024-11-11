n = 25;
rng(3);
A = randn(n);
B = randn(n);
C = eye(n);

tic; [mu1,lambda1] = critical_points_MFRD(A,C,B); ta=toc

target = 0;
neigs = 6;

opts = [];
opts.use_eigs = 1;
opts.target = target;
opts.neigs = 2*neigs;
tic; [mu2,lambda2] = critical_points_project(A,C,B,opts); tb=toc

plot(real(lambda1),imag(lambda1),'.b','MarkerSize',20)
hold on
plot(0,0,'dg')
plot(real(lambda2),imag(lambda2),'or','MarkerSize',9)
hold off
axis([-0.5 0.5 -0.5 0.5])
legend('all','target','KS')

lambda3 = [];
mu3 = [];

for k=1:length(lambda2)
   [mm,ll,x,y,flag,err,step] = ZGV_GaussNewton(A,C,B,mu2(k),lambda2(k));
   lambda3(k,1) = ll;
   mu3(k,1) = mm;
end

solutions = lambda3