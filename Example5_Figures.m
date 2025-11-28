% Example 7.3 from B. Plestenjak: On properties and numerical computation 
% of critical points of eigencurves of bivariate matrix pencils
%
% The timings and error were obtained by running Example5_Experiment for 
% n = 5:5:45

% Bor Plestenjak 2024, revised in 2025

n = 5:5:45;
t1 =  [3.95E-03 4.49E-02 4.17E-01 2.19E+00 1.07E+01 3.83E+01 1.19E+02 2.83E+02 6.08E+02];
t1r = [8.00E-03 8.13E-02 5.04E-01 2.36E+00 1.11E+01 3.85E+01 1.21E+02 2.86E+02 6.08E+02];
t2 =  [1.77E-02 1.32E-01 4.40E-01 1.45E+00 4.31E+00 1.13E+01 2.74E+01 5.61E+01 1.11E+02];
t2r = [2.28E-02 1.33E-01 4.97E-01 1.56E+00 4.54E+00 1.17E+01 2.81E+01 5.75E+01 1.13E+02];
t3 =  [1.56E-02 7.13E-02 1.72E-01 4.06E-01 9.41E-01 1.74E+00 3.16E+00 5.59E+00 9.82E+00];

e1 =  [3.30E-13 2.99E-12 1.51E-11 4.96E-11 5.47E-10 5.20E-10 1.56E-10 1.56E-10 8.18E-10];
e1r = [7.00E-14 6.36E-14 1.22E-13 1.47E-13 2.95E-13 1.21E-13 1.67E-13 3.73E-13 9.36E-14];
e2 =  [8.48E-15 1.52E-14 5.16E-14 1.12E-13 1.54E-13 3.97E-13 2.24E-12 7.61E-12 4.22E-11];
e2r = [1.06E-15 2.11E-15 3.20E-15 5.31E-15 6.30E-15 1.77E-14 9.33E-15 1.33E-14 2.63E-14];
e3 =  [1.95E-15 3.38E-15 9.08E-15 1.29E-14 1.65E-14 3.42E-14 2.37E-14 3.88E-14 3.40E-14];

PlotSettings

semilogy(n,t1,'r*-','MarkerSize',12)
hold on
semilogy(n,t1r,'r*--','MarkerSize',12)
semilogy(n,t2,'gs-','MarkerSize',12)
semilogy(n,t2r,'gs--','MarkerSize',12)
semilogy(n,t3,'bd--','MarkerSize',12)
legend('Alg. 1','Alg. 1 + Alg. 4','Alg. 3','Alg. 3 + Alg. 4','MFRD + Alg. 4','Location','NorthWest')
ylabel('time in seconds')
xlabel('n')

figure
semilogy(n,e1,'r*-','MarkerSize',12)
hold on
semilogy(n,e1r,'r*--','MarkerSize',12)
semilogy(n,e2,'gs-','MarkerSize',12)
semilogy(n,e2r,'gs--','MarkerSize',12)
semilogy(n,e3,'bd--','MarkerSize',12)
hold off
legend('Alg. 1','Alg. 1 + Alg. 4','Alg. 3','Alg. 3 + Alg. 4','MFRD + Alg. 4','Location','NorthWest')
ylabel('error')
xlabel('n')
axis([0 50 1e-16 1e-6])
