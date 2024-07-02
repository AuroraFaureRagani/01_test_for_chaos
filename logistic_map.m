clearvars
N = 2000;
x0 = 0.2;
mu = 3.55;
tau = 1;

% logistic map time series
phi = zeros(1, N);
phi(1) = x0;
for i = 2:1:N
    phi(i) = mu*phi(i-1)*(1-phi(i-1));
end
%%
cc = 0.9;
[p, q, M, D, K1, K2] = pq_system(N, phi, cc, tau);
% plot p vs q
figure()
plot(p,q);
xlabel('p')
ylabel('q')
% plot mean square displacement
figure()
plot(1:numel(M), M, 1:numel(M), D)
xlim([0 200])
title('mean square displacement versus n')
legend('M_c(n)', 'D_c(n)')
xlabel('n')
ylabel('mean square displacement')


%% K vs c
c = linspace(0, 2*pi, 1000);
tau = 0;
%
for i = 1:numel(c)
    
    [p, q, M, D, K_corr(i), K_regr(i)] = pq_system(N, phi, c(i), tau);
%     nn = numel(D);
%     cov = 1/nn*sum(((1:nn) - mean(1:nn)).*(D - mean(D)));
%     aa = ((1:nn) - mean(1:nn));
%     bb = (D - mean(D));
%     K(i) = cov/(1/nn * sqrt(sum(aa.^2)*sum(bb.^2)));
end

figure()
plot(c, K_regr, c, K_corr); grid on;
yline(1, '--')
title('K(c)')
legend('regression', 'correlation')

%% K vs mu
mu = 3.5:0.001:4;
c = 3*pi/5*rand(1, 100)+pi/5;
N = 2000;
x0 = 0.2;
tau = 1;
for j = 1:numel(mu)
    mu(j)
    phi = zeros(1, N);
    phi(1) = x0;
    for n = 2:1:N
        phi(n) = mu(j)*phi(n-1)*(1-phi(n-1));
    end
    for i = 1:numel(c)
        [p, q, M, D, K_corr(i), K_regr(i)] = pq_system(N, phi, c(i), tau);
    end
    Kcorr_N(j) =  median(K_corr);
    Kregr_N(j) = median(K_regr);
end

figure()
plot(mu, Kcorr_N, mu, Kregr_N);
xlabel('\mu')
ylabel('K')

%% K vs N --> non torna

Nvec = 2:1000:2e+4;
c = 3*pi/5*rand(1, 100)+pi/5;
x0 = 0.2;
tau = 1;
%mu = 3.569;
mu = 3.575;
for j = 1:numel(Nvec)
    Nvec(j)
    phi = zeros(1, Nvec(j));
    phi(1) = x0;
    for n = 2:1:Nvec(j)
        phi(n) = mu*phi(n-1)*(1-phi(n-1));
    end
    for i = 1:numel(c)
        [p, q, M, D, K_corr(i), K_regr(i)] = pq_system(Nvec(j), phi, c(i), tau);
    end
    Kcorr_N(j) =  median(K_corr);
    %Kregr_N(j) = median(K_regr);
end
%
figure()
plot(Nvec, Kcorr_N, Nvec, Kregr_N);
legend('corr', 'regr')
xlabel('N')
ylabel('K')