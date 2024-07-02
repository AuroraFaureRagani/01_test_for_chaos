function [p,q, M, D, K_corr, K_regr] = pq_system(N, phi, c, tau)

N0 = floor(N / 10);


%% one dimension
% pq system driven by phi
nvec = 1:1:N;
p = cumsum(phi(nvec).*cos(nvec.*c));
q = cumsum(phi(nvec).*sin(nvec.*c));

% mean square displacement
M = zeros(1, N0);
nlist = 1:1:N0;
for n = nlist
    %M(n) = tau^2/N* sum( (p(n+1:N) - p(1:N-n)).^2 + (q(n+1:N) - q(1:N-n)).^2 );
    M(n) = 1/N* sum( (p(n+1:N) - p(1:N-n)).^2 + (q(n+1:N) - q(1:N-n)).^2 );
end

% modified mean square displacement 
E_phi = mean(phi);
D = M - (E_phi).^2.*(1-cos(c.*(nlist)))./(1-cos(c));

%% compute K (asymptotic growth rate)
% correlation method
C = cov(D, nlist, 1);
K_corr = C(1,2) / sqrt(C(1,1)*C(2,2));


% regression method : Fit a line in the log-log plot using polyfit
if (min(D)>0)
    fit = polyfit(log(nlist), log(D), 1);
    K_regr = fit(1); % slope
else
    a = 1.1;
    D_tilde = D + a*abs(min(D));
    fit_D = polyfit(log(1:numel(D)), log(D_tilde), 1);
    K_regr = fit_D(1);
end
%%
% p(1, :)=phi(1, :)*cos(c);
% q(1, :)=phi(1, :)*sin(c);
% 
% for n = 2:1:N
%     p(n, :) = p(n-1, :) + phi(n, :).*cos(n*c);
%     q(n, :) = q(n-1, :) + phi(n, :).*sin(n*c);
% end
% 
% % mean square displacement
% N0 = ceil(N/50);
% M = zeros(N0, size(phi(1,:)));
% 
% for n = 1:1:N0
%     M(n, :) = 1/N * sum( (p(n+1:N, :) - p(1:N-n,  :)).^2 + (q(n+1:N, :) - q(1:N-n, :)).^2 );
% end
% 
% % compute K
% % Fit a line using polyfit
% fit = polyfit(log(1:numel(M)), log(M), 1);
% K_M = fit(1); % slope
% 
% % modified mean square displacement 
% E_phi = 1/N*sum(phi);
% D = M - (E_phi).^2*(1-cos(c.*(1:numel(M))))./(1-cos(c));
% 
% a = 1.1;
% D_tilde = D + a*max(abs(D));
% fit_D = polyfit(log(1:numel(D)), log(D_tilde), 1);
% K_D = fit_D(1);

end