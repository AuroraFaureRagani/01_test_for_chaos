%% Lorenz96 system
clearvars
F = 4.84;                  % forcing parameter
Dt = 0.01;           % timestep for FAST var

t_final = 3e+4;          % duration of the simulation in natural time units (ultima fig generata con 1e+5, tau = 2.5)
NT = fix(t_final/Dt);   % no. of discrete time steps for SLOW var
nosc = 8;               % no. of oscillators [K]

% Initial conditions
x = zeros([nosc NT]);
x(1:nosc,1) = 2*rand([nosc 1])-1;


NT
tau = 0.1;
t_final/tau
%%
%f = @(i, z, t, nosc)(z(mod(i-2, nosc)+1)*(z(mod(i,nosc)+1) - z(mod(i-3, nosc)+1) ) - z(i) + F);
% Euler integration
tic
k=1;
j=1;
phi = zeros([1, ceil(t_final/tau)]);
    for m = 2:NT
        
        if(mod(m,floor(NT/100))==0)
            disp("----------" + k + "% ----------");
            %disp("step " + m + " of " + NT);
            k=k+1;
        end

        F_m = zeros(nosc);
       
        
       for i=1:nosc
            F_m(i) = F_m(i) + (x(mod(i-2, nosc)+1,m-1)*(x(mod(i,nosc)+1,m-1) - x(mod(i-3, nosc)+1,m-1) ) - x(i,m-1) + F) ;
       end
        

        for i=1:nosc
            % one step euler
            x(i,m) = x(i,m-1) + Dt*F_m(i);
            % one step runge kutta 4
            % x(i,m) = rungekutta4(@(z, t) f(i, [x(1:i-1, m-1); z; x(i+1:end, m-1)], t, nosc), (m-2)*Dt, Dt, x(i, m-1), 1);
            %x(i,m) = rungekutta4(@(t, z) f(i, z, t, nosc), (m-2)*Dt, Dt, x(:, m-1), 1);
        end %i
        
    csi = 0;
    if(mod(m*Dt, 2.5)==0)
        phi(j) = (1+csi)*(x(2, m) + x(3, m) + x(4, m));
        j=j+1;
    end
    end %m
    
    ok = isempty( find( isnan(x(1,:)) | isinf(x(1,:)), 1 ) );   %Check wheter the integration have worked or not

    if not(ok) 
        disp('INTEGRATION DID NOT WORK') 
    end
toc
%
% porzione = 1;
% figure()
% plot((1000*(porzione-1)+1:porzione*1000).*Dt, x(1,1000*(porzione-1)+1:porzione*1000))
% title('time serie')
% xlabel('time')
% ylabel('X_1')


cc = 0.9;
[p, q, M, D, K1, K2] = pq_system(numel(phi), phi, cc, tau);

figure()
plot(1:numel(M), M, 1:numel(M), D); grid on
xlim([0 100])
title('mean square displacement versus n')
legend('M_c(n)', 'D_c(n)')
xlabel('log(n)')
ylabel('mean square displacement')
K1
K2
%% plot p vs q
figure()
plot(p(1:10:end),q(1:10:end)); grid on
title('auxiliary system')
xlabel('p')
ylabel('q')
%% plot mean square displacement
figure()
plot(1:numel(M), M, 1:numel(M), D); grid on
xlim([0 100])
title('mean square displacement versus n')
legend('M_c(n)', 'D_c(n)')
xlabel('log(n)')
ylabel('mean square displacement')
K1
K2
% (F=5.5 : K1 = 0.9979, K2 = 0.9454)
% (F=4.5 : K1 = -0.0375, K2 = 0.0325)
%% K vs c (F fixed)
c = linspace(0, 2*pi, 100);
%c = 3*pi/5*rand(1, 100)+pi/5;

F = 5.5;
for i = 1:numel(c)
    c(i)
    phi = phi_lorenz96(Dt, t_final, tau, nosc, F);
    disp('computing K')
    [p, q, M, D, K_corr(i), K_regr(i)] = pq_system(numel(phi), phi, c(i), tau);
end
%%
figure()
plot(c, K_regr, c, K_corr); grid on;
yline(1, '--')
yline(0, '--')
xlim([0, 2*pi])
xlabel('c')
ylabel('K(c)')
title('K_c (F = 4.5)')
legend('regression', 'correlation')

%% K vs F
F = 4.5:0.01:5.5;
c = 3*pi/5*rand(1, 20)+pi/5;
tau = 2.5;
for j = 1:numel(F)
    F(j)
    tic
    phi = phi_lorenz96(Dt, t_final, tau, nosc, F(j));
    disp('computing K')
    for i = 1:numel(c)
        [p, q, M, D, K_corr(i), K_regr(i)] = pq_system(numel(phi), phi, c(i), tau);
    end

    Kcorr(j) =  median(K_corr);
    Kregr(j) = median(K_regr);
    disp('done')
    toc
end
%%
figure()
plot(F, Kcorr, F, Kregr); grid on, hold on
title('Classifier K')
legend('correlation method', 'regression method')
xlabel('F')
ylabel('K')

%% K vs tau

F = 4.9;
c = 3*pi/5*rand(1, 1)+pi/5;
tau = 0.1:0.1:2.5;
for j = 1:numel(tau)
    tau(j)
    
    phi = phi_lorenz96(Dt, t_final, tau(j), nosc, F, tau(j));
    disp('computing K')
    for i = 1:numel(c)
        [p, q, M, D, K_corr(i), K_regr(i)] = pq_system(numel(phi), phi, c(i), tau);
    end

    Kcorr(j) =  median(K_corr)
    Kregr(j) = median(K_regr);
    disp('done')
end
%%
figure()
plot(tau, Kcorr, tau, Kregr); grid on, hold on
title('Classifier K')
legend('correlation method', 'regression method')
xlabel('\tau_s')
ylabel('K')

%% K vs N

F = 4.5:0.01:5.5;
c = 3*pi/5*rand(1, 5)+pi/5;
tau = 2.5;
t_final = 10.^[2:5];



for k = 1:numel(t_final)
    t_final(k)
    
    for j = 1:numel(F)
        F(j)
        
        phi = phi_lorenz96(Dt, t_final(k), tau, nosc, F(j), F(j));
        disp('computing K')
        for i = 1:numel(c)
            [p, q, M, D, K_corr(i), K_regr(i)] = pq_system(numel(phi), phi, c(i), tau);
        end
        Kcorr(j) =  median(K_corr);
        Kregr(j) = median(K_regr);
        disp('done')
        
    end
    figure(1)
    plot(F, Kcorr, 'DisplayName', "N = " + fix(t_final(k)/tau)); grid on, hold on
    figure(2)
    plot(F, Kregr, 'DisplayName', "N = " + fix(t_final(k)/tau)); grid on, hold on
end
%%
figure(1)
title('Classifier K (correlation method)')
xlabel('F')
ylabel('K')
legend()

figure(2)
title('Classifier K (regression method)')
xlabel('F')
ylabel('K')
legend()
