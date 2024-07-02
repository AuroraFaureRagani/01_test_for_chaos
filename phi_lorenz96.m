function phi = phi_lorenz96(Dt, t_final, tau, nosc, F, print)

    NT = fix(t_final/Dt); 

    x = zeros([nosc NT]);
    x(1:nosc,1) = 2*rand([nosc 1])-1;


    k=1;
    j=1;
    phi = zeros([1, ceil(t_final/tau)]);
    for m = 2:NT
        
        if(mod(m,floor(NT/100))==0)
            disp("param = " + print + "----------" + k + "% ----------");
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
        end %i
        
    csi = 0;
    if(mod(m*Dt, 2.5)==0)
        phi(j) = (1+csi)*(x(2, m) + x(3, m) + x(4, m));
        j=j+1;
    end
    end %m
    
    %ok = isempty( find( isnan(x(1,:)) | isinf(x(1,:)), 1 ) );   %Check wheter the integration have worked or not

%     if not(ok) 
%         disp('INTEGRATION DID NOT WORK') 
%     end
end