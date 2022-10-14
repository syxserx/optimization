function [flux, throughput,YY,C_1_outlet,C_2_outlet,C_1_acm,...
    C_2_acm,flux_prime,C_1_outlet_prime,C_2_outlet_prime, A, T,f]=...
    two_species_evol(a_0,xi,beta,gamma,lambda_1,W)



obj1=W(4);
obj2=1-W(4);

h=0.01; % space grid size
k=0.001; %time grid siz

% lambda_1=1;
lambda_2=lambda_1*gamma;

c_10=xi; % input concentration of the filtering solution
c_20=1-xi;
X=0:h:1; % X grid
n=length(X); % length of X
YY=1-X;

% take function handle as input for initial pore profile
for jj=1:n
    A(jj,1)=a_0(X(jj));   % initial pore shape for all x
end


%***Check the initial pore file for physics constraint ***%
if max(max(A))>1||min(min(A))<=0
    
    disp('physics constraint violated')
    
    f=1;
    
else
    %***If physics constraint check passed, calculate resistance R_m ***%
    for ii=1:n
        Rm(1)=h*trapz(1./(A(:,1).^4));
    end
    
    
    u(1)=1/Rm(1); % darcy velocity of the unit box
    
    
    flux(1)=u(1);
    throughput(1)=0;
    
    
    T(1)=0;
    
    C_1(1,1)=c_10; % concentration of particle 1 at the inlet at time 0
    C_2(1,1)=c_20; % concentration of particle 2 at the inlet at time 0
    
    % calculate velocity in the pore
    for jj=1:n
        up(jj,1)=4*u(1)/(pi*A(jj,1).^2);
    end
    
    
    % calculate the concentration change for each particle
    for jj=2:n
        C_1(jj,1)=(1-lambda_1*h/(A(jj-1,1)*up(jj-1,1)))*C_1(jj-1,1);
        C_2(jj,1)=(1-lambda_2*h/(A(jj-1,1)*up(jj-1,1)))*C_2(jj-1,1);
    end
    
    % initialize for the next time step calculation
    for jj=1:n
        A_old(jj)=A(jj,1);
        C_1_old(jj)=C_1(jj,1);
        C_2_old(jj)=C_2(jj,1);
    end
    
    % record the outlet concentration and accumulative concentration
    C_1_outlet(1)=C_1(n,1);
    C_2_outlet(1)=C_2(n,1);
    C_1_acm(1)=C_1(n,1);
    C_2_acm(1)=C_2(n,1);
    
    
    i=2;
    
    % ***Check if the membrane is out of service i.e. completely fouled
    % with low flux ***%
    while min(u)>0.1*u(1) && min(u)>0.0001 && min(A_old)>0
        
        
        for jj=1:n
            A_new(jj)=A_old(jj)-k*(C_1_old(jj)+beta*C_2_old(jj));
        end
        
        
        % check after fouling if the pore is closed
        if min(A_new)<=0
            disp('block happened');
            u(i)=0;
        else
            %***After pore evol, calc the new resist R_m ***%
            Rm_new=h*trapz(1./(A_new.^4)); % Rm at time kk
            u(i)=1/Rm_new; % calculate the flux
            
            
            % record flux, throughput and pore radius
            flux(i)=u(i);
            throughput(i)=k*trapz(flux(1:i));
            
            % Calcluate the u'(0) with da/dt
            flux_prime(i-1)=trapz((4./A_old.^5).*(A_new-A_old)./k)...
                /(trapz(1./A_old.^4))^2;
            
            % approximate the u'(0)
            %           flux_prime(1)=(flux(2)-flux(1))/k;
            
            for jj=1:n
                A(jj,i)=A_new(jj);
            end
            
            T(i)=k*(i-1);
            
            
            % concentration at the inlet x=0 for particle 1
            C_1_new(1)=c_10;
            % concentration at the inlet x=0 for particle 2
            C_2_new(1)=c_20;
            
            %calc pore velocity
            for jj=1:n
                up(jj,i)=4*u(i)/(pi*A_new(jj).^2);
            end
            
            for ii=1:n
                if min(A_old)<=0
                    
                    for jj=1:n
                        C_1_new(jj)=0;
                        C_2_new(jj)=0;
                    end
                    
                else
                    
                    for jj=2:n
                        C_1_new(jj)=(1-lambda_1*h/(A_new(jj-1)...
                            *up(jj-1,i)))*C_1_new(jj-1);
                        C_2_new(jj)=(1-lambda_2*h/(A_new(jj-1)...
                            *up(jj-1,i)))*C_2_new(jj-1);
                    end
                    
                end
            end
            
            %***record various quantities at step i ***%
            C_1_outlet(i)=C_1_new(n);
            C_2_outlet(i)=C_2_new(n);
            C_1_outlet_prime(i-1)=(C_1_outlet(i)-C_1_outlet(i-1))/k;
            C_2_outlet_prime(i-1)=(C_2_outlet(i)-C_2_outlet(i-1))/k;
            C_1_acm(i)=k*trapz(flux.*C_1_outlet)/throughput(i);
            C_2_acm(i)=k*trapz(flux.*C_2_outlet)/throughput(i);
            
            C_1_old=C_1_new; % update the iteration for the next step
            C_2_old=C_2_new; % update the iteration for the next step
            A_old=A_new;
            i=i+1;
        end
    end
    
    if obj1
        f(1)=W(1)*throughput(end)+W(2)*C_2_acm(end);
        f(2)=W(1)*flux(1)+W(1)*flux_prime(1)...
            +W(2)*C_2_outlet(1)+W(2)*C_2_outlet_prime(1);
        
    elseif obj2
        f(1)=throughput(end)*C_2_acm(end);
        f(2)=flux(1)*C_2_outlet(1)+...
            W(1)*(flux_prime(1)+C_2_outlet_prime(1))+...
            W(2)*flux_prime(1)*C_2_outlet_prime(1);
    end
    %     f(5)=flux(1)*C_2_outlet(1)...
    %         +(flux_prime(1)+C_2_outlet_prime(1));
    %     f(6)=flux(1)*C_2_outlet(1)...
    %         +flux_prime(1)*C_2_outlet_prime(1);
    
    
end


end



%