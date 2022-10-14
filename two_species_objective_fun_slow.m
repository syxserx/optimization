function  f=two_species_objective_fun_slow(p,xi, beta, gamma,lambda_1,R,W)

obj1=W(4);
obj2=1-W(4);

h=0.01; % space grid size
k=0.001; %time grid siz

% lambda_1=1;
lambda_2=lambda_1*gamma;

c_10=xi; % input concentration of the filtering solution
c_20=1-xi;
x1=0:h:p(1);
x2=p(1):h:1;



for jj=1:length(x1)
    A(jj,1)=x1(jj)*(p(2)-1)/p(1)+1; % initial pore shape for all x
end

for jj=length(x1):length(x1)+length(x2)
    A(jj,1)=p(2); % initial pore shape for all x
end

%***Check the initial pore file for physics constraint ***%
if max(A)>1||min(A)<=0
    
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
    
    for jj=1:n
        up(jj,1)=4*u(1)/(pi*A(jj,1).^2);
    end
    
    
    for jj=2:n
        C_1(jj,1)=(1-lambda_1*h/(A(jj-1,1)*up(jj-1,1)))*C_1(jj-1,1);
        C_2(jj,1)=(1-lambda_2*h/(A(jj-1,1)*up(jj-1,1)))*C_2(jj-1,1);
    end
    
    for jj=1:n
        A_old(jj)=A(jj,1);
        C_1_old(jj)=C_1(jj,1);
        C_2_old(jj)=C_2(jj,1);
    end
    
    C_1_outlet(1)=C_1(n,1);
    C_2_outlet(1)=C_2(n,1);
    C_1_acm(1)=C_1(n,1);
    C_2_acm(1)=C_2(n,1);
    
    i=2;
    
    % %***check if the particle removal requirement is violated ***
    if abs(C_1_outlet(1))>=(1-R(1))*c_10
%         ||u(1)<0.01
        f=2;
        disp('particle removal constraint is violated');
    else
     
        %***Check if the membrane is out of service i.e. completely fouled
        %  with low flux ***%
        while min(A_old)>0 && min(u)>0.1*u(1) && min(u)>0.001
            
            for jj=1:n
                A_new(jj)=A_old(jj)-k*(C_1_old(jj)+beta*C_2_old(jj));
            end
            
            if min(A_new)<=0
                disp('block happened');
                u(i)=0;
            else
                %***After pore evolution, calculate the new resistance R_m ***%
                Rm_new=h*trapz(1./(A_new.^4)); % Rm at time kk
                u(i)=1/Rm_new; % add vm to calculate the flux
            end
            
            
            flux(i)=u(i);
            
            throughput(i)=k*trapz(flux(1:i));
            
            T(i)=k*(i-1);
            
            
            C_1_new(1)=c_10; % concentration at the inlet x=0
            C_2_new(1)=c_20;
            
            for jj=1:n
                up(jj,i)=4*u(i)/(pi*A_new(jj).^2);
            end
            
            for jj=2:n
                C_1_new(jj)=(1-lambda_1*h/...
                    (A_new(jj-1)*up(jj-1,i)))*C_1_new(jj-1);
                C_2_new(jj)=(1-lambda_2*h/...
                    (A_new(jj-1)*up(jj-1,i)))*C_2_new(jj-1);
            end
            
            
            %***record various quantities at step i ***%
            C_1_outlet(i)=C_1_new(n);
            C_2_outlet(i)=C_2_new(n);
            C_1_acm(i)=k*trapz(flux.*C_1_outlet)/throughput(i);
            C_2_acm(i)=k*trapz(flux.*C_2_outlet)/throughput(i);
            
            C_1_old=C_1_new; % update the iteration for the next step
            C_2_old=C_2_new; % update the iteration for the next step
            A_old=A_new;
            i=i+1;
            
        end
        
        if obj1
            
         f=-W(1)*throughput(end)-W(2)*C_2_acm(end);
         
        elseif obj2
       
         f=-throughput(end)*C_2_acm(end);
         
        end
        
    end
    
end

