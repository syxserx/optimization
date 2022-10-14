function  f=flux_objective_fun8(p,gamma,lambda,R)

a_0=@(y) polyval(p,y);
%    a_0=@(y) p(1).*(y.^2)+p(2).*y+p(3);
%  a_0_x=@(x,y) (0.5*x+0.5)*(p(1)*y+p(2));
%  a_0_x=@(x,y) (p(1)*x+p(2))*(p(3)*y+p(4));


% gamma=5;
% lambda=.1;
% R=0.999;




c_half_t=1; % input concentration of the filtering solution
h=1/40; % y grid size
hx=h;
hy=h;
k=0.001; %time grid siz

% T=0:k:t; % time grid
X=0:hx:1; % X grid
Y=-0.5:hy:0.5; % this may need to change to make the computation easier
YY=-Y;

nx=length(X);
ny=length(Y);

% A=zeros(length(X),length(Y),length(T));
% Rm=zeros(length(X),length(T));
% Vm=zeros(length(X),length(T));
% Po=zeros(length(X),length(T));
% C=zeros(length(X),length(Y),length(T));



for ii=1:nx
    for jj=1:ny
        A(ii,jj,1)=a_0(YY(jj));   % initial pore shape for all x,y
    end
end

if max(max(A))>1||min(min(A))<=0
    
%      disp(['p1=' num2str(p1), '    p2=' num2str(p2) ' makes the pore larger than the box']);
    
    f=1;
    
else
    
    for ii=1:nx
        
        Rm(ii,1)=0.5*hy*trapezoid_discrete(1./(A(ii,:,1).^4));
        %disp(trapezoid_discrete(A(ii,:,kk)))
    end
    
    % disp(a_0);
    % disp(Rm(1,1));
    
    % disp(Rm(:,1));
    P0=P0_solver(Rm(:,1),hx,gamma);
    % disp(Rm(:,1));
    
    Po(1,1)=1;
    c1=P0(1);
    c2=P0(2);
    Po(nx,1)=c1+c2;
    
    for ii=2:nx-1
        Po(ii,1)=P0(ii+1); % find Po at time 0
    end
    
    for ii=1:nx
        Vm(ii)=(2*Po(ii,1)-c1*(ii-1)*hx-c2)/Rm(ii,1); % add vm to calculate the flux
        p0_p_m(ii,1)=2*Po(ii,1)-c1*(ii-1)*hx-c2; % check if the symmetry comes from p_0^+-p_0^-
    end
    
    flux(1)=0.5*hx*trapezoid_discrete(Vm);
    
    throughput(1)=0;
    
    T(1)=0;
    
    
    
    for ii=1:nx
        
        C(ii,1,1)=c_half_t; % concentration at the inlet at time 0
    end
    
    for ii=1:nx
        
        for jj=2:ny
            
            C(ii,jj,1)=C(ii,jj-1,1)/(1+(hy*lambda*A(ii,jj,1)*Rm(ii,1))/(2*Po(ii,1)-c1*(ii-1)*hx-c2)); % C at time 0
            
        end
        
    end
    
    for ii=1:nx
        for jj=1:ny
            A_old(ii,jj)=A(ii,jj,1);
            C_old(ii,jj)=C(ii,jj,1);
        end
    end
    
    %     C_half(:,1)=C(20,:,1);
    %     C_1(:,1)=C(1,:,1);
    C_outlet(:,1)=C(:,ny,1);
    C_c=Vm.*(C_outlet(:,1))';
    C_avg(1)=0.5*hx*trapezoid_discrete(C_c)/flux(1);
    fltrate_particle_content(1)=0.5*hx*trapezoid_discrete(C_c);
    C_acm(1)=0.5*hx*trapezoid_discrete(C_c)/flux(1);
    c_simple_avg(1)=0.5*hx*trapezoid_discrete(C_outlet(:,1));
    
    i=2;
    
    if C_avg(1)>1-R
        f=2;
    else
        
        while min(flux)>0.001
            
            %      disp(min(flux));
            
            
            
            for ii=1:nx
                if min(A_old(ii,:))<=0
                    
                    for jj=1:ny
                        A_new(ii,jj)=A_old(ii,jj);
                    end
                    
                else
                    for jj=1:ny
                        
                        A_new(ii,jj)=A_old(ii,jj)-k*C_old(ii,jj);  % at t=kk the Radius can be obtained first
                        
                    end
                end
                
            end
            
            if min(A_new)<0
                disp('block happened');
            end
            
            for ii=1:length(X)
                
                Rm_new(ii)=0.5.*hy.*trapezoid_discrete(1./(A_new(ii,:).^4)); % Rm at time kk
                %disp(trapezoid_discrete(A(ii,:,kk)))
            end
            
            P0=P0_solver(Rm_new,hy,gamma);
            
            Po(1)=1;
            c1=P0(1);
            c2=P0(2);
            Po(nx)=c1+c2;
            
            for ii=2:nx-1
                Po(ii)=P0(ii+1); % find Po at time kk
            end
            
            
            for ii=1:nx
                
                p0_p_m(ii,i)=2*Po(ii)-c1*(ii-1)*hx-c2;
                
                if min(A_new(ii,:)<0)
                    Vm(ii)=0;
                else
                    Vm(ii)=(2*Po(ii)-c1*(ii-1)*hx-c2)/Rm_new(ii); % add vm to calculate the flux
                    
                end
            end
            
            flux(i)=0.5*hx*trapezoid_discrete(Vm);
            
            throughput(i)=0.5*k*trapezoid_discrete(flux(1:i));
            
            T(i)=k*(i-1);
            
            for ii=1:nx
                C_new(ii,1)=c_half_t; % concentration at the inlet y_tilda=0.5
            end
            
            for ii=1:nx
                if min(A_old(ii,:))<=0
                    
                    for jj=1:ny
                        C_new(ii,jj)=0;
                    end
                    
                else
                    
                    for jj=2:length(Y)
                        C_new(ii,jj)=C_old(ii,jj-1)/(1+(hy*lambda*A_new(ii,jj)*Rm_new(ii))/(2*Po(ii)-c1*(ii-1)*hx-c2)); % C at time kk
                    end
                    
                end
            end
            
            %             C_half(:,i)=C_new(20,:); % concentration profile for pore no.20 we call it middle
            %             C_1(:,i)=C_new(1,:); % concentration profile for pore no.1
            C_outlet(:,i)=C_new(:,ny); % concentration at the outlet y_tilda=-0.5
            C_c=Vm.*(C_outlet(:,i))';
            C_avg(i)=0.5*hx*trapezoid_discrete(C_c)/flux(i);
            fltrate_particle_content(i)=0.5*hx*trapezoid_discrete(C_c);
            C_acm(i)=0.5*k*trapezoid_discrete(fltrate_particle_content)/throughput(i);
            c_simple_avg(i)=0.5*hx*trapezoid_discrete(C_outlet(:,i));
            
            C_old=C_new; % update the iteration for the next step
            A_old=A_new;
            i=i+1;
            
        end
        
        f=-10*throughput(end)+0*C_avg(1);
    end
    
    
end


% value at each points in the final time

function V=P0_solver(Rm_t,h,gamma)

n=length(Rm_t);

b=zeros(1,n);

b(1)=2*(Rm_t(1)/(h^2)+gamma);
b(2)=-Rm_t(2)/(h^2);

% disp(b);

A=zeros(n);

A(1,1)=-2*Rm_t(1)/h;
A(1,2)=gamma;
A(1,3)=2*Rm_t(1)/(h^2);
A(2,1)=gamma*h;
A(2,2)=gamma;
A(2,3)=-2*(Rm_t(2)/(h^2)+gamma);
A(2,4)=Rm_t(2)/h^2;

A(n-1,1)=gamma*(n-2)*h+Rm_t(n-1)/(h^2);
A(n-1,2)=gamma+Rm_t(n-1)/(h^2);
A(n-1,n-1)=Rm_t(n-1)/(h^2);
A(n-1,n)=-2*(Rm_t(n-1)/(h^2)+gamma);
A(n,1)=gamma*(n-1)*h-2*(Rm_t(n)/(h^2)+gamma);
A(n,2)=gamma-2*(Rm_t(n)/(h^2)+gamma);
A(n,n)=2*Rm_t(n)/(h^2);

for i=3:n-2
    
    A(i,1)=gamma*(i-1)*h;
    A(i,2)=gamma;
    A(i,i)=Rm_t(i)/(h^2);
    A(i,i+1)=-2*(Rm_t(i)/(h^2)+gamma);
    A(i,i+2)=Rm_t(i)/(h^2);
    
end
% disp(A)

V=A\b';

function I=trapezoid_discrete(M)

I=M(1)+M(length(M));

for i=2:length(M)-1
    sum=2.*M(i);
    I=I+sum;
end


%%%
%%
 