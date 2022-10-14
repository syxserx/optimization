% 2022-09-17 v1
%
% 3 inputs:     a_slope
%               a_intercept
%               num_pleats = number of pleats (scalar)
%
%	where pore_profile = a_slope*y + a_intercept on [0,1,] 0 = top, 1 =
%	bottom
%	e.g. a_intercept = 0.9, a_slope = -0.8 means top of pore is 0.9
%	and bottom of pore is 0.1 
%
% 2 outputs:   total_throughput (scalar) at t_f
%              final_cavg (scalar) at t_f
%
% Sample call for linear pore profile 0.9 at top and 0.1 at bottom: 
%   [total_throughput ] = pleated_filter(-0.8,0.9,90)
%
% Sample call for constant pore profile 0.9
%	[total_throughput ] = pleated_filter(0,0.9,90)
%
function f =pleated_filter_obj_fun(p,num_pleats,R)

% close all
% clear all
% clc

% Parameters
tic
time = 1;
phi = 0.5;
l = 0.5;
epstar = 0.5;
lambda = 2;
ga = (1-(1+epstar)^2)/log(1+epstar);
Pec = 1;           % For infinite Pe, let Pec > 100000
% num_pleats = 90;    % N = 2*pi/(2*theta_s) is # of pleats, so for whole cartridge
theta_s = pi/num_pleats;

% Discretization
iterr = 50;
iterz = 50;
itery = 50;
itertime = 50;
dr = (1-l)/iterr;
dy = 1/itery;
dz = 1/iterz;
dt = time/itertime;

% coeffcients for concentration c - diffusion advection equation
alp = -1/Pec;
bet = zeros(iterz+1,iterr+1,itery+1);
gam = zeros(iterz+1,iterr+1,itery+1);

%define function M(x)
M = @(x)  x*(x^2-2+ga*(2*log(x)-1))-(1+epstar)^2/x*((1+epstar)^2-2+ga*(2*log(1+epstar)-1));

%define function K and K'
% Note: To allow for easily input K parameters, K is not a function handle
% K = @(x) 1; %1-0.5*x;
K = 360/(2*pi)*theta_s;
KP = @(x) 0; %-0.5;

Gamma = 10*(1/K^2); % Now is a parameter
% return

% Pre-allocate memory
c = zeros(iterz+1,iterr+1,itery+1,itertime+1); % Concentration
% BC for concentration at the top of membrane
c(:,:,itery+1,:) = 1;
Ac = zeros(itery+1);
Bc = zeros(itery+1,1);
Bc(itery+1,1) = 1;

% Defining a(r,z,y,t) (pore radius)
% Initial pore radius
a0 = @(y) polyval(p,y);
a = zeros(iterz+1,iterr+1,itery+1,itertime+1);
aprime = zeros(iterz+1,iterr+1,itery+1);
for j=1:iterz+1
    for  i=1:iterr+1
        for k=1:itery+1
            a(j,i,k,1) = a0((k-1)*dy);
        end
    end
end

% % test linear profile 9/17/22
% aplot = zeros(k,1);
% for k = 1:itery+1
%     aplot(k) = a(1,1,k,1);
% end
% aplot
% return

% Integral: D = int_l^r dr'/r'k(r') using trapezoidal rule
D = zeros(iterr+1,1);
for i=1:iterr+1
    S=0;
    for j=1:i
        %S = S + 2 / ((l+(j-1)*dr) * K(l+(j-1)*dr));
        S = S + 2 / ((l+(j-1)*dr) * K);
    end
%     D(i) = (S - 1/(l*K(l)) - 1/((l+(i-1)*dr) * K((l+(i-1)*dr)))) * dr/2;
    D(i) = (S - 1/(l*K) - 1/((l+(i-1)*dr) * K)) * dr/2;
end

% Defining A & B coefficient of the BVP for P
A = zeros(iterr+1,1);
B = zeros(iterr+1,1);

for i=1:iterr+1
    %A(i) = (l+(i-1)*dr) * K(l+(i-1)*dr); % A = r*k
    A(i) = (l+(i-1)*dr) * K; % A = r*k
    %B(i) = (l+(i-1)*dr) * KP(l+(i-1)*dr) +  K(l+(i-1)*dr); % B = r*k'+k
    B(i) = (l+(i-1)*dr) * KP(l+(i-1)*dr) +  K; % B = r*k'+k
end


%BB right hand side vector (AA*XX=BB)
BB=zeros((iterr+6)*(iterz+3),1);
for j=1:iterz+3
    BB((iterr+4)+(iterr+6)*(j-1),1)=1;
end

% Initializing Psp and Psm
Psp = zeros(iterz+1,iterr+1,itertime+1);
Psm = zeros(iterz+1,iterr+1,itertime+1);

% Initializing B8, B9, F (F_zz = B8)
B8 = zeros(iterz+1,itertime+1);
B9 = zeros(iterz+1,itertime+1);
F = zeros(iterz+1,itertime+1);

% Initializing pressures Pe & Ph
Pe = zeros(iterz+1,itertime+1);
Ph = zeros(iterz+1,itertime+1);

% Initializing Vm and Vp
Vm = zeros(iterz+1,iterr+1,itertime+1);
Vp = zeros(iterz+1,iterr+1,itery+1,itertime+1);

% Initializing intA and C
intA = zeros(iterz+1,iterr+1,itertime+1);
C = zeros(iterz+1,iterr+1);

% Initializing vectors for the results 
porosity = zeros(iterz+1,iterr+1,itertime+1);
%Cbar = zeros(iterz+1,iterr+1,itertime+1);
cavg = zeros(itertime+1,1);
flux = zeros(itertime+1,1);
cvol = zeros(itertime+1,1);
throughput = zeros(itertime+1,1);
cacm = zeros(itertime+1,1);

logic = 1;
J = itertime;

% ---------------------------Time loop starts here-------------------------
for m=1:itertime+1
    if logic == 1
        
        % Z loop starts here
        for j=1:iterz+1
            
            % Solves intA and gets C coefficient
            % Integral: (int_(-1/2)^(1/2) dy/a^4(y)) using trapezoidal rule
            for i=1:iterr+1
                if sum(a(j,i,:,m)) ~= 0 %i.e. if the pore is not clogged at any point, since we set all cross sections of the pore zero, if at least one cross section is clogged (see line 246)
                    S(i) = 0;
                    for k=1:itery+1
                        S(i) = S(i) + 2/a(j,i,k,m)^4;
                        mmmmm(m,i) = S(i);
                    end
                    intA(j,i,m) = (S(i) - 1/a(j,i,1,m)^4 - 1/a(j,i,itery+1,m)^4)*dy/2;
                    C(j,i) = (-2) * Gamma/intA(j,i,m);
                else % sum(a(j,i,:,m)) = 0 i.e. the pore is not clogged at any point
                    C(j,i) = 0; % Since S & a=0, then C = 0; (Resistance = inf)
                end
            end
            
            
            % -------------------------------------------------------------
            % This for loop creates the submatrix that can be used for both Ps+
            % subMatrix (iterr+6)x(iterr+6)
            % interior nodes for Ps+, B8 and B9
            
            % Initializing AA matrix (AA*XX=BB)
            AA = zeros((iterr+1+2+3)*(iterz+1+2)); % 2 because of ghost points, 3 because of B8, B9 & FF
            subMatrix = zeros(iterr+6,iterr+6);
            XX = zeros((iterr+1+2+3)*(iterz+1+2),1);
            
            % 4/17/21 Fixed coeffs for B_8 & B_9
            for i = 1:iterr+1
                subMatrix(i,i) = A(i)/dr^2 - B(i)/(2*dr);  % coeff p(i-1)
                subMatrix(i,i+1) = (-2)*A(i)/dr^2 + C(j,i);  % coeff p(i)
                subMatrix(i,i+2) = A(i)/dr^2 + B(i)/(2*dr);  % coeff p(i+1)
                subMatrix(i,iterr+4) = -C(j,i)/2*D(i);         % coeff B_8
                subMatrix(i,iterr+5) = -C(j,i)/2;              % coeff B_9
            end
         end
        % Z loop ends
        
        % BC1: Left ghost point - row 1
        subMatrix(iterr+2,1) = 1;
        subMatrix(iterr+2,3) = -1;
        
        % BC2: Right ghost point - row iterr+3
        subMatrix(iterr+3,iterr+1) = -1/(2*dr);
        subMatrix(iterr+3,iterr+3) = 1/(2*dr);
        %subMatrix(iterr+3,iterr+4) = -1/K(1);
        subMatrix(iterr+3,iterr+4) = -1/K;
        
        % BC3: - row iterr+4
        subMatrix(iterr+4,iterr+2) = 1;
        subMatrix(iterr+4,iterr+6) = -16*phi/M(1);
        
        % BC4: - row iterr+5
        subMatrix(iterr+5,2) = -1;
        subMatrix(iterr+5,iterr+5) = 1;
        subMatrix(iterr+5,iterr+6) = 16*phi/l^4;
        
        % -------------------------------------------------------------
        
        % Replicate submatrix into AA for 1+(iterr+6)*(j-1):(iterr+6)*j
        for jj = 1:iterz+3
            AA(1+(iterr+6)*(jj-1):(iterr+6)*jj,1+(iterr+6)*(jj-1):(iterr+6)*jj) = subMatrix(1:iterr+6,1:iterr+6);
        end
        
        %BVP for F
        for jj = 1:iterz+1
            AA((iterr+6)*jj,(iterr+6)*jj) = -1/dr^2;
            AA((iterr+6)*jj,(iterr+6)*(jj+1)) = 2/dr^2;
            AA((iterr+6)*jj,(iterr+6)*(jj+2))= -1/dr^2;
            AA((iterr+6)*jj,(iterr+6)*(jj+1)-2) = 1;
        end
        
        %BVP for F: Coeff of BC: F^1=0
        AA((iterr+6)*(iterz+2),(iterr+6)*2) = 1;
        
        %BVP for F: Coeff of BC: F^0-F^2=0
        AA((iterr+6)*(iterz+3),iterr+6) = 1;
        AA((iterr+6)*(iterz+3),(iterr+6)*3) = -1;
        
        
        %BVP for P: F^{m+1} Coeff in BCs 3 and 4
        for jj=1:iterz+3
            AA((iterr+4)+(iterr+6)*(jj-1),(iterr+6)*(iterz+2))=AA((iterr+4)+(iterr+6)*(jj-1),(iterr+6)*(iterz+2))+16*phi/M(1);
            AA((iterr+5)+(iterr+6)*(jj-1),(iterr+6)*(iterz+2))=AA((iterr+5)+(iterr+6)*(jj-1),(iterr+6)*(iterz+2))-16*phi/l^4;
        end
        
        
        %Solving AA*XX=BB
        XX = AA\BB;
        
        %Assigning XX elements to PPsp, BB8, BB9 and FF
        for jj=1:iterz+3
            for i=1:iterr+3
                PPsp(jj,i)=XX((jj-1)*(iterr+6)+i,1);
            end
            BB8(jj)=XX((jj-1)*(iterr+6)+iterr+4,1);
            BB9(jj)=XX((jj-1)*(iterr+6)+iterr+5,1);
            FF(jj)=XX((jj-1)*(iterr+6)+iterr+6,1);
        end
        
        % Removing ghost points and put into matrix with time
        Psp(1:iterz+1,1:iterr+1,m)=PPsp(2:iterz+2,2:iterr+2);
        B8(1:iterz+1,m)=BB8(2:iterz+2);
        B9(1:iterz+1,m)=BB9(2:iterz+2);
        F(1:iterz+1,m)=FF(2:iterz+2);
        
        % From writeup - Section 4.6, eqn(98)
        for jj=1:iterz+1
            for i=1:iterr+1
                Psm(jj,i,m)=-Psp(jj,i,m)+B8(jj,m)*D(i)+B9(jj,m);
            end
        end
        
        % Concentration - a & c
        for j = 1:iterz+1
            % From writup - Section 4.3, eqn(74) & eqn(75)
            for i=1:iterr+1
                Vm(j,i,m) = (Psm(j,i,m)-Psp(j,i,m))/intA(j,i,m); % intA is 0 at certain points, so Vm becomes Inf
                for k=1:itery+1
                    Vp(j,i,k,m) = (4*Vm(j,i,m)/(pi*a(j,i,k,m)^2));
                end
            end
            
            % Advection Diffusion Eqn with Pe -> Infinity
            % Eqn 76 for c(r,z,y,t)
            % Start at top-down , y = 1/2 (itery+1)
            if Pec > 100000
                for i=1:iterr+1
                    for k=1:itery
                        c(j,i,itery+1-k,m) = c(j,i,itery+2-k,m)/ ...
                            (1+(lambda*dy)/(abs(Vp(j,i,itery+1-k,m))*a(j,i,itery+1-k,m)));
                    end
                end
            else
                %
                % Modified 4/19/21
                % Advection Diffusion Eqn with Pe
                % alp c'' + bet c' + gam c = 0
                % alp = 1/Pec
                % Gc*c(i-1) + Fc*c(i)+ Ec*c(i+1)

                % Calculate a prime 
                for i=1:iterr+1
                    for k=2:itery+1
                        aprime(j,i,k) = (a(j,i,k,m) - a(j,i,k-1,m))/dy;
                    end
                end

                % bet - coeff for c' 
                for i=1:iterr+1
                    for k=2:itery+1
                        bet(j,i,k) = Vp(j,i,k,m) - 2/Pec*aprime(j,i,k)/a(j,i,k,m);
                    end
                end

                % gam - coeff for c
                for i=1:iterr+1
                    for k=1:itery+1
                        gam(j,i,k) = lambda/a(j,i,k,m);
                    end
                end

                % Matrix problem: Ac*c=Bc, solve for c
                % Fill the matrix Ac
                for i=1:iterr+1
                    for k=2:itery
                        Ac(k,k-1) = alp/dy^2 - bet(j,i,k)/(2*dy);
                        Ac(k,k) = -2*alp/dy^2 + gam(j,i,k);
                        Ac(k,k+1) = alp/dy^2 + bet(j,i,k)/(2*dy);
                    end
                    % BC
                    % dc/dy = 0 at y = -1/2
                    Ac(1,1) = -2*alp/dy^2 + gam(j,i,1);
                    Ac(1,2) = 2*alp/dy^2;
                    % c = 1 at y = 1/2
                    Ac(itery+1,itery+1) = 1;

                    % Ac*c = Bc
                    c(j,i,:,m) = Ac\Bc;
                end
            
            end % end for if Pe > 1 
            

            % Pore radius a(r,z,y,t)
            for i=1:iterr+1
                for k=1:itery+1
                    a(j,i,k,m+1) = a(j,i,k,m)-c(j,i,k,m)*dt; % Eqn 42a
                    if a(j,i,k,m+1) <= 0 % Pore radius cannot be negative
                        a(j,i,1:(itery+1),(m+1):(itertime+1)) = 0; % Making negative values = 0
                        break
                    end
                end
            end
           
         end
        
        % Pressures in hollow and empty region - function of z and t
        for jj=1:iterz+1
            Ph(jj,m) = -16*phi/l^4*(F(jj,m)-F(iterz+1,m));
            Pe(jj,m) = 16*phi/M(1)*(F(jj,m)-F(iterz+1,m)) + 1;
        end
        
        % If pore radius at later times == 0, logic = 0 and loops end
        % When entire membrane is clogged, counter = 0
        counter = 0;
        for j=1:iterz+1
            for i=1:iterr+1
                for k=1:itery+1
                    counter = counter + a(j,i,k,m+1);
                end
            end
        end
        
        if counter == 0
            logic = 0;
            J = m; % Final time that we stop the loop at
        end
        
          % 9/21/22
        % Stops code when current flux is 10% of initial flux
        sf = 0;
        sc = 0;
        for j=2:iterz
            for i=2:iterr
                sf = sf + 2*abs(Vm(j,i,m));
                sc = sc + 2*abs(Vm(j,i,m)).*c(j,i,1,m);
            end
        end

        % flux = \int_0^one \int_l^one v_m dr dz
        flux(m) = (sf + sum(abs(Vm(1,:,m))) + sum(abs(Vm(iterz+1,:,m))) ...
                + sum(abs(Vm(:,1,m))) + sum(abs(Vm(:,iterr+1,m))) ...
                - sum(abs(Vm(1,1,m))) - sum(abs(Vm(1,iterr+1,m))) ...
                - sum(abs(Vm(iterz+1,1,m))) - sum(abs(Vm(iterz+1,iterr+1,m))))*(dr/2)*(dz/2);

        cvol(m) = (sc + sum(abs(Vm(1,:,m)).*c(1,:,1,m)) + sum(abs(Vm(iterz+1,:,m).*c(iterz+1,:,1,m))) ...
                + sum(abs(Vm(:,1,m)).*c(:,1,1,m)) + sum(abs(Vm(:,iterr+1,m).*c(:,iterr+1,1,m))) ...
                - sum(abs(Vm(1,1,m).*c(1,1,1,m))) - sum(abs(Vm(1,iterr+1,m).*c(1,iterr+1,1,m))) ...
                - sum(abs(Vm(iterz+1,1,m).*c(iterz+1,1,1,m))) - sum(abs(Vm(iterz+1,iterr+1,m).*c(iterz+1,iterr+1,1,m))))*(dr/2)*(dz/2);

        cavg(m) = cvol(m)/flux(m);
           
        if abs(flux(m) - 0.1*flux(1)) < 10^(-3)
            logic = 0;
            disp('pore closed')
        end
        % -----------------------------------------------------------------
       
        % -----------------------------------------------------------------
        % 9/21/22
        % Stops code when initial average concentration in the filtrate
        % (at the outlet) is greater 10%
        if cavg(1) > 1-R
            logic = 0;
            f=2
            disp('particle removal not satisfied')
        end
        % -----------------------------------------------------------------
    
        % Results - added 4/26/21
        % porosity int_{-1/2}^{1/2} pi a^2 dy
        % porosity = zeros(iterz+1,iterr+1,itertime+1);
        % Cbar = zeros(iterz+1,iterr+1,itertime+1);
        
        for j=1:iterz+1
            for i=1:iterr+1
                s = 0;
                Cs = 0;
                for k=1:itery+1
                    s = s + 2*pi/4*a(j,i,k,m)^2;
                    Cs = Cs + 2*c(j,i,k,m);
                end
                
                porosity(j,i,m) = (s - pi/4*a(j,i,1,m)^2 - pi/4*a(j,i,itery+1,m)^2)*dy/2; 
                Cbar(j,i,m) = (Cs - c(j,i,1,m) - c(j,i,itery+1,m))*dy/2;
            end
        end
        
    end % If logic == 1 ends
    
end
% --------------------------- time loop ends ------------------------------


tf = J*dt;


% Psp
% Psm
% c
% a
toc

% mmmmm
% Includes only the times of values when pore radius is not zero.
aa(:,:,:,:) = a(:,:,:,1:J);
cc(:,:,:,:) = c(:,:,:,1:J);
Pspp(:,:,:) = Psp(:,:,1:J);
Psmm(:,:,:) = Psm(:,:,1:J);
Phh(:,:) = Ph(:,1:J);
Pee(:,:) = Pe(:,1:J);
Pporosity(:,:,:) = porosity(:,:,1:J);
avg_perm(:,:,:) = 1./intA(:,:,1:J);

% Results - added 4/26/21
% Capture Concentration at y = -1/2 (at k = 1 in the code)
%               int_0^1 int_l^1 vm(z,r,t)*c(z,r,y=-1/2,t) dr dz
% cavg(t) = -----------------------------------------------------
%               int_0^1 int_l^1 vm(z,r,t) dr dz
%
%               int_0^t int_0^1 int_l^1 vm(z,r,t)*c(z,r,y=-1/2,t) dr dz dt
% cacm(t) = -----------------------------------------------------
%                  int_0^t int_0^1 int_l^1 vm(z,r,t) dr dz dt

% % (pi/theta_s) is # of pleats, so for whole cartridge
% theta_s = 2*pi/90;
Vm = (pi/theta_s)*Vm;

for m = 1:J-1
            sf = 0;
            sc = 0;
    for j=2:iterz
        for i=2:iterr
            sf = sf + 2*abs(Vm(j,i,m));
            sc = sc + 2*abs(Vm(j,i,m)).*c(j,i,1,m);
        end
    end
    
    % flux = \int_0^one \int_l^one v_m dr dz 
    flux(m) = (sf + sum(abs(Vm(1,:,m))) + sum(abs(Vm(iterz+1,:,m))) ...
            + sum(abs(Vm(:,1,m))) + sum(abs(Vm(:,iterr+1,m))) ...
            - sum(abs(Vm(1,1,m))) - sum(abs(Vm(1,iterr+1,m))) ...
            - sum(abs(Vm(iterz+1,1,m))) - sum(abs(Vm(iterz+1,iterr+1,m))))*(dr/2)*(dz/2); 
        
    cvol(m) = (sc + sum(abs(Vm(1,:,m)).*c(1,:,1,m)) + sum(abs(Vm(iterz+1,:,m).*c(iterz+1,:,1,m))) ...
            + sum(abs(Vm(:,1,m)).*c(:,1,1,m)) + sum(abs(Vm(:,iterr+1,m).*c(:,iterr+1,1,m))) ...
            - sum(abs(Vm(1,1,m).*c(1,1,1,m))) - sum(abs(Vm(1,iterr+1,m).*c(1,iterr+1,1,m))) ...
            - sum(abs(Vm(iterz+1,1,m).*c(iterz+1,1,1,m))) - sum(abs(Vm(iterz+1,iterr+1,m).*c(iterz+1,iterr+1,1,m))))*(dr/2)*(dz/2);
    
	cavg(m) = cvol(m)/flux(m);
    
    st = 0;
    scc = 0;
    for mm =1:m
        st = st + 2*flux(mm);
        scc = scc + 2*cvol(mm);
    end
    
    % throughput = \int_0^t \int_0^one \int_l^one v_m dr dz dt
    % integral of flux over time
    throughput(m) = (st - flux(1) - flux(m))*dt/2;
    
    cacm(m) = (scc - cvol(1) - cvol(m))*dt/2 / throughput(m);
    
end


% 7/14/2021
% Get rid of zero value in the flux vector
% size(flux)
% size(throughput)
flux = flux(flux~=0);
% throughput = throughput(throughput~=0);
throughput(length(flux)+1:end) = [];
cavg(length(flux)+1:end) = [];
cacm(length(flux)+1:end) = [];

total_throughput=max(throughput);
initial_cavg = cavg(1);
final_cavg = cavg(end);

% if final_cavg>1-R
%     f=2;
% else
    f=-total_throughput
% end

end







