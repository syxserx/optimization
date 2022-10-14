
%% Constant pressure
close all,clear all,clc;

rng default % For reproducibility
ms = MultiStart;
gs=GlobalSearch;
options = optimset('Largescale','off','Display','iter');

% k=0.001;
% h=1/40;
% xi=0.5;
% beta=0.5; 
% lambda_1=1;
% R=[0.5, 0.9];

% u(0)+u'(0)+c2(0)+c2'(0)
   
M=[
 % W(1), W(2), F(1)/S(0), obj1(1)/obj2(0), R_1, xi, beta, lambda_1, sp;
  
  %obj1 Fast
 
%    0, 1, 1, 1, 0.99, 0.5, 0.1, .1, 10000;
%    1, 0, 1, 1, 0.99, 0.5, 0.1, .1, 10000;
%    0.5, 0.5, 1, 1, 0.99, 0.5, 0.1, .1, 10000;
   
  %obj2 Fast
   0, 0, 0, 0, 0.99, 90, 0.1, 1, 10;
%    0, 0, 0, 0, 0.99, 90, 0.1, 1, 1000;
%    0, 0, 0, 0, 0.99, 90, 0.1, 1, 10000;

  %obj1 Slow
 
%    1, 0, 0, 1, 0.99, 0.5, 0.1, .1, 10000;
%    0, 1, 0, 1, 0.99, 0.5, 0.1, .1, 10000;
%    0.5, 0.5, 0, 1, 0.99, 0.5, 0.1, .1, 10000;
    
  %obj2 Slow
%    0, 0, 0, 0, 0.99, 0.5, 0.1, 1, 10000;
   
  ];


for ii=1:size(M,1)
    W(1)=M(ii,1);
    W(2)=M(ii,2);
    W(3)=M(ii,3);
    W(4)=M(ii,4);
    R=M(ii,5);
    num_pleats=M(ii,6);
    beta=M(ii,7);
    lambda_1=M(ii,8);
    sp=M(ii,9);
    

n=2;
x0=0.5.*ones(1,n);
lb=-1*ones(1,n);
ub=1.*ones(1,n);



problem = createOptimProblem('fmincon','x0',x0,'nonlcon',...
    @(p)unitdisk1(p),'objective',@(p) ...
  pleated_filter_obj_fun_1(p,num_pleats,R,10),...
    'lb',lb,'ub',ub,'options',optimset);

ms.UseParallel='always';

tic;[xmin,fmin,flag,outpt,allmins] = run(ms,problem,sp); runtime=toc;

% fig=figure('Visible','off')
% disp(allmins)
% disp(xmin)
% disp(fmin)

% figure(1)
% plot(arrayfun(@(x)x.Fval,allmins),'k*')
% xlabel('Solution number')
% ylabel('Function value')
% title(['num_pleats=' num2str(xi), '  R=' num2str(R), '  p=' num2str(xmin)...
%     ' fmin=' num2str(fmin) ], 'interpreter','latex')
    
%get the directory of your input files:
[pathstr,namestr,extstr] =...
    fileparts('/home/y/ys379/Test');

% Create the folder you want the file to be saved to 
FolderDestination= strcat('beta=', num2str(beta),...
    ', R=', num2str(R), 'sp=',num2str(sp));
mkdir ('/home/y/ys379/Test',FolderDestination)

% create new file name for each run
 Filename=strcat('R=', num2str(R),',','sp=',num2str(sp));
 
%use that when you save
matfile = fullfile('/home/y/ys379/Test',FolderDestination,...
    [Filename,'.mat']);
% pdffile = fullfile('/Users/sunyx/Desktop',FolderDestination,...
%     [Filename,'.pdf']);
% figfile = fullfile('/Users/sunyx/Desktop', FolderDestination,...
%     [Filename,'.fig']);

save(matfile);
% saveas(figure(1),pdffile);
% saveas(figure(1),figfile);

end

%% constant flux
% close all, clear;
% 
% rng default % For reproducibility
% ms = MultiStart;
% gs=GlobalSearch;
% options = optimset('Largescale','off','Display','iter');
% 
% % k=0.001;
% % h=1/40;
% % xi=1;
% 
% % J*c2acm optimization
% 
% 
% lambda_1=10;
% R=[0.99, 0.9];
% % N=1000;
% 
% M=[% xi, beta, N, d1,  d0,   sp, 
%       0.1, 0.1, 1000, 0.3, 0.3, 10000;
%       0.5, 0.1, 1000, 0.3, 0.3, 10000;
%       0.9, 0.1, 1000, 0.3, 0.3, 10000;
%       0.1, 0.1,  500, 0.3, 0.3, 10000;
%       0.5, 0.1,  500, 0.3, 0.3, 10000;
%       0.9, 0.1,  500, 0.3, 0.3, 10000;
%       0.1, 0.1,  100, 0.3, 0.3, 10000;
%       0.5, 0.1,  100, 0.3, 0.3, 10000;
%       0.9, 0.1,  100, 0.3, 0.3, 10000;
%   ];
% 
% 
% for ii=1:size(M,1)
%     xi=M(ii,1);
%     beta=M(ii,2);
%     N=M(ii,3);
%     x0(1)=M(ii,4);
%     x0(2)=M(ii,5);
%     sp=M(ii,6);
%     
% gamma=beta;  
%     
% n=2;
% % x0=0.3.*ones(1,n);
% % lb=-1.*ones(1,n);
% % ub=1.*ones(1,n);
% lb=[-1,0];
% ub=[1,1];
% 
% 
% problem = createOptimProblem('fmincon','x0',x0,'nonlcon',...
%     @(p)unitdisk1(p),'objective',@(p) ...
%     two_species_objective_fun_cflux_1(p,xi, beta, gamma,lambda_1,R, N),...
%     'lb',lb,'ub',ub,'options',optimset);
% 
% ms.UseParallel='always';
% 
% tic;[xmin,fmin,flag,outpt,allmins] = run(ms,problem,sp); runtime=toc;
% 
% % fig=figure('Visible','off')
% % disp(allmins)
% % disp(xmin)
% % disp(fmin)
% 
% 
% 
% 
% figure(1)
% plot(arrayfun(@(x)x.Fval,allmins),'k*')
% xlabel('Solution number')
% ylabel('Function value')
% title([' $\xi=$' num2str(xi), ' $\beta=$' num2str(beta),...
%     ' $\gamma=$' num2str(gamma) '  R=' num2str(R), '  p=' num2str(xmin)...
%     ' fmin=' num2str(fmin) ], 'interpreter','latex')
% 
%     
% %get the directory of your input files:
% [pathstr,namestr,extstr] =...
%     fileparts('/Users/sunyx/Documents/MATLAB/two_species');
% 
% % Create the folder you want the file to be saved to 
% FolderDestination= strcat('cflux_J_c2acm',', $xi=$',...
%     num2str(xi),', beta=',num2str(beta),...
%     ', N=',num2str(N), ', x0=[',num2str(x0(1)),num2str(x0(2)), ']');
% mkdir ('/Users/sunyx/Desktop',FolderDestination)
% 
% % create new file name for each run
%  Filename=strcat('xi=',num2str(xi),'beta=',num2str(beta));
%  
% %use that when you save
% matfile = fullfile('/Users/sunyx/Desktop',FolderDestination,...
%     [Filename,'.mat']);
% pdffile = fullfile('/Users/sunyx/Desktop',FolderDestination,...
%     [Filename,'.pdf']);
% figfile = fullfile('/Users/sunyx/Desktop', FolderDestination,...
%     [Filename,'.fig']);
% 
% save(matfile);
% saveas(figure(1),pdffile);
% saveas(figure(1),figfile);
% 
% end
% 
% %% Multi-species
% 
% close all,clear all,clc;
% 
% rng default % For reproducibility
% ms = MultiStart;
% gs=GlobalSearch;
% options = optimset('Largescale','off','Display','iter');
% 
% xi=[0.7,0.15,0.15];
% beta=[1,0.1,0.5]; 
% lambda_1=1;
% R=[0.99, 0.5, 0.9];
%     
% n=2;
% x0=0.3.*ones(1,n);
% lb=-1.*ones(1,n);
% ub=1.*ones(1,n);
% sp=10000;
% 
%  problem = createOptimProblem('fmincon','x0',x0,'nonlcon',...
%      @(p)unitdisk1(p),'objective',@(p) ...
%      two_species_objective_fun_fast_m(p,xi, beta,lambda_1,R),...
%      'lb',lb,'ub',ub,'options',optimset);
%  
%  ms.UseParallel='always';
%  
%  tic;[xmin,fmin,flag,outpt,allmins] = run(ms,problem,sp); runtime=toc;
% 
% figure(1)
% plot(arrayfun(@(x)x.Fval,allmins),'k*')
% xlabel('Solution number')
% ylabel('Function value')
% title([' $\xi=$' num2str(xi), ' $\beta=$' num2str(beta),...
%     '  R=' num2str(R), '  p=' num2str(xmin)...
%     ' fmin=' num2str(fmin) ], 'interpreter','latex')
%     
% %get the directory of your input files:
% [pathstr,namestr,extstr] =...
%     fileparts('/Users/sunyx/Documents/MATLAB/two_species');
% 
% % Create the folder you want the file to be saved to 
% FolderDestination= strcat('xi=', num2str(xi),...
%     ', beta=', num2str(beta),...
%     ', R=', num2str(R),...
%     ', sp=',num2str(sp));
% mkdir ('/Users/sunyx/Desktop',FolderDestination)
% 
% % create new file name for each run
%  Filename=strcat('beta=', num2str(beta),', R=', num2str(R));
%  
% %use that when you save
% matfile = fullfile('/Users/sunyx/Desktop',FolderDestination,...
%     [Filename,'.mat']);
% pdffile = fullfile('/Users/sunyx/Desktop',FolderDestination,...
%     [Filename,'.pdf']);
% figfile = fullfile('/Users/sunyx/Desktop', FolderDestination,...
%     [Filename,'.fig']);
% 
% save(matfile);
% saveas(figure(1),pdffile);
% saveas(figure(1),figfile);

% end

