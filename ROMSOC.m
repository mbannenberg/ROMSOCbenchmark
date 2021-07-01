% ------------------------------------------------------------------------------
% Reduced Order Multirate Simulation Of Circuits BDF-1
% 
% 
% This script parses a netlist *.cir file into matrices and function
% vectors and casts them in the form:
%
%       F = @(x_p,x,t) A*x + E*x_p + func_p(x) + func_r(t);
%
%
%
% Copyright 2021 Marcus Bannenberg (BUW, bannenberg@uni-wuppertal.de)
% ------------------------------------------------------------------------------

clear all

% Add parser m_files
addpath(genpath('Parser'));
% Add components m_files
addpath(genpath('Components'));
% Add numerical m_files
addpath(genpath('Numerical'));
% Add circuits
addpath(genpath('Circuits'));


% Select circuit
fname = 'multirate_example_long_nonlin.cir';

% Parse the netlist to matrices and function vectors
[E, A, func_p, func_r, x0] = parse_netlist(fname);

% Set simulation parameters
t_0 = 0; t_end = 0.004;
N = 100; m = 20;
t = linspace(t_0,t_end,m*N+1);
tol = 1e-8;

M = numel(x0);
idx_A = [1:3 M]';
idx_L = [4:M-1]';

% Compute reference solution
try
    B = strcat(regexprep(cellfun(@func2str, func_p, 'uni', 0), '^@\(x\)', ''), ';');
    func_p_ode =  str2func(strcat('@(x) [', B{:}, ']'));
    
    B = strcat(regexprep(cellfun(@func2str, func_r, 'uni', 0), '^@\(t\)', ''), ';');
    func_r_ode =  str2func(strcat('@(t) [', B{:}, ']'));

    f = @(t,x) -A*x - func_p_ode(x) - func_r_ode(t);
    opts = odeset('Mass',E,'RelTol',1e-13);
    disp('Solving the circuit with ode15s');
    tic
    [t_ref,y_ref] = ode15s(f,t,x0,opts);
    toc
    y_ref = y_ref';
    
    
    Y = zeros(size(y_ref));
    for i = 1:numel(t)
        Y(:,i) = func_p_ode(y_ref(:,i));
    end


    X = Y;
    X_diff = X(:,2:end)-X(:,1:end-1); 
    epsilon = epsilon_procedure(X_diff);
    mor_object.U_gap = mess(X_diff,0.1);
    mor_object.g = size(mor_object.U_gap,2);

    [mor_object.S, mor_object.M, mor_object.U_gap_S, mor_object.idx] = gappypod(mor_object.U_gap);
    
    X = Y;
    X_diff = y_ref(idx_L,2:end) - y_ref(idx_L,1:end-1);
    epsilon = epsilon_procedure(X_diff);
    u = mess(X_diff,epsilon);
    r_mess = size(u,2);                                                     
    mor_object.r = r_mess;
    mor_object.dim_org = size(y_ref,2); 
    mor_object.U_r = u;

    mor_object.Phi = mor_object.U_r.';
    mor_object.dim_A = numel(idx_A); mor_object.dim_L = numel(idx_L);

    Phi_full = zeros(mor_object.r+numel(idx_A),numel(idx_A)+numel(idx_L));
    Phi_full(idx_A(idx_A < idx_L(1)),idx_A(idx_A < idx_L(1))) = eye(size(idx_A(idx_A<idx_L(1)),1));
    Phi_full(idx_L(1):idx_L(1) + mor_object.r-1,idx_L(1):idx_L(1) + numel(idx_L)-1) = mor_object.Phi;
    
    % second block starting point
    sbsp = numel(idx_A(idx_A < idx_L(1))) + mor_object.r + 1;
    Phi_full(sbsp:end,idx_A(idx_A>idx_L(end))) = eye(numel(idx_A(idx_A>idx_L(end))));
    mor_object.Phi = Phi_full;
    U_r_full = zeros(numel(idx_A)+numel(idx_L),mor_object.r+numel(idx_A));
    U_r_full(idx_A(idx_A < idx_L(1)),idx_A(idx_A < idx_L(1))) = eye(size(idx_A(idx_A<idx_L(1)),1));
    U_r_full(idx_L(1):idx_L(1) + numel(idx_L)-1,idx_L(1):idx_L(1) + mor_object.r-1) = mor_object.U_r;
    
    % second block starting point
    sbsp = numel(idx_A(idx_A < idx_L(1))) + mor_object.r + 1;
    U_r_full(idx_A(idx_A>idx_L(end)),sbsp:end) = eye(numel(idx_A(idx_A>idx_L(end))));
    mor_object.U_r = U_r_full;    
catch
    disp('Reference solver failed');
end

% Array of the number of time steps mulitplied by a factor 2
N_list = N*2.^[0:3];
for i = 1:numel(N_list)
    N = N_list(i);
    t_mr = linspace(t_0,t_end,m*N+1);
    t = linspace(t_0,t_end,N+1);
    tol = 1e-8;
    disp('Solving the circuit with ROMRBDF-1');
    tic
    [t_mr, y_romrbdf] = ROMRBDF(E, A, func_p, func_r, x0,t_mr,1,tol,m,mor_object);
    romr_timer(i) = toc
    
    disp('Solving the circuit with MRBDF-1');
    tic
    [t_mr, y_mrbdf] = MRBDF(E, A, func_p, func_r, x0,t_mr,1,tol,m);
    mr_timer(i) = toc
    

    disp('Solving the circuit with BDF-1');
    tic
    [t_bdf, y_bdf] = BDF(E, A, func_p, func_r, x0,t,1,tol);
    timer(i) = toc
    

    error_ROMRBDF(:,i) = abs(y_romrbdf(:,end) - y_ref(:,end));
    error_MRBDF(:,i) = abs(y_mrbdf(:,end) - y_ref(:,end));
    error_BDF(:,i) = abs(y_bdf(:,end) - y_ref(:,end));
end

% Set colors  for plotting
color_blue = [0, 0.4470, 0.7410];
color_red = [0.8500, 0.3250, 0.0980];
color_yellow = [0.9290, 0.6940, 0.1250];

% Create figure
h = figure();
loglog(N_list,error_BDF(3,:)','o-','LineWidth',2,'Color',color_blue,'MarkerSize',8); hold on;
loglog(N_list,(error_BDF(3,1)./(2.^[1:numel(N_list)]))','LineWidth',2,'color','black');
loglog(N_list,error_ROMRBDF(3,:)','+-','LineWidth',2,'Color',color_red,'MarkerSize',8); 
legend('BDF','O(H)','ROMR');
grid on;
title('Convergence of the error');
xlabel('Number of macro steps');
ylabel('|u^{ref}_{3} - u_3|');
axis([N_list(1)*0.9 N_list(end)*1.1 error_ROMRBDF(3,end)*0.5 error_BDF(3,1)*1.5]);
set(gca, 'FontName', 'Times New Roman','FontSize',14);

figure()
loglog(timer,error_BDF(3,:)','o-','LineWidth',2,'Color',color_blue,'MarkerSize',8); hold on;
loglog(1e2*(2.^[0:numel(timer)-1]),(6e-5./(2.^[0:numel(timer)-1]))','LineWidth',2,'color','black');
loglog(mr_timer,error_MRBDF(3,:)','o-','LineWidth',2,'Color',color_yellow,'MarkerSize',8); 
loglog(romr_timer,error_ROMRBDF(3,:)','+-','LineWidth',2,'Color',color_red,'MarkerSize',8);
grid on;
title('Computation time vs the error');
xlabel('Computation time in seconds');
ylabel('|u^{ref}_{3} - u_3|');
axis([romr_timer(1)*0.9 mr_timer(end)*1.1 error_ROMRBDF(3,end)*0.5 error_BDF(3,1)*1.5]);
set(gca, 'FontName', 'Times New Roman','FontSize',14);
legend('BDF','O(H)','MR','ROMR');
