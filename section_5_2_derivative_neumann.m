%% implementation of Algorithm 2 in section 5.2.
format compact short infsup

N_veig = 64; % default value: 64
N_eigs = 96; % default value: 96


for i= 1:100
    tStart = tic;

    I=I_hull(theta(i),theta(i-1));
    [F_lower,F_upper,Err]=derivative_calc(I,N_veig,N_eigs);
    
    fprintf('N_veig=%d, N_eig=%d\n',N_veig,N_eigs);
    fprintf('i=%d, F_lower=%.8f, F_upper=%.8f, Err=%.8f\n\n',i,F_lower,F_upper,Err);
    
    file_name='results/derivative_neumann.csv';
    writematrix([i,F_lower,F_upper,Err],file_name,'WriteMode','append');
    
    tEnd = toc(tStart); one_loop = int32((100-i)*tEnd);
    fprintf('estimated time: %d minutes %d seconds\n',idivide(one_loop,60),mod(one_loop,60));
end

function theta_i = theta(i)
    theta_i = I_pi/3*(1-i/I_intval('10')^5);
end

function [F_inf,F_sup,Err]=derivative_calc(I,N_veig,N_eigs)
    
    Ch=I_intval('0.1893')/N_veig;

    I0=I_intval(I_inf(I));
    I1=I_intval(I_sup(I));
    c01= I_intval(min(I_inf(cos(I1-1)/cos(I0-1)),I_inf(cos(I1+1)/cos(I0+1))));

%%  create meshes and matrices
%   create mesh for rigorous calculation
    tri_I_intval = [I_intval('0'),I_intval('0');I_intval('1'),I_intval('0');cos(I0),sin(I0)]; % triangular domain
    mesh_veig=get_mesh(tri_I_intval,N_veig);
    
    [A_cr,B_cr]=create_cr_AB(mesh_veig);
    [A_cr_neu,B_cr_neu]=add_neu_cr(A_cr,B_cr,mesh_veig);

%   create mesh for approximate calculation
    mesh_eigs=get_mesh(tri_I_intval,N_eigs);
    [A_cg_eigs,B_cg_eigs]=create_cg_AB(mesh_eigs);
    [A_cg_neu_eigs,B_cg_neu_eigs,C,deleted_basis]=add_neu_cg(A_cg_eigs,B_cg_eigs,mesh_eigs);
    fA_neu=I_mid(A_cg_neu_eigs); fB_neu=I_mid(B_cg_neu_eigs);
    
%%  rigorously estimate bounds of eigenvalues
    
%   lower bound for lambda1
    cr_lambdas=I_veig(sparse(A_cr_neu),sparse(B_cr_neu),1:2);
    l_lam1_0=cr_lambdas(1)/(1+Ch^2*cr_lambdas(1));
    l_lam1=I_intval(I_inf(c01*l_lam1_0));
    
%   lower bound for lambda2
    l_lam2_0=cr_lambdas(2)/(1+Ch^2*cr_lambdas(2));
    l_lam2=I_intval(I_inf(c01*l_lam2_0));

%%  calculate an approximate eigenfunction
    [v,~]=eigs(sparse(fA_neu),sparse(fB_neu),1,'sm');
    x_hat_0=x_recover_neu(v,C,deleted_basis);
    norm_x_hat_0=x_hat_0*B_cg_eigs*x_hat_0';

%%  calculate approximate derivative
    a = (cos(I)-cos(I0))/sin(sin(I0));
    b = sin(I)/sin(I0);
    xx_0=u_x2(x_hat_0,mesh_eigs); xy_0=u_xu_y(x_hat_0,mesh_eigs); yy_0=u_y2(x_hat_0,mesh_eigs);
    xx=b*xx_0; xy=-a*xx_0+xy_0; yy=a^2*b^(-1)*xx_0-2*a*b^(-1)*xy_0+b^(-1)*yy_0; norm_x_hat=b*norm_x_hat_0;

    F_0=(-2/tan(I0)*yy_0+2*xy_0)/norm_x_hat_0;
    F=(-2/tan(I)*(yy)+2*xy)/norm_x_hat;

%%  calculate approximate eigenvalue
    lam1h_0=I_intval((x_hat_0*A_cg_eigs*x_hat_0')/(x_hat_0*B_cg_eigs*x_hat_0'));
    lam1h=c01*lam1h_0;
    u_lam1h=I_intval(I_sup(lam1h));

%%  calculate error
    rho=l_lam2;
    
    eta = sqrt(l_lam1+u_lam1h-2*l_lam1*sqrt((rho-u_lam1h)/(rho-l_lam1)));
    Err=2*(2*sqrt(u_lam1h)/tan(I)+sqrt(2*u_lam1h))*eta;
    F = hull(F-Err,F+Err);
    
%%  I_intval to double
    F_inf=I_inf(F); F_sup = I_sup(F); Err=I_sup(Err);

end