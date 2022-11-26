format long
N=64;

file_name='results/simplicity_neu.csv';

theta_l=I_intval('999')/I_intval('1000')*I_pi/3;
theta_u=I_pi/3;
theta=I_hull(theta_l,theta_u);
[lambda1_upper,lambda2_lower]=make_eigen_bounds(theta, N);
result=[lambda1_upper,lambda2_lower];
info_result=[lambda1_upper,lambda2_lower]
writematrix(result,file_name,'WriteMode','append');

function [upper_bound_lambda1_I,lower_bound_lambda2_I]=make_eigen_bounds(theta, N)
    theta_l=I_intval(I_inf(theta));
    theta_u=I_intval(I_sup(theta));
    
    % mesh generation
    tri_I_intval = [I_intval('0'),I_intval('0');I_intval('1'),I_intval('0');cos(theta_l),sin(theta_l)]; % triangular domain
    mesh=get_mesh(tri_I_intval,N);
    
    % conforming FEM
    [A,B]=create_cg_AB(mesh);
    [A,B]=add_neu_cg(A,B,mesh);
    lambda1=I_veig(sparse(A),sparse(B),1);
    upper_bound_lambda1=sup(lambda1(1));

    % upper bound estimation
    upper_bound_lambda1_I=I_sup(cos(theta_l)^2/cos(theta_u)^2*upper_bound_lambda1);
    
    % CR FEM
    [A,B]=create_cr_AB(mesh);
    [A,B]=add_neu_cr(A,B,mesh);
    cr_lambdas2=I_veig(sparse(A),sparse(B),1:2);
    cr_lambda2=cr_lambdas2(2);

    % lower bound estimation
    Ch=I_intval('0.1893')/N;
    lower_bound_lambda2=cr_lambda2(1)/(1+Ch^2*cr_lambda2(1));
    lower_bound_lambda2_I=I_inf(sin(theta_l)^2/sin(theta_u)^2*lower_bound_lambda2);
    
end