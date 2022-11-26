format compact short infsup

N=128;
theta=I_pi/I_intval('3');

% mesh generation
tri_I_intval = [I_intval('0'),I_intval('0');I_intval('1'),I_intval('0');cos(theta),sin(theta)]; % triangular domain
mesh=get_mesh(tri_I_intval,N);

% conforming FEM
[A_cg,B_cg]=create_cg_AB(mesh);
[A_cg_neu,B_cg_neu,~,~]=add_neu_cg(A_cg,B_cg,mesh);
lambda1=I_veig(sparse(A_cg_neu),sparse(B_cg_neu),1);
upper_bound_lambda1=sup(lambda1(1));

file_name='results/pi3_upper_bound_neumann.csv';
writematrix(upper_bound_lambda1,file_name,'WriteMode','append');

disp(upper_bound_lambda1);