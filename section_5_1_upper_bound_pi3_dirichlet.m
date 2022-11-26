format long

N=256;
theta=I_pi/3;

% mesh generation
tri_I_intval = [I_intval('0'),I_intval('0');I_intval('1'),I_intval('0');cos(theta),sin(theta)]; % triangular domain
mesh=get_mesh(tri_I_intval,N);

% conforming FEM
[A,B]=create_cg_AB(mesh);
[A,B]=add_dc_cg(A,B,false,mesh);
lambda1=I_veig(sparse(A),sparse(B),1);
upper_bound_lambda1=sup(lambda1(1));

file_name='results/pi3_upper_bound_dirichlet.csv';
writematrix(upper_bound_lambda1,file_name,'WriteMode','append');

disp(upper_bound_lambda1);