format compact short infsup

file_name='results/lower_bounds_dirichlet.csv';

for i=1:170 % 1~170
    tStart = tic;
    lambda1_lower=make_eigen_bounds(theta(i-1),theta(i),N_mesh(i));
    info_result=[i,lambda1_lower]
    writematrix([i,lambda1_lower],file_name,'WriteMode','append');
    tEnd = toc(tStart); one_loop = int16((170-i)*tEnd);
    fprintf('estimated time: %d minutes %d seconds\n',idivide(one_loop,60),mod(one_loop,60));
end

function N_i=N_mesh(i)
    if 0<=i && i<=150
        N_i=40;
    elseif 151<=i && i<=160
        N_i=50;
    elseif 161<=i && i<=170
        N_i=65;
    else
        N_i=3;
    end
end

function theta_i=theta(i)
    if i==0
        theta_i=0;
    elseif 1<=i && i<=45
        theta_i=i*I_intval('0.02')*I_pi/3;
    elseif 46<=i && i<=50
        theta_i=(I_intval('0.90')+(i-45)*I_intval('0.01'))*I_pi/3;
    elseif 51<=i && i<=90
        theta_i=(I_intval('0.95')+(i-50)*I_intval('0.001'))*I_pi/3;
    elseif 91<=i && i<=170
        theta_i=(I_intval('0.99')+(i-90)*I_intval('0.0001'))*I_pi/3;
    else
        theta_i=I_intval(1);
    end
end

function lambda1_lower=make_eigen_bounds(theta_l,theta_u,N)
    
    % mesh generation
    tri_I_intval = [I_intval('0'),I_intval('0');I_intval('1'),I_intval('0');cos(theta_u),sin(theta_u)]; % triangular domain
    cl = I_intval(min(I_inf(cos(theta_l-1)/cos(theta_u-1)),I_inf(cos(theta_l+1)/cos(theta_u+1))));
    mesh=get_mesh(tri_I_intval,5);
    
    % approximate eigenvalue from CR FEM
    [A,B]=create_cr_AB(mesh);
    [A,B]=add_dc_cr(A,B,false,mesh)
    cr_lambda1=I_veig(A,B,1);
    
    % lower bound estimation
    Ch=I_intval('0.1893')/N;
    lower_bound=cr_lambda1(1)/(1+Ch^2*cr_lambda1(1));
    lambda1_lower=I_inf(cl*lower_bound);
end