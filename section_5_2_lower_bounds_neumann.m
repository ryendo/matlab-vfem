format long

file_name='results/lower_bounds_neumann.csv';

for i=1:320 %1~320
    tStart = tic;
    lambda1_lower=make_eigen_bounds(theta(i-1),theta(i),N_mesh(i));
    info_result=[i,lambda1_lower]
    writematrix([i,lambda1_lower],file_name,'WriteMode','append');
    tEnd = toc(tStart); one_loop = int16((320-i)*tEnd);
    fprintf('estimated time: %d minutes %d seconds\n',idivide(one_loop,60),mod(one_loop,60));
end

function N_i=N_mesh(i)
    if 0<=i && i<=312
        N_i=40;
    elseif 313<=i && i<=316
        N_i=50+2*(i-313);
    elseif 317<=i && i<=320
        N_i=56+10*(i-317);
    else
        N_i=3;
    end
end

function theta_i=theta(i)
    if i==0
        theta_i=0;
    elseif 1<=i && i<=40
        theta_i=i*I_intval('0.02')*I_pi/3;
    elseif 41<=i && i<=230
        theta_i=(I_intval('0.80')+(i-40)*I_intval('0.001'))*I_pi/3;
    elseif 231<=i && i<=320
        theta_i=(I_intval('0.99')+(i-230)*I_intval('0.0001'))*I_pi/3;
    else
        theta_i=I_intval(1);
    end
end

function lambda1_lower=make_eigen_bounds(theta_l,theta_u,N)
    
    % mesh generation
    tri_I_intval = [I_intval('0'),I_intval('0');I_intval('1'),I_intval('0');cos(theta_u),sin(theta_u)]; % triangular domain
    cl = I_intval(min(I_inf(cos(theta_l-1)/cos(theta_u-1)),I_inf(cos(theta_l+1)/cos(theta_u+1))));
    mesh=get_mesh(tri_I_intval,N);
    
    % approximate eigenvalue from CR FEM
    [A,B]=create_cr_AB(mesh);
    [A,B]=add_neu_cr(A,B,mesh);

    cr_lambda1=I_veig(sparse(A),sparse(B),1);
    
    % lower bound estimation
    Ch=I_intval('0.1893')/N;
    lower_bound=cr_lambda1(1)/(1+Ch^2*cr_lambda1(1));
    lambda1_lower=I_inf(cl*lower_bound);
end