function [sNA,sNB]=add_neu_cr(A,B,mesh)
    
    A=sparse(A); B=sparse(B);
    [edges_n,~] = size(mesh.edges);
    [H1,H2,null_dof_idx] = sub_m_constraint(mesh);
    ind=1:edges_n;
    ind(null_dof_idx)=[];
    P=sparse(I_eye(edges_n,edges_n)) - H2*H1;
    NA = P'*A*P; NA = NA(ind,ind); %掛け算に時間かかる
    NB = P'*B*P; NB = NB(ind,ind);
    sNA = I_hull(NA,NA');
    sNB = I_hull(NB,NB');

end

function [OutM1, H, null_dof_idx] = sub_m_constraint(mesh)
    edges_n = size(mesh.edges,1);
    num = edges_n;
    OutM1 = I_intval(zeros(3,num)); 
    OutM2 = I_intval(zeros(3,num)); 
    null_dof_idx = [0,0,0];

    for domain_edge_ind=1:3
        edge_bdry_list = is_edge_on_each_bdry(domain_edge_ind,mesh);
        for k=1:edges_n
            if edge_bdry_list(k)
                g_index = k;
                OutM1( domain_edge_ind, g_index ) = -1;
                if null_dof_idx(domain_edge_ind) == 0
                    OutM2( domain_edge_ind, g_index ) = -1; 
                    null_dof_idx(domain_edge_ind) = g_index;
                end
            end           
            
        end
    end    
    H=OutM2';
end

