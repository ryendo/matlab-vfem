function recovered_x=x_recover_dc(x,mesh)
    
    x=I_intval(x);
    [nodes_n,~]=size(mesh.nodes);
    recovered_x=I_zeros(nodes_n,1);
    is_node_on_bdry_list = is_node_on_bdry(mesh);
    
    idx=1;
    for i=1:nodes_n
        if is_node_on_bdry_list(i)
            recovered_x(i)=I_intval(0);
        else
            recovered_x(i)=x(idx);
            idx=idx+1;
        end
    end
    
end