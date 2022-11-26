function [A,B]=create_cg_AB(mesh)

    %create A and B
    [nodes_n,~] = size(mesh.nodes);
    [elements_n,~] = size(mesh.elements);

    A=I_intval(I_zeros(nodes_n,nodes_n));
    B=I_intval(I_zeros(nodes_n,nodes_n));
    for idx=1:elements_n
        mesh.edges = mesh.nodes( mesh.elements(idx,[3,1,2]), : ) - mesh.nodes( mesh.elements(idx,[2,3,1]), :); %3 by 2
        S =  abs(mesh.edges(1,:)/2*[mesh.edges(2,2); -mesh.edges(2,1)]);
        e_e = mesh.edges*mesh.edges';

        global_idx = mesh.elements(idx,1:3);

        % contribute local matirix to A and B
        A(global_idx, global_idx) = A(global_idx, global_idx) + e_e/(4*S);
        B(global_idx, global_idx) = B(global_idx, global_idx) + S/I_intval(12)*(I_ones(3,3)+I_eye(3,3));
        
    end
        
end