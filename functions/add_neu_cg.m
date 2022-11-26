function [sNA,sNB,C,deleted_basis]=add_neu_cg(A,B,mesh)
    A=sparse(A); B=sparse(B);
    [nodes_n,~] = size(mesh.nodes);
    F1 = is_node_on_each_bdry(1,mesh); F2 = is_node_on_each_bdry(2,mesh); F3 = is_node_on_each_bdry(3,mesh);
    deleted_basis = sort(find(F1+F2+F3==2,3));
    C = I_intval([F1;F2;F3]); C(:,deleted_basis)=C(:,deleted_basis)/2; C=C*2;
    
    %P:null basis of C
    sym_C=(C(1,:)+C(2,:)+C(3,:))/2;
    C(1,:)=sym_C-C(1,:); C(2,:)=sym_C-C(2,:); C(3,:)=sym_C-C(3,:);
    P=I_eye(nodes_n,nodes_n);
    P(deleted_basis(1),:)=-C(1,:); P(deleted_basis(1),deleted_basis(1))=0;
    P(deleted_basis(2),:)=-C(2,:); P(deleted_basis(2),deleted_basis(2))=0;
    P(deleted_basis(3),:)=-C(3,:); P(deleted_basis(3),deleted_basis(3))=0;
    P(:,deleted_basis)=[];

    NA = P'*A*P; NB = P'*B*P;
    
    sNA = I_hull(NA,NA');
    sNB = I_hull(NB,NB');
    
end
