function recovered_x=x_recover_neu(x,C,deleted_basis)
    
    x=I_intval(x);
    [~,nodes_n]=size(C);
    
    idx=1;
    for i=1:nodes_n
        if ismember(i,deleted_basis)
            recovered_x(i)=I_intval(0);
        else
            recovered_x(i)=x(idx);
            idx=idx+1;
        end
    end
    
    recovered_x(deleted_basis(1))=-C(1,:)*recovered_x';
    recovered_x(deleted_basis(2))=-C(2,:)*recovered_x';
    recovered_x(deleted_basis(3))=-C(3,:)*recovered_x';

end