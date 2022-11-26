function output=u_y2(eigenvec,mesh)
    norm2=0;
    eigenvec=I_intval(eigenvec);
    [elements_n,~]=size(mesh.elements);
    for idx=1:elements_n
        x1=mesh.nodes(mesh.elements(idx,1),1);
        y1=mesh.nodes(mesh.elements(idx,1),2);
        x2=mesh.nodes(mesh.elements(idx,2),1);
        y2=mesh.nodes(mesh.elements(idx,2),2);
        x3=mesh.nodes(mesh.elements(idx,3),1);
        y3=mesh.nodes(mesh.elements(idx,3),2);
        K_measure = abs((x1-x2)*(y3-y2)-(y1-y2)*(x3-x2))/2;

        node_id=zeros(1,3);
        r=I_intval(zeros(1,3));
        for i=1:3
            node_id(i)=mesh.elements(idx,i);
            r(i)=eigenvec(node_id(i));
        end
        
        Ly = 1/(2*K_measure)*[x3-x2;x1-x3;x2-x1];
        local_l2_norm = K_measure*(r*Ly)^2;
        norm2=norm2+local_l2_norm;
    end
    output=norm2;
end