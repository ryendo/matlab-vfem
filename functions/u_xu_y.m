function output=u_xu_y(eigenvec,mesh)
    dot=0;
    eigenvec=I_intval(eigenvec);
    [elements_n,~]=size(mesh.elements);
    for idx=1:elements_n
%         element_ijk=mesh.elements(idx,:);
%         xy_list=mesh.nodes(element_ijk,:);
%       おなじコードをできるだけまとめる

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

% 行列作成のコードを参考にして、edge の形でかく。（各点の座標をつかわない）        
        Lx = 1/(2*K_measure)*[y2-y3;y3-y1;y1-y2];
        Ly = 1/(2*K_measure)*[x3-x2;x1-x3;x2-x1];
        local_l2_dot = K_measure*(r*Lx)*(r*Ly);
        dot=dot+local_l2_dot;
    end
    output=dot;
end