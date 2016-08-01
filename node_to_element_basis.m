function [mua, mus] = node_to_element_basis(hMesh, idx, mua_node, mus_node)

for i = 1:hMesh.ElementCount()
    coord = (hMesh.Element(i).Data());
    %x_coord = mean(coord(:,1));
    %y_coord = mean(coord(:,2));
    coords(i,:) = mean(coord);
    temp = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
    mua(i) = ([mua_node(idx(i,1)), mua_node(idx(i,2)), mua_node(idx(i,3))]*temp);
    mus(i) = ([mus_node(idx(i,1)), mus_node(idx(i,2)), mus_node(idx(i,3))]*temp);
end
