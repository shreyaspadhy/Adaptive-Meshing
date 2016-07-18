function [res, coords] = calc_residue(hMesh, phi, idx, kap0, neighbors, mua0, ne, qq, Q)

for i = 1:ne
    coord = (hMesh.Element(i).Data());
    %x_coord = mean(coord(:,1));
    %y_coord = mean(coord(:,2));
    coords(i,:) = mean(coord);
    temp = hMesh.Element(i).ShapeDer(coords(i,:)', 'global');
    grad_el(i,:) = [phi(idx(i,1,qq)), phi(idx(i,2,qq)), phi(idx(i,3,qq))]*temp';
end

grad_el = grad_el./3;

grad_el(:,1) = kap0(:).*grad_el(:,1);
grad_el(:,2) = kap0(:).*grad_el(:,2);

%neighbors = hMesh.ElementNeighbours();
div = zeros(ne,1);

for i = 1:ne
    for j = 1:length(neighbors(i,:))
        if(neighbors(i,j) ~= 0)
            ds = coords(neighbors(i,j),:) - coords(i,:);
            diff = grad_el(neighbors(i,j),:) - grad_el(i,:);
            div(i) = div(i) + dot(diff,ds)./(norm(ds)*norm(ds));
        end
    end
end

second_term = zeros(ne, 1);

for i = 1:ne
    temp2 = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
    second_term(i) = mua0(i).*([phi(idx(i,1,qq)), phi(idx(i,2,qq)), phi(idx(i,3,qq))]*temp2);
end

rhs = zeros(ne,1);

for i = 1:ne
    temp3 = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
    rhs(i) = ([Q(idx(i,1,qq)), Q(idx(i,2,qq)), Q(idx(i,3,qq))]*temp3);
end

res = abs(-div + second_term - rhs);