function [res, res_nh] = calc_vol_residue_3D(hMesh, phi, idx, mua0, mus0, Q)
%Calculates an error residue at the center of each element
% Inputs
     %hMesh: a toastMesh object
     %phi  : a n_nodes x n_sources matrix containing phi values per source
     %idx  : a n_el x 5 matrix containing element data (4 vertices+1 level)
     %mua0 : a n_nodes x 1 matrix containing mua0
     %mus0 : a n_nodes x 1 matrix containing mus0
     %Q    : Qvec, a n_nodes x n_sources matrix containing source data

% Output
     %res  : a n_el x 1 matrix containing error residue at each el center
     %res_nh : a n_vodes x 1 matrix containing error residue at each node
     
% TO RUN: 
     %[res, res_nh] = calc_vol_residue_3D(toast_mesh, ff.exPhi, ff.element, ff.optiProp.o_mua, ff.optiProp.o_mus, Ex_qvec);
neighbors = hMesh.ElementNeighbours();
ne = length(neighbors);
nh = length(Q);
coords = zeros(ne,3);

phi = mean(phi,2);
Q = mean(Q,2);

for i = 1:nh
    kap0(i) = 1./3.*(mua0(i)+mus0(i));
end

for i = 1:ne
    coord = (hMesh.Element(i).Data())
    %x_coord = mean(coord(:,1));
    %y_coord = mean(coord(:,2));
    coords(i,:) = mean(coord);
    temp = hMesh.Element(i).ShapeDer(coords(i,:)', 'global');
    grad_el(i,:) = [kap0(idx(i,1)).*phi(idx(i,1)), kap0(idx(i,2)).*phi(idx(i,2)), kap0(idx(i,3)).*phi(idx(i,3)), kap0(idx(i,4)).*phi(idx(i,4))]*temp';
end

% grad_el(:,1) = kap0(:).*grad_el(:,1);
% grad_el(:,2) = kap0(:).*grad_el(:,2);

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
    second_term(i) = ([mua0(idx(i,1)).*phi(idx(i,1)), mua0(idx(i,2)).*phi(idx(i,2)), mua0(idx(i,3)).*phi(idx(i,3)), mua0(idx(i,4)).*phi(idx(i,4))]*temp2);
end

rhs = zeros(ne,1);

for i = 1:ne
    temp3 = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
    rhs(i) = ([Q(idx(i,1)), Q(idx(i,2)), Q(idx(i,3)), Q(idx(i,4))]*temp3);
end

res = abs(-div + second_term - rhs);

sizes = hMesh.ElementSize();

for i = 1:ne
    res(i) = res(i).*sizes(i);
end

res_nh = zeros(nh,1);
for i = 1:nh
    I = find(idx == i)
    for j = 1:length(I)
        a = mod(I(j),ne);
        if(a == 0)
            a = ne;
        end
        res_nh(i) = res_nh(i) + res(a);
    end
    res_nh(i) = res_nh(i)./length(I);
end