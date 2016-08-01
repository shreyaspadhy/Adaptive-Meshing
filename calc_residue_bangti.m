function [res] = calc_residue_bangti(hMesh, phi, idx, kap0, neighbors, mua0, mua_node, ne, Q, vert)

jump_total = zeros(ne,1);

phi = sum(phi, 2)./100;
Q = sum(Q, 2)./100;

figure(5); hMesh.Display(full(phi));
figure(6); hMesh.Display(kap0);

grad_mua = zeros(ne,2);
for i = 1:ne
    temp_el = hMesh.Element(i).ShapeDer(mean(hMesh.Element(i).Data())', 'global');
    grad_mua(i,:) = [mua_node(idx(i,1)), mua_node(idx(i,2)), mua_node(idx(i,3))]*temp_el';
    for j = 1:length(neighbors(i,:))
        if(neighbors(i,j) ~= 0)
            
            vertices = intersect(idx(i,:), idx(neighbors(i,j),:));
            
            vert_mid = (vert(vertices(1),:) + vert(vertices(2),:))./2;
            
            temp_el = hMesh.Element(i).ShapeDer(vert_mid', 'global');
            temp_nbr = hMesh.Element(neighbors(i,j)).ShapeDer(vert_mid', 'global');
            
            grad_mid_el = [phi(idx(i,1)), phi(idx(i,2)), phi(idx(i,3))]*temp_el';
            grad_mid_nbr = [phi(idx(neighbors(i,j),1)), phi(idx(neighbors(i,j),2)), phi(idx(neighbors(i,j),3))]*temp_nbr';
            
            
            jump = abs(norm(kap0(i).*grad_mid_el - kap0(neighbors(i,j)).*grad_mid_nbr));
            
            side_length = sqrt((vert(vertices(2),2) - vert(vertices(1),2))^2 + (vert(vertices(2),1) - vert(vertices(1),1))^2);
            
            jump = jump.*jump.*side_length;
            jump_total(i) = jump_total(i) + jump;
        end
    end
    
end
abs_grad_mua = sqrt(grad_mua(:,1).*grad_mua(:,1) + grad_mua(:,2).*grad_mua(:,2))
figure(7); hMesh.Display(abs_grad_mua);

for i = 1:ne
    coord = (hMesh.Element(i).Data());
    coords(i,:) = mean(coord);
    temp = hMesh.Element(i).ShapeDer(coords(i,:)', 'global');
    grad_el(i,:) = [phi(idx(i,1)), phi(idx(i,2)), phi(idx(i,3))]*temp';
end
% 
grad_el = grad_el./3;

grad_el(:,1) = kap0(:).*grad_el(:,1);
grad_el(:,2) = kap0(:).*grad_el(:,2);
% 
% %neighbors = hMesh.ElementNeighbours();
div = zeros(ne,1);

hT = zeros(ne,1);
% 
for i = 1:ne
    for j = 1:length(neighbors(i,:))
        if(neighbors(i,j) ~= 0)
            ds = coords(neighbors(i,j),:) - coords(i,:);
            diff = grad_el(neighbors(i,j),:) - grad_el(i,:);
            div(i) = div(i) + dot(diff,ds)./(norm(ds)*norm(ds));
            hT(i) = sqrt(hMesh.Element(i).Size());
        end
    end
end

rhs = zeros(ne,1);

for i = 1:ne
    temp3 = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
    rhs(i) = ([Q(idx(i,1)), Q(idx(i,2)), Q(idx(i,3))]*temp3);
end

second_term = zeros(ne, 1);

for i = 1:ne
    temp2 = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
    second_term(i) = mua0(i).*([phi(idx(i,1)), phi(idx(i,2)), phi(idx(i,3))]*temp2);
end

res1 = abs(-div + second_term - rhs);

sizes = hMesh.ElementSize();
res1 = res1.*res1.*sizes;
%res = res1 + jump_total;
%res = jump_total;

res = res1./max(res1) + jump_total./max(jump_total) + abs_grad_mua./max(abs_grad_mua);
%res = res1./max(res1) + jump_total./max(jump_total);
%res = jump_total;
%+ jump_total;
% 

% 
% rhs = zeros(ne,1);
% 
% for i = 1:ne
%     temp3 = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
%     rhs(i) = ([Q(idx(i,1,qq)), Q(idx(i,2,qq)), Q(idx(i,3,qq))]*temp3);
% end
% 
% res = abs(-div + second_term - rhs);
%     temp = hMesh.Element(i).ShapeDer(coords(i,:)', 'global');
%     grad_el(i,:) = [phi(idx(i,1,qq)), phi(idx(i,2,qq)), phi(idx(i,3,qq))]*temp';
% end
% 
% grad_el = grad_el./3;
% 
% grad_el(:,1) = kap0(:).*grad_el(:,1);
% grad_el(:,2) = kap0(:).*grad_el(:,2);
% 
% %neighbors = hMesh.ElementNeighbours();
% div = zeros(ne,1);
% 
% for i = 1:ne
%     for j = 1:length(neighbors(i,:))
%         if(neighbors(i,j) ~= 0)
%             ds = coords(neighbors(i,j),:) - coords(i,:);
%             diff = grad_el(neighbors(i,j),:) - grad_el(i,:);
%             div(i) = div(i) + dot(diff,ds)./(norm(ds)*norm(ds));
%         end
%     end
% end
% 
% second_term = zeros(ne, 1);
% 
% for i = 1:ne
%     temp2 = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
%     second_term(i) = mua0(i).*([phi(idx(i,1,qq)), phi(idx(i,2,qq)), phi(idx(i,3,qq))]*temp2);
% end
% 

% 
% res = abs(-div + second_term - rhs);