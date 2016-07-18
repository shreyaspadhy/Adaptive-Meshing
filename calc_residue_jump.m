function [res] = calc_residue_jump(hMesh, phi, idx, kap0, neighbors, mua0, ne, qq, Q, vert)

jump_total = zeros(ne,1);

% for i = 1:ne
%     coord = (hMesh.Element(i).Data());
%     coords(i,:) = mean(coord);
% end


for i = 1:ne
    for j = 1:length(neighbors(i,:))
        if(neighbors(i,j) ~= 0)
            
            vertices = intersect(idx(i,:,qq), idx(neighbors(i,j),:,qq));
            
            vert_mid = (vert(vertices(1),:) + vert(vertices(2),:))./2;
            
            temp_el = hMesh.Element(i).ShapeDer(vert_mid', 'global');
            temp_nbr = hMesh.Element(neighbors(i,j)).ShapeDer(vert_mid', 'global');
            
            grad_mid_el = [phi(idx(i,1,qq)), phi(idx(i,2,qq)), phi(idx(i,3,qq))]*temp_el';
            grad_mid_nbr = [phi(idx(neighbors(i,j),1,qq)), phi(idx(neighbors(i,j),2,qq)), phi(idx(neighbors(i,j),3,qq))]*temp_nbr';
            
            jump = abs(norm(kap0(i).*grad_mid_el - kap0(neighbors(i,j)).*grad_mid_nbr));
            
            side_length = sqrt((vert(vertices(2),2) - vert(vertices(1),2))^2 + (vert(vertices(2),1) - vert(vertices(1),1))^2);
            
            jump = jump.*side_length;
            jump_total(i) = jump_total(i) + jump;
        end
    end
    
end

res = jump_total;

for i = 1:ne
    coord = (hMesh.Element(i).Data());
    coords(i,:) = mean(coord);
    temp = hMesh.Element(i).ShapeDer(coords(i,:)', 'global');
    grad_el(i,:) = [phi(idx(i,1,qq)), phi(idx(i,2,qq)), phi(idx(i,3,qq))]*temp';
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
    rhs(i) = ([Q(idx(i,1,qq)), Q(idx(i,2,qq)), Q(idx(i,3,qq))]*temp3);
end

second_term = zeros(ne, 1);

for i = 1:ne
    temp2 = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
    second_term(i) = mua0(i).*([phi(idx(i,1,qq)), phi(idx(i,2,qq)), phi(idx(i,3,qq))]*temp2);
end

res1 = abs(-div + second_term - rhs);

res = res1 + jump_total;

sizes = hMesh.ElementSize();

res = res.*sizes;
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