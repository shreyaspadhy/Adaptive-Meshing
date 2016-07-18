function [h, mua] = refine_mesh_sierpinski(vert,el,re, mua0, mua, neighbors)
% refine the set 're' in list 'el'
% linear triangles only, so far

nh = length(vert);
ne = length(el);
mua = mua0;

for  i = 1:length(re)
    %if(re(j) ~= 0)
       % A = find(neighbors(re(j),:) == 0);
%     end
    if(re(j) ~= 0)
        ee = el(re(j),:);
        nbrs = neighbors(re(j),:);
        
        vnew = zeros(3,2);
        vind = zeros(3,1);
        flag = 0;
        flag2 = 0;
        
        for j = 1:3
            if(nbrs(j) ~= 0)
                cmn_vert = intersect(ee, idx(nbrs(j),:,qq));
                if(length(cmn_vert) ~= 2)
                    display('Neighbor does not share two common vertices');
                else
                    mid = (vert(cmn_vert(1),:)+vert(cmn_vert(2),:))/2;
        
        for i = 1:3
            vert1 = ee(mod(i,3)+1)
            vert2 = ee(mod(i+1,3)+1)
            
            
            for k = 1:length(nbrs)
                if(nbrs(i) ~= 0)
                    if(nbrs(k) ~= 0)
                    nbrs(k)
                    el(nbrs(k))
                    I = find(el(nbrs(k),:) == vert1)
                    J = find(el(nbrs(k),:) == vert2)
                    if(~isempty(I) && ~isempty(J))
                        el_nbr = nbrs(k);
                        vert3 = el(el_nbr, 6 - I - J);
                        flag2 = 1;
                    end
                    end
                else
                    flag = 1;
                end
            end
        
        vnew(i,:) = (vert(vert1,:)+vert(vert2,:))/2;
        vert(nh+1,:) = vnew(i,:);
        nh = nh+1;
        vind(i) = nh;
        
        if(flag ~= 1 && flag2 == 1)
            flag2 = 0;
            el(ne,:)
            [vert1, vert3, nh]
            el(ne+1,:) = [vert1, vert3, vind(i)];
            mua = [mua; mua0(el_nbr)];
            el(ne+2,:) = [vert2, vert3, vind(i)];
            mua = [mua; mua0(el_nbr)];
            ne = ne+2;
%             for l = 1:length(nbrs)
%                 K = find(re == nbrs(l))
%                 re(K) = 0;
%             end
        else
            flag = 0;
        end
        end
        %h = toastMesh(vert,el,15*ones(ne,1));
        %neighbors = h.ElementNeighbours();
    
    
    el(ne+1,:) = [vind(1), vind(2), vind(3)];
    mua = [mua; mua0(re(j))];
    ne = ne+1;
    el(ne+1,:) = [vind(1), vind(2), ee(2)];
    mua = [mua; mua0(re(j))];
    ne = ne+1;
    el(ne+1,:) = [vind(3), vind(2), ee(3)];
    mua = [mua; mua0(re(j))];
    ne = ne+1;
    el(ne+1,:) = [vind(1), vind(3), ee(1)];
    mua = [mua; mua0(re(j))];
    ne = ne+1;
    
    for l = 1:length(nbrs)
        if(nbrs(l) ~= 0)
         K = find(re == nbrs(l))
         re(K) = 0;
        end
    end
    end
end
h = toastMesh(vert,el,15*ones(ne,1));