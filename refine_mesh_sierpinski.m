function [h] = refine_mesh_sierpinski(vert,idx,re, mua0, neighbors)
% refine the set 're' in list 'el'
% linear triangles only, so far

nh = length(vert);
ne = length(idx);

flag = zeros(length(re),2);

for i = 1:length(re)
    flag(i,:) = [re(i), 0];
end

for i = 1:length(re)
    if(flag(i,2) ~= 1)
        mid_vert = zeros(3,2);
        mid_vert_ind = zeros(3,1);
        
        elem = idx(re(i),:);
        elem_ind = re(i)
        
        %Add the three new vertices
        for j = 1:3
            
            vert1 = elem(mod(j-1,3)+1);
            vert2 = elem(mod(j,3)+1);
            vert3 = elem(mod(j+1,3)+1);
            mid_vert(j,:) = (vert(vert1,:) + vert(vert2,:))/2;
            
            vert(nh+1,:) = mid_vert(j,:);
            mid_vert_ind(j) = nh+1;
            nh = nh+1;
            flag1 = 0;
            
            for k = 1:3
                if(neighbors(elem_ind,k) ~= 0 && flag1 ~= 1)
                    idx(neighbors(elem_ind,k),:)
                    elem
                    com_vert_temp = intersect(elem, idx(neighbors(elem_ind, k),:))
                    
                    idx(neighbors(elem_ind,k),:)
                    I = find(idx(neighbors(elem_ind,k),:) == vert1);
                    J = find(idx(neighbors(elem_ind,k),:) == vert2);
                    flag1 = 0;
                    
                    if(~isempty(I) && ~isempty(J))
                        com_vert = [idx(neighbors(elem_ind,k),I), idx(neighbors(elem_ind,k),J)];
                        opp_vert = [idx(neighbors(elem_ind,k),6 - I - J)];
                        index = k;
                        flag1 = 1;
                        no_add = 0;
%                         for l = 1:3
%                             if(idx(neighbors(elem_ind,k),l) ~= com_vert(1) && idx(neighbors(elem_ind,k),l) ~= com_vert(2))
%                                 opp_vert = idx(neighbors(elem_ind, k), l, qq);
%                             end
%                         end
                    elseif(flag1 ~= 1)
                            no_add = 1;
                    end
                end
            end
            
            if(no_add ~= 1)
            idx(neighbors(elem_ind,index),:) = [com_vert(1); mid_vert_ind(j); opp_vert];
            
            idx(ne+1,:) = [com_vert(2); mid_vert_ind(j); opp_vert];
            ne = ne+1;
            
            remove_ind = find(re == neighbors(elem_ind,index))
            if(remove_ind ~= 0)
                flag(remove_ind,2) = 1;
            end
            else
                no_add = 0;
            end
            
        end
        
        idx(ne+1,:) = [idx(re(i),1); mid_vert_ind(1); mid_vert_ind(3)];
        
        idx(ne+2,:) = [idx(re(i),2); mid_vert_ind(1); mid_vert_ind(2)];
        
        idx(ne+3,:) = [idx(re(i),3); mid_vert_ind(2); mid_vert_ind(3)];
        
        ne = ne+3;
        
        idx(re(i),:) = [mid_vert_ind(1); mid_vert_ind(2); mid_vert_ind(3)];
        
        h = toastMesh(vert,idx,15*ones(ne,1));
        
        neighbors = h.ElementNeighbours();
    end
end
        
sizes = h.ElementSize();
[vertnew, idxnew] = h.Data();

for i = 1:ne
    if(sizes(i) <= 0)
        %sizes(i)
        temp = idxnew(i, :);
        idxnew(i,:) = [temp(2), temp(1), temp(3)];
        %h = toastMesh(vertnew,idxnew,15*ones(ne,1));
        %sizes = h.ElementSize();
        %sizes(i)
    end
end

h = toastMesh(vertnew,idxnew,15*ones(ne,1));

    