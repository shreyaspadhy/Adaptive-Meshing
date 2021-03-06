function [h] = refine_mesh_sierpinski_3D(vert,idx,re, neighbors)
% refine the set 're' in list 'el'
% linear triangles only, so far

nh = length(vert);
ne = length(idx);

flag = zeros(length(re),2);

for i = 1:length(re)
    flag(i,:) = [re(i), 0];
end

new_vert = zeros(10,3);
new_vert_ind = zeros(10,1);

for i = 1:length(re)
    if(flag(i,2) ~= 1)
        mid_vert = zeros(3,2);
        mid_vert_ind = zeros(3,1);
        
        elem = idx(re(i),:);
        elem_ind = re(i);
        
        new_vert(1,:) = vert(elem(1),:); new_vert_ind(1) = elem(1);
        new_vert(2,:) = vert(elem(2),:); new_vert_ind(2) = elem(2);
        new_vert(3,:) = vert(elem(3),:); new_vert_ind(3) = elem(3);
        new_vert(4,:) = vert(elem(4),:); new_vert_ind(4) = elem(4);
        new_vert(5,:) = (new_vert(1,:) + new_vert(2,:))/2;
        vert(nh+1,:) = new_vert(5,:); new_vert_ind(5) = nh+1;
        new_vert(6,:) = (new_vert(1,:) + new_vert(3,:))/2;
        vert(nh+2,:) = new_vert(6,:); new_vert_ind(6) = nh+2;
        new_vert(7,:) = (new_vert(1,:) + new_vert(4,:))/2;
        vert(nh+3,:) = new_vert(7,:); new_vert_ind(7) = nh+3;
        new_vert(8,:) = (new_vert(2,:) + new_vert(3,:))/2;
        vert(nh+4,:) = new_vert(8,:); new_vert_ind(8) = nh+4;
        new_vert(9,:) = (new_vert(3,:) + new_vert(4,:))/2;
        vert(nh+5,:) = new_vert(9,:); new_vert_ind(9) = nh+5;
        new_vert(10,:) = (new_vert(2,:) + new_vert(4,:))/2;
        vert(nh+6,:) = new_vert(10,:); new_vert_ind(10) = nh+6;
        nh = nh+6;
        
        idx(ne+1,:) = [new_vert_ind(1), new_vert_ind(5), new_vert_ind(6), new_vert_ind(7)];
        idx(ne+2,:) = [new_vert_ind(2), new_vert_ind(5), new_vert_ind(8), new_vert_ind(10)];
        idx(ne+3,:) = [new_vert_ind(3), new_vert_ind(6), new_vert_ind(8), new_vert_ind(9)];
        idx(ne+4,:) = [new_vert_ind(4), new_vert_ind(7), new_vert_ind(9), new_vert_ind(10)];
        idx(ne+5,:) = [new_vert_ind(6), new_vert_ind(7), new_vert_ind(9), new_vert_ind(10)];
        idx(ne+6,:) = [new_vert_ind(6), new_vert_ind(8), new_vert_ind(9), new_vert_ind(10)];
        idx(ne+7,:) = [new_vert_ind(5), new_vert_ind(6), new_vert_ind(8), new_vert_ind(10)];
        %idx(ne+8,:) = [new_vert_ind(5), new_vert_ind(6), new_vert_ind(7), new_vert_ind(10)];
        ne = ne+7;
        
        %Add the three new vertices
        for j = 1:4
            
            vert1 = elem(mod(j-1,4)+1);
            vert2 = elem(mod(j,4)+1);
            vert3 = elem(mod(j+1,4)+1);
            
            flag1 = 0;
            
            for k = 1:4
                if(neighbors(elem_ind,k) ~= 0 && flag1 ~= 1)
                    
                    com_vert_temp = intersect(elem, idx(neighbors(elem_ind, k),:));
                    
                    I = find(idx(neighbors(elem_ind,k),:) == vert1);
                    J = find(idx(neighbors(elem_ind,k),:) == vert2);
                    K = find(idx(neighbors(elem_ind,k),:) == vert3);
                    flag1 = 0;
                    
                    if(~isempty(I) && ~isempty(J) && ~isempty(K))
                        com_vert = [idx(neighbors(elem_ind,k),I), idx(neighbors(elem_ind,k),J), idx(neighbors(elem_ind,k),K)];
                        opp_vert = [idx(neighbors(elem_ind,k),10 - I - J - K)];
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
                    else
                        no_add = 1;
                    end
                end
            end
            
            if(no_add ~= 1 && neighbors(elem_ind,index) ~= 0)
                if(j == 1)
                    elem_ind
                    index
                    neighbors(elem_ind,index)
                    [com_vert(1), new_vert_ind(5), new_vert_ind(6), opp_vert]
                    idx(neighbors(elem_ind,index),:) = [com_vert(1), new_vert_ind(5), new_vert_ind(6), opp_vert];
                    
                    idx(ne+1,:) = [com_vert(2); new_vert_ind(5); new_vert_ind(8); opp_vert];
                    idx(ne+2,:) = [com_vert(3); new_vert_ind(6); new_vert_ind(8); opp_vert];
                    idx(ne+3,:) = [new_vert_ind(5); new_vert_ind(6); new_vert_ind(8); opp_vert];
                    ne = ne+3;
                    
                    remove_ind = find(re == neighbors(elem_ind,index));
                    if(remove_ind ~= 0)
                        flag(remove_ind,2) = 1;
                    end
                elseif(j == 2)
                    idx(neighbors(elem_ind,index),:) = [com_vert(1); new_vert_ind(8); new_vert_ind(10); opp_vert];
                    
                    idx(ne+1,:) = [com_vert(2); new_vert_ind(8); new_vert_ind(9); opp_vert];
                    idx(ne+2,:) = [com_vert(3); new_vert_ind(9); new_vert_ind(10); opp_vert];
                    idx(ne+3,:) = [new_vert_ind(8); new_vert_ind(9); new_vert_ind(10); opp_vert];
                    ne = ne+3;
                    
                    remove_ind = find(re == neighbors(elem_ind,index));
                    if(remove_ind ~= 0)
                        flag(remove_ind,2) = 1;
                    end
                elseif(j == 3)
                    idx(neighbors(elem_ind,index),:) = [com_vert(1); new_vert_ind(6); new_vert_ind(9); opp_vert];
                    
                    idx(ne+1,:) = [com_vert(2); new_vert_ind(7); new_vert_ind(9); opp_vert];
                    idx(ne+2,:) = [com_vert(3); new_vert_ind(6); new_vert_ind(7); opp_vert];
                    idx(ne+3,:) = [new_vert_ind(6); new_vert_ind(7); new_vert_ind(9); opp_vert];
                    ne = ne+3;
                    
                    remove_ind = find(re == neighbors(elem_ind,index));
                    if(remove_ind ~= 0)
                        flag(remove_ind,2) = 1;
                    end
                elseif(j == 4)
                    idx(neighbors(elem_ind,index),:) = [com_vert(1); new_vert_ind(7); new_vert_ind(10); opp_vert];
                    
                    idx(ne+1,:) = [com_vert(2); new_vert_ind(5); new_vert_ind(7); opp_vert];
                    idx(ne+2,:) = [com_vert(3); new_vert_ind(5); new_vert_ind(10); opp_vert];
                    idx(ne+3,:) = [new_vert_ind(5); new_vert_ind(7); new_vert_ind(10); opp_vert];
                    ne = ne+3;
                    
                    remove_ind = find(re == neighbors(elem_ind,index));
                    if(remove_ind ~= 0)
                        flag(remove_ind,2) = 1;
                    end
                    
                end
            else
                no_add = 0;
            end
            
        end
        
        
        idx(re(i),:) = [new_vert_ind(5), new_vert_ind(6), new_vert_ind(7), new_vert_ind(10)];
        
        h = toastMesh(vert,idx,3*ones(ne,1));
        
        neighbors = h.ElementNeighbours();
    end
end

sizes = h.ElementSize();
[vertnew, idxnew] = h.Data();

for i = 1:ne
    if(sizes(i) <= 0)
        %sizes(i)
        temp = idxnew(i, :);
        idxnew(i,:) = [temp(2), temp(1), temp(3), temp(4)];
        %h = toastMesh(vertnew,idxnew,15*ones(ne,1));
        %sizes = h.ElementSize();
        %sizes(i)
    end
end

h = toastMesh(vertnew,idxnew,3*ones(ne,1));

