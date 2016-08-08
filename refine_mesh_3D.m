function [h] = refine_mesh_3D(vert,el,re)
% refine the set 're' in list 'el'
% linear triangles only, so far

nh = length(vert);
ne = length(el);

for j = 1:length(re)
    ee = el(re(j),:);
    newv = sum( vert(ee,:))/4;
    vert(nh+1,:) = newv;
    nh = nh+1;
    el(re(j),:)   = [ee(1), ee(2), ee(3), nh];
    el(ne+1,:) = [ee(2), ee(3), ee(4), nh];
    el(ne+2,:) = [ee(3), ee(4), ee(1), nh];  
    el(ne+3,:) = [ee(1), ee(2), ee(4), nh];
    ne = ne + 3;
end
h = toastMesh(vert,el,3*ones(ne,1));

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