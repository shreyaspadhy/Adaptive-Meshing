function [h, mua] = refine_mesh_2015(vert,el,re, mua0, mua)
% refine the set 're' in list 'el'
% linear triangles only, so far

nh = length(vert);
ne = length(el);
mua = mua0;

for j = 1:length(re)
    ee = el(re(j),:);
    newv = sum( vert(ee,:))/3;
    vert(nh+1,:) = newv;
    nh = nh+1;
    el(re(j),:)   = [ee(1), ee(2),nh];
    el(ne+1,:) = [ee(2), ee(3),nh];
    mua = [mua; mua0(re(j))];
    el(ne+2,:) = [ee(3), ee(1),nh];     
    mua = [mua; mua0(re(j))];
    ne = ne + 2;
end
h = toastMesh(vert,el,15*ones(ne,1));