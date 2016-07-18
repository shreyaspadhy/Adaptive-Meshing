function [cond] = condition_of_mesh(hMesh)

[vert, idx] = hMesh.Data();

for i = 1:length(idx)
    a = sqrt((vert(idx(i,3),2) - vert(idx(i,2),2))^2 + (vert(idx(i,3),1) - vert(idx(i,2),1))^2);
    b = sqrt((vert(idx(i,3),2) - vert(idx(i,1),2))^2 + (vert(idx(i,3),1) - vert(idx(i,1),1))^2);
    c = sqrt((vert(idx(i,2),2) - vert(idx(i,1),2))^2 + (vert(idx(i,2),1) - vert(idx(i,1),1))^2);
    
    cosel(1) = (a^2 + b^2 - c^2)/(2.*a.*b);
    cosel(2) = (b^2 + c^2 - a^2)/(2.*b.*c);
    cosel(3) = (c^2 + a^2 - b^2)/(2.*a.*c);
    
    cond(i) = max(cosel);
end