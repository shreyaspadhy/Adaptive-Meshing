% Adaptive Meshing with local residue calculation

meshpath = '/home/spadhy/Desktop/toastpp/test/2D/meshes/';
qmname = 'circle25_1x1.qm';

%% Read Mesh and Parameters on element basis

hMesh = toastMesh('circle_blob_2.msh','gmsh');
ne = hMesh.ElementCount;
nh = hMesh.NodeCount;
%regidx = hMesh.Region;
%regno = unique(regidx);
%blobel = find(regidx == regno(2)); % assuming% basis2 = toastBasis(h2,grd);

refind = 1.4; % refractive index
c0 = 0.3; % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.01; % background absorption [1/mm]
mus_bkg = 1; % background scattering [1/mm];


bmua = imread('demo_matlab_fwd2_mua.png');
bmus = imread('demo_matlab_fwd2_mus.png');
bmua = double(bmua)./255.*0.02 + 0.01;
bmus = double(bmus)./255.*1.0 + 1.0;
grd = size(bmua);

basis = toastBasis(hMesh,grd);
mua_node = basis.Map('B->M',bmua);
mus_node = basis.Map('B->M',bmus);

hMesh.ReadQM([meshpath qmname]);
[vert,idx] = hMesh.Data();
figure(1); clf;
hMesh.Display();
% Mesh has been made and properties populated

ref = ones(ne,1)*refind;
mua0 = zeros(ne,1);
mus0 = zeros(ne,1);

%% Transfer from node basis to element basis
for i = 1:ne
    coord = (hMesh.Element(i).Data());
    %x_coord = mean(coord(:,1));
    %y_coord = mean(coord(:,2));
    coords(i,:) = mean(coord);
    temp = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
    mua0(i) = ([mua_node(idx(i,1)), mua_node(idx(i,2)), mua_node(idx(i,3))]*temp);
    mus0(i) = ([mus_node(idx(i,1)), mus_node(idx(i,2)), mus_node(idx(i,3))]*temp);
end

%mua0 = ones(ne,1)*mua_bkg;
%mus0 = ones(ne,1)*mus_bkg;
figure(10);
subplot(2,1,1); hMesh.Display(mua0);
subplot(2,1,2); hMesh.Display(mus0);

kap0 = 1./(3.*(mua0+mus0));

%% Sources
Q = hMesh.Qvec ( 'Neumann', 'Gaussian', 2);
nQ = size(Q,2);

%% Solve
smat = dotSysmat(hMesh, mua0, mus0, ones(ne,1),0, 'EL');
phi = smat \Q;

%% Loop over number of refinements
for ad = 1:4
    %% Residue Related Variables
    res = zeros(ne,1);
    grad = zeros(ne,2);
    grad_el = zeros(ne,2);
    coords = zeros(ne,2);
    
    %% Loop over sources
    for qq = 1:1:nQ
        figure(2);clf;
        subplot(2,2,1); hMesh.Display(full(phi(:,qq)));title('Field');
        minp = min(phi(:,qq));
        subplot(2,2,2); hMesh.Display(log(abs(full(phi(:,qq)))+abs(1.01*minp)));title('Field on log scale');
        
        neighbors = hMesh.ElementNeighbours();
        %% Compute local residue estimate
        [res] = calc_residue_jump(hMesh, phi, idx, kap0, neighbors, mua0, ne, qq, Q, vert);
        
        minr = min(res);
        subplot(2,2,3); hMesh.Display(full(res));title('Error Estimate');
        subplot(2,2,4); hMesh.Display(log(abs(res)+abs(1.01*minr)));title('Error Estimate on log scale');
        pause(0.5);
    end
    
    %% Find elements that need refinement
    dns = log(res);
    dns = (dns-min(dns))/(max(dns)-min(dns))*1.8+0.2;
    
    %Choose the top 10% of residues to refine
    res_sort = zeros(ne,2);
    for i = 1:ne
        res_sort(i,1) = res(i);
        res_sort(i,2) = i;
    end
    
    res_sort2 = sortrows(res_sort,1);
    small = 0;
    num = round(ne/10);
    
    re = zeros(num,1);
    for p = 1:1:num
        re(p) = res_sort2(ne-p,2);
        res_sort2(ne-p,2)
    end
    
    %% Average res values of each surrounding element to a vertex
    res_nh = zeros(nh,1);
    for i = 1:nh
        I = find(idx == i)
        for j = 1:length(I)
            a = mod(I(j),ne);
            if(a == 0)
                a = 1400;
            end
            res_nh(i) = res_nh(i) + res(a);
        end
        res_nh(i) = res_nh(i)./length(I);
    end
    
    figure(22);
    minr_nh = min(res_nh);
    subplot(2,1,1); hMesh.Display(log(abs(res_nh)+abs(1.01*minr_nh)));
    subplot(2,1,2); hMesh.Display(res_nh);
    
    %mua = zeros(ne,1);
    
    %% Refine using gmsh
    
    %dns = -log(res_nh);
    %dns = (dns-min(dns))/(max(dns)-min(dns))*1.8+0.2;
    
           %hMesh.Write('tmp.pos','gmsh-pos',dns);
    % %
    % %
           %system('gmsh -bgm tmp.pos -merge circle.geo -2 -o adaptive.msh');
    % %
          %h2 = toastMesh('adaptive.msh','gmsh');
    % %     h2.Display;
    %% Refine using Meshing Strategies (i), (ii), (iii)
        %[h2, mua] = refine_mesh_2015(vert,idx,re, mua0, mua);
    [h2] = refine_mesh_sierpinski(vert, idx, re, mua0, neighbors,qq);
    
    
    %% Populate new mesh properties
    
    nh2 = h2.NodeCount();%(log(abs(full(phi2(:,qq)))+abs(1.01*minp)))
    ne2 = h2.ElementCount();
    
    h2.ReadQM([meshpath qmname]); % not great programming. Should not need to read twice...
    
    Q2 = h2.Qvec ('Neumann', 'Gaussian', 2);
    
    [vert2, idx2] = h2.Data();
    figure(14);clf;
    h2.Display();
    
    sizes = h2.ElementSize();
    
    %% Work on refined mesh
    
    basis2 = toastBasis(h2,grd);
    mua_node2 = basis2.Map('B->M', bmua);
    mus_node2 = basis2.Map('B->M', bmus);
    
    mua = zeros(ne2,1);
    mus = zeros(ne2,1);
    
    for i = 1:ne2
        coord = (h2.Element(i).Data());
        %x_coord = mean(coord(:,1));
        %y_coord = mean(coord(:,2));
        coords2(i,:) = mean(coord);
        temp2 = h2.Element(i).ShapeFun(coords2(i,:)', 'global');
        mua(i) = ([mua_node2(idx2(i,1,qq)), mua_node2(idx2(i,2,qq)), mua_node2(idx2(i,3,qq))]*temp2);
        mus(i) = ([mus_node2(idx2(i,1,qq)), mus_node2(idx2(i,2,qq)), mus_node2(idx2(i,3,qq))]*temp2);
    end
    %mua = basis2.Map('B->M', bmua);
    %mua = mua_bkg*ones(ne2,1);
    %mus = mus_bkg*ones(ne2,1);
    
    kap = 1./(3.*(mua+mus));
    
    figure(11);
    subplot(2,1,1); h2.Display(mua);
    subplot(2,1,2); h2.Display(mus);
    
    smat2 = dotSysmat(h2,mua,mus,ones(ne2,1),0, 'EL');
    
    phi2 = smat2 \Q2;
    
    for qq = 1:4:nQ
        figure(5);clf;
        subplot(2,2,1); h2.Display(full(phi2(:,qq)));title('Field');
        minp = min(phi2(:,qq));
        subplot(2,2,2); h2.Display(log(abs(full(phi2(:,qq)))+abs(1.01*minp)));title('Field on log scale');
        %subplot(2,2,3); hMesh.Display(mua0);title('\mu_a Original');
        
        
        %% Calculate residue on refined mesh
        res2 = zeros(ne2, 1);
        %coords2 = zeros(ne2,2);
        grad_el2 = zeros(ne2,2);
        neighbors2 = h2.ElementNeighbours();
        
        [res2] = calc_residue_jump(h2, phi2, idx2, kap, neighbors2, mua, ne2, qq, Q2, vert2);
        
        %% Display stuff
        %res2 = full(smat2*phi2(:,qq) - Q2(:,qq));
        minr2 = min(res2);
        subplot(2,2,3); h2.Display(full(res2));title('Residual');
        
        subplot(2,2,4); h2.Display(log(abs(res2)+abs(1.01*minr2)));title('Residual on log scale');
        %subplot(2,3,6); h2.Display(mus);title('\mu_s Original');
        
        figure(15);
        subplot(2,1,1); hMesh.Display(full(res), 'range', [0, max(res)]);
        subplot(2,1,2); h2.Display(full(res2), 'range', [0, max(res)]);
        pause(0.5);
        
        res_nh2 = zeros(nh2,1);
        for i = 1:nh2
            I = find(idx2 == i)
            for j = 1:length(I)
                a = mod(I(j),ne2);
                if(a == 0)
                    a = ne2;
                end
                res_nh2(i) = res_nh2(i) + res2(a);
            end
            res_nh2(i) = res_nh2(i)./length(I);
        end
        
        figure(22);
        minr_nh2 = min(res_nh2);
        subplot(2,1,1); h2.Display(log(abs(res_nh2)+abs(1.01*minr_nh2)));
        subplot(2,1,2); h2.Display(res_nh2);
        %         figure(5);clf;
        %         subplot(2,1,1); hMesh.Display();
        %         subplot(2,1,2); h2.Display(mua);
    end
    
    hMesh = h2;
    ne = ne2;
    nh = nh2;
    vert = vert2;
    idx = idx2;
    mua0 = mua;
    mus0 = mus;
    kap0 = kap;
    phi = phi2;
    neighbors = neighbors2;
    coords = coords2;
    Q = Q2;
end

% count = 0;
% maxr = max(res2);
%
% for i = 1:length(res2)
%     if (res2(i) > 0.03*1e-10)
%         count = count+1;
%     end
% end
%
% j = 1;
% re2 = zeros(count,1);
% for i = 1:length(res2)
%     if(res2(i) > 0.03*1e-10)
%         re2(j) = i;
%         j = j + 1;
%     end
% end
% mua3 = zeros(ne2,1);
% [h3, mua3] = refine_mesh_2015(vert2,idx2,re2, mua, mua3);
% nh3 = h3.NodeCount();%(log(abs(full(phi2(:,qq)))+abs(1.01*minp)))
% ne3 = h3.ElementCount();
%
% h3.ReadQM([meshpath qmname]); % not great programming. Should not need to read twice...
%
% Q3 = h3.Qvec ('Neumann', 'Gaussian', 2);
%
% [vert3, idx3] = h3.Data();
% figure(11);
% h3.Display();

%%
%toastDeleteBasis(hBasis);
%toastDeleteMesh(hMesh);