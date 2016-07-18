% Adaptive Meshing with local residue calculation

meshpath = '/home/spadhy/Desktop/toastpp/test/2D/meshes/';
qmname = 'circle25_1x1.qm';

%% Read Mesh and Parameters on element basis

hMesh = toastMesh('circle_standard.msh','gmsh');
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


bmua = imread('demo_matlab_fwd2_mua.png');-
bmus = imread('demo_matlab_fwd2_mus.png');
bmua = double(bmua)./255.*0.02 + 0.01;
bmus = double(bmus)./255.*1.0 + 1.0;
grd = size(bmua);

basis = toastBasis(hMesh,grd);
mua0 = basis.Map('B->M',bmua);
mus0 = basis.Map('B->M',bmus);

ref = ones(nh,1)*refind;
%mua0 = ones(nh,1)*mua_bkg;
%mus0 = ones(nh,1)*mus_bkg;
%mua0(blobel) = mua_bkg*2;
kap0 = 1./(3.*(mua0+mus0));

hMesh.ReadQM([meshpath qmname]);
[vert,idx] = hMesh.Data();

figure(1); clf;
hMesh.Display(mua0);

smat = dotSysmat(hMesh, mua0, mus0, ones(nh,1),0);

%% Reconstruction parameters
% 
% mua = 0.01*ones(ne,1);
% mus = 1*ones(ne,1);
% kap = 1./(3.*(mua+mus));

%% Sources
Q = hMesh.Qvec ( 'Neumann', 'Gaussian', 2);
nQ = size(Q,2);

%% Solve
phi = smat \Q;

%% Loop over number of refinements
for ad = 1:4
    %% Residue Related Variables
    res = zeros(ne,1);
    grad = zeros(ne,2);
    grad_el = zeros(ne,2);
    coords = zeros(ne,2);
    neighbors = hMesh.ElementNeighbours();
    
    %% Loop over sources
    for qq = 1:1:nQ
        figure(2);clf;
        subplot(2,2,1); hMesh.Display(full(phi(:,qq)));title('Field');
        minp = min(phi(:,qq));
        subplot(2,2,2); hMesh.Display(log(abs(full(phi(:,qq)))+abs(1.01*minp)));title('Field on log scale');
        
        %% Compute local residue estimate
        [res, coords] = calc_residue_linear_inverse(hMesh, phi, idx, kap0, neighbors, mua0, ne, qq);

        minr = min(res);
        subplot(2,2,3); hMesh.Display(full(res));title('Error Estimate');
        subplot(2,2,4); hMesh.Display(log(abs(res)+abs(1.01*minr)));title('Error Estimate on log scale');
        pause(0.5);
    end
    
    %% Find elements that need refinement
    %dns = log(res);
    %dns = (dns-min(dns))/(max(dns)-min(dns))*1.8+0.2;
    
    % Choose the top 10% of residues to refine
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
    
    
            
%     index
%     
%     count = 0;
%     maxr = max(res);
%     
%     for i = 1:length(res)
%         if (res(i) > res(index))
%             count = count+1;
%         end
%     end
%     
%     j = 1;
%     re = zeros(count,1);
%     for i = 1:length(res)
%         if(res(i) > res(index))
%             re(j) = i;
%             j = j + 1;
%         end
%     end

    %% Refine using gmsh
     dns = log(res);
     dns = (dns-min(dns))/(max(dns)-min(dns))*3.8+0.2;
     
%     dns_sort = zeros(ne,2);
%     for i = 1:ne
%         dns_sort(i,1) = dns(i);
%         dns_sort(i,2) = i;
%     end
%     
%     dns_sort2 = sortrows(dns_sort,1);
%     small = 0;
%     num = round(ne/10);
%     
%     de = zeros(num,1);
%     for p = 1:1:num
%         dns(dns_sort2(ne-p,2)) = dns(dns_sort2(ne-p,2))*2;
%     end
     
     
     hMesh.Write('tmp.pos','gmsh-pos',dns);
     
     figure(3);clf;
     %hMesh.Display(dns);
     
     system('gmsh -merge -bgm tmp.pos circle.geo -2 -o circle_blob_3.msh');
     %system('gmsh circle_blob_2.msh -bgm tmp.pos -algo meshadapt -2 -o circle_blob_2.msh');

     h2 = toastMesh('circle_blob_3.msh','gmsh');
%     h2.Display;
    %% locally refine
    %[h2, mua] = refine_mesh_2015(vert,idx,re, mua0, mua);
    %[h2, mua] = refine_mesh_sierpinski(vert, idx, re, mua0, mua,  neighbors);
    nh2 = h2.NodeCount();%(log(abs(full(phi2(:,qq)))+abs(1.01*minp)))
    ne2 = h2.ElementCount();
    
    h2.ReadQM([meshpath qmname]); % not great programming. Should not need to read twice...
    
    Q2 = h2.Qvec ('Neumann', 'Gaussian', 2);
    
    [vert2, idx2] = h2.Data();
    figure(4);clf;
    h2.Display();
    
    %% Work on refined mesh
    
    basis2 = toastBasis(h2,grd);
    
    mua = basis2.Map('B->M', bmua);
    mus = basis2.Map('B->M', bmus);
  
    
    kap = 1./(3.*(mua+mus));
    
    smat2 = dotSysmat(h2,mua,mus,ones(nh2,1),0);
    
    phi2 = smat2 \Q2;
    
    for qq = 1:4:nQ
        figure(5);clf;
        subplot(2,2,1); h2.Display(full(phi2(:,qq)));title('Field');
        minp = min(phi2(:,qq));
        subplot(2,2,2); h2.Display(log(abs(full(phi2(:,qq)))+abs(1.01*minp)));title('Field on log scale');
        %subplot(2,2,3); hMesh.Display(mua0);title('\mu_a Original');
        
        
        %% Calculate residue on refined mesh
        res2 = zeros(ne2, 1);
        coords2 = zeros(ne2,2);
        grad_el2 = zeros(ne2,2);
        neighbors2 = h2.ElementNeighbours();
        
        [res2, coords2] = calc_residue_linear_inverse(h2, phi2, idx2, kap, neighbors2, mua, ne2, qq);
        
        %% Display stuff
        %res2 = full(smat2*phi2(:,qq) - Q2(:,qq));
        minr2 = min(res2);
        subplot(2,2,3); h2.Display(full(res2));title('Residual');
        
        subplot(2,2,4); h2.Display(log(abs(res2)+abs(1.01*minr2)));title('Residual on log scale');
        %subplot(2,3,6); h2.Display(mus);title('\mu_s Original');
        
        pause(0.5);
        
%         figure(5);clf;
%         subplot(2,1,1); hMesh.Display();
%         subplot(2,1,2); h2.Display(mua);
    end
    
    hMesh = h2;
    ne = ne2;
    vert = vert2;
    idx = idx2;
    mua0 = mua;
    mus0 = mus;
    kap0 = kap;
    phi = phi2;
    neighbors = neighbors2;
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