% 
% simple attempt to do adaptive meshing 
% this version is "Toast 2"
%

meshpath = '/home/spadhy/Desktop/toastpp/test/2D/meshes/';

%meshname = 'hires_tri3.msh';
%meshname = 'lores_tri3.msh';

%meshname = 'tinycircle.msh';

%me shname  =  ' circ le25 _ 32. msh';
%meshname = 'circle25_64.msh';
%meshname = 'heads99/neonat_res10_fixed.opt';
qmname = 'circle25_1x1.qm';    

%bx =  64; by =   64;

%% Create Mesh and set attributes
% rad = 25;   % mesh radius [mm]
% nsect = 6;  % number of sectors
% nring = 10; % number of rings
% nbnd = 2;   % number of boundary rings
% 
% [vtx,idx,eltp] = mkcircle(rad,nsect,nring,nbnd);
%             % create the mesh geometry
% hMesh = toastMesh (vtx,idx,eltp);
%             % create the mesh object
% 
%            
% % get mesh
% %hMesh = toastMesh([meshpath meshname]);
% hMesh.ReadQM([meshpath qmname]);
% nh = hMesh.NodeCount();
% ne = hMesh.ElementCount();
% [vert,idx] = hMesh.Data();
% 
% BB= hMesh.BoundingBox();
% bmin = BB(1,:);
% bmax = BB(2,:);
% 
% figure(4); clf;
% subplot(2,2,1); hMesh.Display();
% 
% hBasis = toastBasis(hMesh,[bx,by]);

%% create sysmat
% mua0 = 0.01;
% mus0 = 1;
% kap0 = 1/(3*(mua0+mus0));
% smat = dotSysmat(hMesh,mua0*ones(nh,1),mus0*ones(nh,1),ones(nh,1),0);

%% read in source from file
% mua0 = toastNim([meshpath 'tgt_mua_ellips_tri10.nim']);
% mus0 = toastNim([meshpath 'tgt_mus_ellips_tri10.nim']);
% smat = dotSysmat(hMesh, mua0, mus0, ones(ne,1),0, 'EL');

%% read in parameters from image
% bmua = imread('demo_matlab_fwd2_mua.png');
% bmus = imread('demo_matlab_fwd2_mus.png');
% 
% bmua = double(bmua)./255.*0.02 + 0.01;
% bmus = double(bmus)./255.*1.0 + 1.0;
% 
% grd = size(bmua);
% basis = toastBasis(hMesh,grd);
% 
% mua0 = basis.Map('B->M', bmua);
% mus0 = basis.Map('B->M', bmus);
% 
% kap0 = 1./(3.*(mua0 + mus0));
% 
% smat = dotSysmat(hMesh, mua0, mus0, ones(nh,1),0);
% 
% mua0_coarse = mua0;

%% Parameters on element basis

hMesh = toastMesh('circle_blob_2.msh','gmsh');
ne = hMesh.ElementCount;
nh = hMesh.NodeCount;
regidx = hMesh.Region;
regno = unique(regidx);
blobel = find(regidx == regno(2)); % assuming that the second surface marks the inclusion

refind = 1.4; % refractive index
c0 = 0.3; % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.01; % background absorption [1/mm]
mus_bkg = 1; % background scattering [1/mm];
ref = ones(ne,1)*refind;
mua0 = ones(ne,1)*mua_bkg;
mus0 = ones(ne,1)*mus_bkg;
mua0(blobel) = mua_bkg*2;
kap0 = 1./(3.*(mua0+mus0));
hMesh.ReadQM([meshpath qmname]);
[vert,idx] = hMesh.Data();
grd = [128,128];

smat = dotSysmat(hMesh, mua0, mus0, ones(ne,1),0, 'EL');

%% Reconstruction parameters

mua = 0.01*ones(ne,1);
mus = 1*ones(ne,1);
kap = 1./(3.*(mua+mus));

%% Sources

Q = hMesh.Qvec ( 'Neumann', 'Gaussian', 2);
nQ = size(Q,2);


%% Solve
phi = smat \Q;

%% Residue Related Variables
res = zeros(ne,1);
grad = zeros(ne,2);
grad_el = zeros(ne,2);
coords = zeros(ne,2);

%% Loop over sources

for qq = 1:1:nQ
figure(1);clf;
subplot(2,2,1); hMesh.Display(full(phi(:,qq)));title('Field');
minp = min(phi(:,qq));
subplot(2,2,2); hMesh.Display(log(abs(full(phi(:,qq)))+abs(1.01*minp)));title('Field on log scale');

%% global residual
%res = full(smat*phi(:,qq) - Q(:,qq));

%% compute basic a-posteriori estimate
% for i = 1:ne
%     coord = (hMesh.Element(i).Data());
%     x_coord = mean(coord(:,1));
%     y_coord = mean(coord(:,2));
%     coords(i,:) = [x_coord, y_coord];
% end

for i = 1:ne
    coord = (hMesh.Element(i).Data());
    %x_coord = mean(coord(:,1));
    %y_coord = mean(coord(:,2));
    coords(i,:) = mean(coord);
    temp = hMesh.Element(i).ShapeDer(coords(i,:)', 'global');
    grad_el(i,:) = [phi(idx(i,1,qq)), phi(idx(i,2,qq)), phi(idx(i,3,qq))]*temp';
    
end

grad_el = grad_el./3;

grad_el(:,1) = kap0(:).*grad_el(:,1);
grad_el(:,2) = kap0(:).*grad_el(:,2);

neighbors = hMesh.ElementNeighbours();
div = zeros(ne,1);

for i = 1:ne
    for j = 1:length(neighbors(i,:))
        if(neighbors(i,j) ~= 0)
            ds = coords(neighbors(i,j),:) - coords(i,:);
            diff = grad_el(neighbors(i,j),:) - grad_el(i,:);
            dot(diff,ds)
            norm(ds)
            div(i) = div(i) + dot(diff,ds)./norm(ds);
        end
    end
end

second_term = zeros(ne, 1);

for i = 1:ne
    temp2 = hMesh.Element(i).ShapeFun(coords(i,:)', 'global');
    second_term(i) = mua(i).*([phi(idx(i,1,qq)), phi(idx(i,2,qq)), phi(idx(i,3,qq))]*temp2);
end

res = abs(-div + second_term);
%res = -(grad(:,1)+grad(:,2)) + mua0.*phi;
%res = -div(:,2) + mua.*phi
    
%res_el_plot = basis.Map('M->B', res_el);
%imagesc(res_el_plot);
minr = min(res);
subplot(2,2,3); hMesh.Display(full(res));title('Error Estimate');
subplot(2,2,4); hMesh.Display(log(abs(res)+abs(1.01*minr)));title('Error Estimate on log scale');
pause(0.5);
end

count = 0;
maxr = max(res_el);

for i = 1:length(res)
    if (res(i) > 0.2*1e-10)
        count = count+1;
    end
end

j = 1;
re = zeros(count,1);
for i = 1:length(res)
    if(res(i) > 0.2*1e-10)
        re(j) = i;
        j = j + 1;
    end
end

    
    

%% globally refine
h2 = refine_mesh_2015(vert,idx,re);
nh2 = h2.NodeCount();%(log(abs(full(phi2(:,qq)))+abs(1.01*minp)))
ne2 = h2.ElementCount();

h2.ReadQM([meshpath qmname]); % not great programming. Should not need to read twice...

Q2 = h2.Qvec ('Neumann', 'Gaussian', 2);

figure(5);
h2.Display();

basis2 = toastBasis(h2,grd);

mua0 = basis2.Map('B->M', bmua);
mus0 = basis2.Map('B->M', bmus);

smat2 = dotSysmat(h2,mua0,mus0,ones(nh2,1),0);
phi2 = smat2 \Q2;

for qq = 1:4:nQ
figure(2);clf;
subplot(2,3,1); h2.Display(full(phi2(:,qq)));title('Field');
minp = min(phi2(:,qq));
subplot(2,3,2); h2.Display(log(abs(full(phi2(:,qq)))+abs(1.01*minp)));title('Field on log scale');
subplot(2,3,3); hMesh.Display(mua0_coarse);title('\mu_a Original');


res2 = full(smat2*phi2(:,qq) - Q2(:,qq));
minr2 = min(res2);
subplot(2,3,4); h2.Display(full(res2));title('Residual');

subplot(2,3,5); h2.Display(log(abs(res2)+abs(1.01*minr2)));title('Residual on log scale');
subplot(2,3,6); h2.Display(mus0);title('\mu_s Original');

pause(0.5);

figure(3);clf;
subplot(2,1,1); hMesh.Display();
subplot(2,1,2); h2.Display();
end

%%
%toastDeleteBasis(hBasis);
%toastDeleteMesh(hMesh);