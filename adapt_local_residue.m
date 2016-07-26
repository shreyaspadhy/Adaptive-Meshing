% Adaptive Meshing with local residue calculation
function adapt_local_residue

%meshpath = '/home/spadhy/Desktop/toastpp/test/2D/meshes/';
meshpath = 'C:\Users\shrey\Desktop\TOAST\toast\toast\Adaptive-Meshing\';
%qmname = 'circle25_32x32.qm';

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
grd_img = size(bmua);

basis_img = toastBasis(hMesh,grd_img);
mua_node = basis_img.Map('B->M',bmua);
mus_node = basis_img.Map('B->M',bmus);

%grid.ElRef()
%hMesh.ReadQM([meshpath qmname]);
[vert,idx] = hMesh.Data();
% figure(1); clf;
% hMesh.Display();
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
subplot(2,1,1); hMesh.Display(mua0); title('\mu_a on coarse mesh');
subplot(2,1,2); hMesh.Display(mus0); title('\mu_s on coarse mesh');

kap0 = 1./(3.*(mua0+mus0));

%% Sources
% Q = hMesh.Qvec ( 'Neumann', 'Gaussian', 2);
%nQ = size(Q,2);
% mvec = hMesh.Mvec('Gaussian',2,refind);

rad = 25; % mesh radius [mm]
nopt = 10;
for i=1:nopt
  phiq = (i-1)/nopt*2*pi;
  qpos(i,:) = rad*[cos(phiq), sin(phiq)];
  phim = (i-0.5)/nopt*2*pi;
  mpos(i,:) = rad*[cos(phim), sin(phim)];
end
hMesh.SetQM(qpos,mpos);
Q = real(hMesh.Qvec('Neumann','Gaussian',2));
mvec = real(hMesh.Mvec('Gaussian',2,refind));
nQ = nopt;

%% Solve
smat = dotSysmat(hMesh, mua0, mus0, ones(ne,1),0, 'EL');
phi = smat \Q;

res_source = zeros(ne,nQ);
%% Loop over number of refinements
for ad = 1:2
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
        [res_source(:,qq)] = calc_residue_jump(hMesh, phi, idx, kap0, neighbors, mua0, ne, qq, Q, vert);
        res = res_source(:,qq);
        minr = min(res);
        subplot(2,2,3); hMesh.Display(full(res));title('Error Estimate');
        subplot(2,2,4); hMesh.Display(log(abs(res)+abs(1.01*minr)));title('Error Estimate on log scale');
        pause(0.5);
    end
    
    res = sum(res_source,2);
    %% Find elements that need refinement
    dns = log(res);
    dns = (dns-min(dns))/(max(dns)-min(dns))*1.8+0.2;
    
    %Choose the top 20% of residues to refine
    res_sort = zeros(ne,2);
    for i = 1:ne
        res_sort(i,1) = res(i);
        res_sort(i,2) = i;
    end
    
    res_sort2 = sortrows(res_sort,1);
    small = 0;
    num = 4*round(ne/10);
    
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
    subplot(2,1,2); hMesh.Display(log(abs(res_nh)+abs(1.01*minr_nh))); title('Log of residue at each element center');
    subplot(2,1,1); hMesh.Display(res_nh); title('Residue at each element center');
    
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
    
    %h2 = toastMesh('circle_blob.msh', 'gmsh');
    
    %% Populate new mesh properties
    
    nh2 = h2.NodeCount();%(log(abs(full(phi2(:,qq)))+abs(1.01*minp)))
    ne2 = h2.ElementCount();
    
    %h2.ReadQM([meshpath qmname]); % not great programming. Should not need to read twice...
    
    %     Q2 = h2.Qvec ('Neumann', 'Gaussian', 2);
    %     mvec2 = h2.Mvec('Gaussian',2, refind);
    
  
    h2.SetQM(qpos,mpos);
    Q2 = real(h2.Qvec('Neumann','Gaussian',2));
    mvec2 = real(h2.Mvec('Gaussian',2,refind));
    
    [vert2, idx2] = h2.Data();
    figure(14);clf;
    h2.Display(); title('Refined mesh');
    
    sizes = h2.ElementSize();
    
    %% Work on refined mesh
    
    basis2 = toastBasis(h2,grd_img);
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
        mua(i) = ([mua_node2(idx2(i,1)), mua_node2(idx2(i,2)), mua_node2(idx2(i,3))]*temp2);
        mus(i) = ([mus_node2(idx2(i,1)), mus_node2(idx2(i,2)), mus_node2(idx2(i,3))]*temp2);
    end
    %mua = basis2.Map('B->M', bmua);
    %mua = mua_bkg*ones(ne2,1);
    %mus = mus_bkg*ones(ne2,1);
    
    kap = 1./(3.*(mua+mus));
    
    figure(11);
    subplot(2,1,1); h2.Display(mua);
    subplot(2,1,2); h2.Display(mus);
    
    res2_source = zeros(ne2, nQ);
    res2d = zeros(ne2,1);
    
    smat2 = dotSysmat(h2,mua_node2,mus_node2,refind*ones(nh2,1),0);
    
    phi2 = smat2 \Q2;
    data_model_f = reshape(log(mvec2'*phi2), [], 1); %This is the actual data generated by forward solver on refined mesh
    for qq = 1:1:nQ
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
        
        [res2_source(:,qq)] = calc_residue_jump(h2, phi2, idx2, kap, neighbors2, mua, ne2, qq, Q2, vert2);
        res2 = res2_source(:,qq);
        %% Display stuff
        %res2 = full(smat2*phi2(:,qq) - Q2(:,qq));
        minr2 = min(res2);
        subplot(2,2,3); h2.Display(full(res2));title('Residual');
        
        subplot(2,2,4); h2.Display(log(abs(res2)+abs(1.01*minr2)));title('Residual on log scale');
        %subplot(2,3,6); h2.Display(mus);title('\mu_s Original');
    end
    res2d = sum(res2_source,2);
        figure(15); title('Comparing rediues between coarse and fine mesh');
        subplot(2,1,1); hMesh.Display(full(res), 'range', [0, max(res)]);
        subplot(2,1,2); h2.Display(full(res2d), 'range', [0, max(res)]);
        pause(0.5);
        
        res_nh2 = zeros(nh2,1);
        for i = 1:nh2
            I = find(idx2 == i)
            for j = 1:length(I)
                a = mod(I(j),ne2);
                if(a == 0)
                    a = ne2;
                end
                res_nh2(i) = res_nh2(i) + res2d(a);
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
    mvec = mvec2;
    
    %% Solve the Inverse Problem on refined mesh
    
%     mesh_f = toastMesh('adaptive.msh','gmsh');
%     basisf = toastBasis(mesh_f,grd_img);
%     mua_node_f = basisf.Map('B->M',bmua);
%     mus_node_f = basisf.Map('B->M',bmus);
%     ref_f = refind*ones(mesh_f.NodeCount,1);
%     mesh_f.SetQM(qpos,mpos);
%     Qf = real(mesh_f.Qvec('Neumann','Gaussian',2));
% 
%     mvecf = real(mesh_f.Mvec('Gaussian',2,refind));    
%     sysMatf = dotSysmat(mesh_f, mua_node_f, mus_node_f, ref_f);
%     data_model_f = reshape(log(mvecf'*(sysMatf\Qf)), [], 1);
    
    
    mua_r = ones(nh,1)*mua_bkg;
    %mus_r = ones(nh,1)*mus_bkg;
    mus_r = mus_node2;
    ref = ones(nh,1)*refind;
    cm = 0.3/refind;
    
    smat = dotSysmat(hMesh, mua_r, mus_r, ref);
    proj = reshape(log(mvec'*(smat\Q)), [], 1);
    sd = ones(size(proj));
    %data_model = reshape(data_model_f, [], 1);
    data_model = data_model_f;
    grd_inv = [256, 256];
    basis3 = toastBasis(hMesh, grd_inv);
    
    bmua_r = basis3.Map('M->B', mua_r);
    bcmua = bmua_r*cm;
    scmua = basis3.Map('B->S', bcmua);
    x = scmua;
    logx = log(x);
    slen = length(x);
    
    regul = toastRegul('TV', basis3, logx, 1e-4, 'Beta', 0.01);
    
    err0 = toastObjective(proj, data_model, sd, regul, logx);
    err = err0;
    errp = inf;
    itr = 1;
    step = 0.1;
    fprintf('Iteration %d, objective %f\n', 0, err);
    
    %% Inverse Solver Loop
    
    itrmax = 100; %CG iteration limit
    tolCG = 1e-6;
    
    while (itr <= itrmax) && (err > tolCG*err0) && (errp-err > tolCG)
        
        errp = err;
        
        r = -toastGradient(hMesh, basis3, Q, mvec, mua_r, mus_r, ref, 0, ...
            data_model, sd, 'method', 'cg', 'tolerance', 1e-12);
        r = r(1:slen);
        r = r.*x;
        r = r - regul.Gradient(logx);
        
        %CG search direction update
        if itr > 1
            delta_old = delta_new;
            delta_mid = r'*s;
        end
        s = r;
        if itr == 1
            d = s;
            delta_new = r' * d;
        else
            delta_new = r' * s;
            beta = (delta_new - delta_mid) / delta_old;
            if mod(itr, 20) == 0 || beta <= 0
                d = s;
            else
                d = s + d*beta;
            end
        end
%d = r;
        
        % Line search along update direction
        step = toastLineSearch(logx, d, step, err, @objective);
        logx = logx + d*step; %Update solution estimate
        
        mua_r = basis3.Map('S->M', exp(logx)/cm);
        
        proj = reshape(log(mvec'*(dotSysmat(hMesh, mua_r, mus_r, ref)\Q)), [], 1);
        
        err = toastObjective(proj, data_model, sd, regul, logx);
        
        fprintf('Iteration %d, objective %f\n', itr, err);
        itr = itr+1;
        
        figure(50);
    end
             hMesh.Display(mua_r);
   
end

function p = objective(logx)  
mua_ = basis3.Map('S->M',(exp(logx)/cm));
proj_ = reshape(log(mvec' * (dotSysmat(hMesh, mua_, mus_r, ref)\Q)), [], 1);
p = toastObjective(proj_, data_model, sd, regul, logx);    
end

end
 
