function inverse_solver_3D

%% Load Mesh Data
close all;
refind = 1.4; % refractive index
c0 = 0.3; % speed of light in vacuum [mm/ps]
cm = c0/refind; % speed of light in the medium [mm/ps]
mua_bkg = 0.01; % background absorption [1/mm]
mus_bkg = 1; % background scattering [1/mm];

% Load parameters
hmua = toastMesh('C:\Users\shrey\Desktop\TOAST\toast\toast\Adaptive-Meshing\cyl4_blobs.msh');
%Load initial inverse mesh
h0 = toastMesh('C:\Users\shrey\Desktop\TOAST\toast\toast\Adaptive-Meshing\cyl2.msh');
ne = h0.ElementCount;
nh = h0.NodeCount;
[vert,idx] = h0.Data();
ref = ones(ne,1)*refind;

h0.ReadQM('cyl_5ring.qm');
Q = real(h0.Qvec('Neumann','Gaussian',2));
mvec = real(h0.Mvec('Gaussian',2,refind));
nQ = size(Q,2);

grd_img = [32 32 32];

% Store parameters in nodal basis on h0
nim = toastNim('C:\Users\shrey\Desktop\TOAST\toast\toast\Adaptive-Meshing\mua_tgt_cyl4.nim');
mua0 = nim.Values;
nim = toastNim('C:\Users\shrey\Desktop\TOAST\toast\toast\Adaptive-Meshing\mus_tgt_cyl4.nim');
mus0 = nim.Values;

basis_inv = toastBasis(hmua, grd_img);
bmua = basis_inv.Map('M->B', mua0);
bmus = basis_inv.Map('M->B', mus0);

kap0 = 1./(3./(mua0+mus0));
% Mesh has been made and properties populated
%% Transfer from node basis to element basis
% mua0 = zeros(ne,1);
% mus0 = zeros(ne,1);
% coords = zeros(ne,2);
% 
% for i = 1:ne
%     coord = (h0.Element(i).Data());
%     coords(i,:) = mean(coord);
%     temp = h0.Element(i).ShapeFun(coords(i,:)', 'global');
%     mua0(i) = ([mua_node(idx(i,1)), mua_node(idx(i,2)), mua_node(idx(i,3))]*temp);
%     %mus0(i) = ([mus_node(idx(i,1)), mus_node(idx(i,2)), mus_node(idx(i,3))]*temp);
% end
% 
% %mua0 = ones(ne,1)*mua_bkg;
% mus0 = ones(ne,1)*mus_bkg;
% kap0 = 1./(3.*(mua0+mus0));
% 
% figure(1); title('Recovery Parameters');
% subplot(3,1,1); h0.Display(mua0); title('\mu_a on coarse mesh');
% subplot(3,1,2); h0.Display(mus0); title('\mu_s on coarse mesh');
% subplot(3,1,3); h0.Display(kap0); title('\kappa on coarse mesh');
% 
% fprintf('Inverse mesh populated and parameters stored');
%% Generate exact data from fine mesh
dataMesh = toastMesh('C:\Users\shrey\Desktop\TOAST\toast\toast\Adaptive-Meshing\cyl2.msh');
dataMesh.ReadQM('cyl_5ring.qm');
Q_dat = real(h0.Qvec('Neumann','Gaussian',2));
mvec_dat = real(h0.Mvec('Gaussian',2,refind));

basis_dat = toastBasis(dataMesh,grd_img);
mua_dat = basis_dat.Map('B->M',bmua);
%mus_dat = basis_dat.Map('B->M',bmus);
mus_dat = mus_bkg*ones(dataMesh.NodeCount(),1);
ref_dat = refind*ones(dataMesh.NodeCount(),1);

mua_homog = mua_bkg*ones(dataMesh.NodeCount(),1);

Q_dat = real(dataMesh.Qvec('Neumann','Gaussian',2));
mvec_dat = real(dataMesh.Mvec('Gaussian',2,refind));

% Generate exact data
smat = dotSysmat(dataMesh, mua_dat, mus_dat, ref_dat);
phi = smat \Q_dat;
y0_inhomog = reshape(log(mvec_dat'*(phi)), [], 1);
%y0_inhomog  = y0_inhomog + 0.025.*randn(size(y0_inhomog));


% Generate homogenous data
smat2 = dotSysmat(dataMesh, mua_homog, mus_dat, ref_dat);
phi2 = smat2 \Q_dat;
y0_homog = reshape(log(mvec_dat'*(phi2)), [], 1);

%% Refinement Loop
mua_r = ones(nh,1)*mua_bkg;
num_ref = 3;

for num_iter = 1:num_ref   %Number of refinements
    %% Initialize guess on inverse mesh
    
    mua_r_homog = ones(nh,1)*mua_bkg;
    %mus_r = ones(nh,1)*mus_bkg;
    %mus_r = mus_node;
    mus_r = ones(nh,1)*mus_bkg;
    ref = ones(nh,1)*refind;
   
    %mua_r = mua_r_homog;
    
    % Generate homogenous data on coarse mesh
    smat = dotSysmat(h0, mua_r_homog, mus_r, ref);
    y1_homog = reshape(log(mvec'*(smat\Q)), [], 1);
    data_model = y0_inhomog - y0_homog+y1_homog;
    %data_model = data_model + 0.1.*randn(size(data_model));
    
    % Generate improved guess on coarse mesh
    smat = dotSysmat(h0, mua_r, mus_r, ref);
    y1_inhomog = reshape(log(mvec'*(smat\Q)), [], 1);
    proj = y1_inhomog;
    sd = ones(size(proj));
    
    %Changing the inverse basis discretization
%     if(num_iter == 1)
%         grd_inv = [256, 256];
%     elseif(num_iter == 2)
%         grd_inv = [512, 512];
%     elseif(num_iter == 3)
%         grd_inv = [1024, 1024];
%     end
    
    grd_inv = [32, 32, 32];
    
    basis3 = toastBasis(h0, grd_inv);
    
    % Create solution vector
    bmua_r = basis3.Map('M->B', mua_r);
    bcmua = bmua_r*cm;
    scmua = basis3.Map('B->S', bcmua);
    x = scmua;
    logx = log(x);
    slen = length(x);
    
    % Set regularization parameters
    alpha = 1e-2;
    regul = toastRegul('TV', basis3, logx, alpha, 'Beta', 0.002);
    
    err0 = toastObjective(proj, data_model, sd, regul, logx);
    err = err0;
    errp = inf;
    itr = 1;
    step = 0.1;
    fprintf('Iteration %d, objective %f\n', 0, err);
    
    %% Inverse Solver Loop
    
    itrmax = 15; %CG iteration limit
    tolCG = 1e-6;
    
    while (itr <= itrmax) && (err > tolCG*err0) && (errp-err > tolCG)
        
        errp = err;
        
        r = -toastGradient(h0, basis3, Q, mvec, mua_r, mus_r, ref, 0, ...
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
        
        % Line search along update direction
        step = toastLineSearch(logx, d, step, err, @objective);
        logx = logx + d*step; %Update solution estimate
        
        mua_r = basis3.Map('S->M', exp(logx)/cm);
        
        proj = reshape(log(mvec'*(dotSysmat(h0, mua_r, mus_r, ref)\Q)), [], 1);
        
        err = toastObjective(proj, data_model, sd, regul, logx);
        
        fprintf('Iteration %d, objective %f\n', itr, err);
        itr = itr+1;
        
        figure(50);
    end
    h0.Display(mua_r);
    
    %% Refine h0 depending on the recovered mua_r
    
    
%     
%     
%     
%     for i = 1:ne
%         coord = (h0.Element(i).Data());
%         %x_coord = mean(coord(:,1));
%         %y_coord = mean(coord(:,2));
%         coords(i,:) = mean(coord);
%         temp = h0.Element(i).ShapeFun(coords(i,:)', 'global');
%         mua0(i) = ([mua_r(idx(i,1)), mua_r(idx(i,2)), mua_r(idx(i,3))]*temp);
%         mus0(i) = ([mus_r(idx(i,1)), mus_r(idx(i,2)), mus_r(idx(i,3))]*temp);
%     end
%     kap0 = 1./(3.*(mua0+mus0));
    neighbors = h0.ElementNeighbours();
    
    smat3 = dotSysmat(h0, mua_r, mus_r, ref);
    phi_res = smat3\Q;
    
    figure(3); h0.Display(mua0);
    basis4 = toastBasis(h0, grd_img);
    bmua_rec = basis4.Map('M->B',mua_r);
    
    
    [res, res_nh] = calc_vol_residue_3D(h0, phi_res, idx, mua_r, mus_r, Q);
    
    figure(2); h0.Display(res);
    
    res_inv = calc_residue_bangti(h0, phi_res, idx, kap0, neighbors, mua0, mua_r, ne, Q, vert);
    
    figure(2); h0.Display(res_inv); title('Inverse Residue');
    %% Set h0 as the refined mesh
    
    if(num_iter ~= num_ref)
        %Choose the top percent_residue of residues to refine
        res_sort = zeros(ne,2);
        for i = 1:ne
            res_sort(i,1) = res_inv(i);
            res_sort(i,2) = i;
        end
        
        res_sort2 = sortrows(res_sort,1);
        small = 0;
        percent_residue = 40;
        num = (percent_residue/10)*round(ne/10);
        
        re = zeros(num,1);
        for p = 1:1:num
            re(p) = res_sort2(ne-p,2);
        end
        
        [h2] = refine_mesh_sierpinski(vert, idx, re, mua0, neighbors);
        
        figure(10); h2.Display();
        
        
        h0 = h2;
        ne = h0.ElementCount;
        nh = h0.NodeCount;
        [vert,idx] = h0.Data();
        ref = ones(ne,1)*refind;
        
        h0.SetQM(qpos,mpos);
        Q = real(h0.Qvec('Neumann','Gaussian',2));
        mvec = real(h0.Mvec('Gaussian',2,refind));
        nQ = nopt;
        
        % Store parameters in nodal basis on h0
        basis_img = toastBasis(h0,grd_img);
        mua_node = basis_img.Map('B->M',bmua_rec);
        mua_r = mua_node;
        %mus_node = basis_img.Map('B->M',bmus);
        mus_node = mus_bkg*ones(nh,1);
    end
    % Mesh has been made and properties populated
    
end

    function p = objective(logx)
        mua_ = basis3.Map('S->M',(exp(logx)/cm));
        proj_ = reshape(log(mvec' * (dotSysmat(h0, mua_, mus_r, ref)\Q)), [], 1);
        p = toastObjective(proj_, data_model, sd, regul, logx);
    end

end