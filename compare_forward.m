h1 = toastMesh('circle_blob.msh', 'gmsh');
ne1 = h1.ElementCount;
nh1 = h1.NodeCount;

refind = 1.4;

rad = 25; % mesh radius [mm]
nopt = 10;
for i=1:nopt
  phiq = (i-1)/nopt*2*pi;
  qpos(i,:) = rad*[cos(phiq), sin(phiq)];
  phim = (i-0.5)/nopt*2*pi;
  mpos(i,:) = rad*[cos(phim), sin(phim)];
end

h1.SetQM(qpos, mpos);
Q1 = real(h1.Qvec('Neumann','Gaussian',2));
mvec1 = real(h1.Mvec('Gaussian',2,refind));

bmua = imread('demo_matlab_fwd2_mua.png');
bmus = imread('demo_matlab_fwd2_mus.png');
bmua = double(bmua)./255.*0.02 + 0.01;
bmus = double(bmus)./255.*1.0 + 1.0;
grd_img = size(bmua);

basis_img = toastBasis(h1,grd_img);
mua_1 = basis_img.Map('B->M',bmua);
mus_1 = basis_img.Map('B->M',bmus);
ref_1 = ones(nh1,1)*refind;

smat1 = dotSysmat(h1, mua_1, mus_1, ref_1, 0);

data1 = log(mvec1'*(smat1\Q1));

h2 = toastMesh('circle_blob.msh', 'gmsh');
ne2 = h2.ElementCount;
nh2 = h2.NodeCount;

h2.SetQM(qpos, mpos);
Q2 = real(h2.Qvec('Neumann','Gaussian',2));
mvec2 = real(h2.Mvec('Gaussian',2,refind));

basis_img_2 = toastBasis(h2, grd_img);
mua_2 = basis_img_2.Map('B->M', bmua);
mus_2 = basis_img_2.Map('B->M', bmus);
ref_2 = ones(nh2,1)*refind;

smat2 = dotSysmat(h2, mua_2, mus_2, ref_2, 0);

data2 = log(mvec2'*(smat2\Q2));