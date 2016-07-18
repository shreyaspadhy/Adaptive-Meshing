
%% make a simple DOT problem
meshpath = '/home/spadhy/Desktop/toastpp/test/2D/meshes/';
meshname = 'circle25_32.msh';
qmname = 'circle25_32x32.qm';        % QM file
refind = 1.4;                        % refractive index
bx =64; by = 64;                    % solution basis: grid dimension
freq = 0;                            % modulation frequency [MHz]
noiselevel = 2e-2;                   % multiplicative noise level
READIM = false; %true;%false; %
wname = 'haar';%'db4';%'coif4';
Nw = 5;
ExplicitWM = false;%true;            % this builds wavelets as a matrix

WST = true;                          % use softthresholding on wavelets
WHT = false;                         % use hard thresholding on wavelets
SGD = false;                         % Steepest Descent
ISTA= true;                         % "normal" (i.e. pixel) ISTA
FISTA= true;                         % "normal" (i.e. pixel) FISTA
ADMM= true;                         % "normal" (i.e. pixel) ADMM
%wf = 2e-2;                           % relative threshold for wavelet compression
DIRECT_PRIMAL_SOLVER = false;        % solve ADMM primal step with backslash
NISTAit  = 1000;
NFISTAit = 200;
NADMMit = 10;
LSQRit = 40; % solve LSQR with small number of iterations
LSQRtol = 1e-6; % tolerence of LSQR
% ISTA FLAGS
ISTA_FLAGS.FISTA= true;             
ISTA_FLAGS.pos = true;
ISTA_FLAGS.Wavelets =true;
ISTA_FLAGS.Iterates = true;
% ======================================================================
% End user-defined parameters
% ======================================================================

%% Initialisations
%toastCatchErrors();

% Set up some variables
blen = bx*by;
c0 = 0.3;
cm = c0/refind;

% Read a TOAST mesh definition from file.
hMesh = toastMesh ([meshpath meshname]);
hMesh.ReadQM ([meshpath qmname]);
n = hMesh.NodeCount ();
dmask = hMesh.DataLinkList();

%% Set up homogeneous initial parameter estimates
mua = ones(n,1) * 0.025;
mus = ones(n,1) * 2;
ref = ones(n,1) * refind;
kap = 1./(3*(mua+mus));

% Set up the mapper between FEM and solution bases
hBasis = toastBasis(hMesh,[bx,by]);
%solmask = hBasis.SolutionMask ();
%solmask2 = [solmask solmask+blen];

% Generate source vectors
qvec = hMesh.Qvec ('Neumann', 'Gaussian', 2);

% Generate measurement vectors
mvec = hMesh.Mvec ('Gaussian', 2);

%% Initial data set f[x0]
smat = dotSysmat (hMesh,mua, mus, ref, freq);
lgamma = reshape (log(mvec.' * (smat\qvec)), [], 1);
lgamma = lgamma(dmask);
proj = [real(lgamma);imag(lgamma)];

% data scaling
sd_lnmod = ones(size(lgamma));% * norm(mdata-real(lgamma)); 
sd_phase = ones(size(lgamma));% * norm(pdata-imag(lgamma));
sd = [sd_lnmod;sd_phase];

%% initial parameter estimates in solution basis
bmua = hBasis.Map('M->B', mua);
bmus = hBasis.Map('M->B', mus);
bmua_itr(1,:) = bmua;
bmus_itr(1,:) = bmus;
bkap = hBasis.Map('M->B', kap);
bcmua = bmua*cm;
bckap = bkap*cm;
%scmua = bcmua(solmask);
%sckap = bckap(solmask);
scmua = hBasis.Map('B->S', bcmua);
sckap = hBasis.Map('B->S', bckap);


x = [scmua;sckap];
logx = log(x);
p = length(x);
step = 1.0; % initial step length for line search



%%
fprintf (1,'Calculating Jacobian\n');

A = toastJacobian (hMesh, hBasis, qvec, mvec, mua, mus, ref, freq, 'direct');
m = size(A,1); n = size(A,2);
ndat = m/2; nsol = n/2;
Ja = A(1:ndat,1:nsol); % this is the "DC-mua" part of the problem.
%Ja = diag(1./exp(lgamma)) * Ja; % data rescaling; 
%% observations 
tgtmuaim = zeros(bx,by);
%
if READIM
    im = imresize(imread('ot.pgm','pgm'),[bx,by])/255;
    tgtmuaim(floor(0.27*bx):ceil(0.73*bx),floor(0.27*by):ceil(0.73*by)) = im(floor(0.27*bx):ceil(0.73*bx),floor(0.27*by):ceil(0.73*by));
else
    tgtmuaim(floor(0.3*bx):ceil(0.7*bx),floor(0.3*by):ceil(0.48*by)) = 1;
    tgtmuaim(floor(0.4*bx):ceil(0.6*bx),floor(0.48*by):ceil(0.7*by)) = 1;
end
x1a = hBasis.Map('B->S',tgtmuaim(:));

y = Ja*x1a;
y = (1 + noiselevel*randn(size(y))).*y;

lambda = sum(sum(Ja.^2));
alpha = noiselevel^2*lambda;
tau = 2/lambda;
T = alpha*tau;    % threshold for SoftThresholding
%% Initialise regularisation
ttau = alpha;
beta = 1e-2;
%hReg = toastRegul ('LAPLACIAN', hBasis, logx, tau);
%hReg = toastRegul('TK1',hBasis,logx,tau);
hReg = toastRegul ('TV', hBasis, logx, ttau, 'Beta', beta);
%% initial guess - min energy

tic;
%x0 = pcg(Ja'*Ja + alpha*speye(nsol),Ja'*y,1e-6,100);   
x0 = pcg(@(x) (Ja' * (Ja*x) + alpha*x),Ja'*y,1e-6,100);   
toc;
disp('solved using pcg');

figure(2);clf
subplot(3,3,1);imagesc(tgtmuaim);colorbar;title('\mu_a target');

x1im = reshape(hBasis.Map('S->B',x0(1:nsol)),bx,by);
subplot(3,3,3);imagesc(x1im);colorbar;title('\mu_a solution from pcg');

%% now gradient descent (slow!)
if(SGD)
    xsd = zeros(size(x0));
    Niter = 5000;
    tic
    for k = 1:Niter
        xsd = xsd + tau*Ja'*(y - Ja*xsd);
        if ~mod(k,100)
            disp(['SGD : iteration ', num2str(k)]);
        end
    end
    toc;
    disp('solved using gradient descent');
    
    xsdim = reshape(hBasis.Map('S->B',xsd(1:nsol)),bx,by);
    figure(2);
    subplot(3,3,7);imagesc(xsdim);colorbar;title('Landweber');
end;

%% set up operators for Shrinkage - pixel case

pops.W    = @(x) x;
pops.Winv = @(c) c;
pops.WTr  = @(c) c;
pops.B    = @(x)x;
pops.Binv = @(x)x;
pops.A    = @(x)Ja*x;
pops.ATr  = @(y) Ja'*y;
pops.ntc  = nsol;

pops.ndat = ndat;
pops.nsol = nsol;
pops.Dims = hBasis.Dims;

%% now ISTA (slow!)
if (ISTA)
    xista = zeros(nsol,1);
    Niter = NISTAit;
 %   T = alpha*tau;
    ISTA_FLAGS.Rayleighstep = true;
    ISTA_FLAGS.pos = false;
    ISTA_FLAGS.FISTA = false;
    ISTA_FLAGS.Wavelets = false;
    ISTA_FLAGS.Iterates = true;
    tic;
%    [x,c,postista] = DotWavDeconvByISTA(hBasis,Ja,y,alpha,Niter,tau,ISTA_FLAGS);
    [x,c,postista] = LinearShrinkage(pops,y,alpha,NISTAit,tau,ISTA_FLAGS);
    toc;
    disp('--------------  solved using ISTA  --------------');
    xista = x(:,end);
    for k = 1:NISTAit
        lerrista(k) = norm(y - Ja*x(:,k));
        xistim = reshape(hBasis.Map('S->B',x(:,k)),bx,by);
        xerrista(k) = norm(xistim-tgtmuaim );
        perrista(k) = norm(c(:,k),1);        
    end
    figure(2);
    subplot(3,3,4);imagesc(xistim);colorbar;title('ISTA');
end
%% now FISTA (faster?)
if (FISTA)
    Niter = NFISTAit;
    ISTA_FLAGS.pos = false;
    ISTA_FLAGS.FISTA = true;
    ISTA_FLAGS.Wavelets = false;
    ISTA_FLAGS.Iterates = true;
    tic;
%    [x,c,postfista] = DotWavDeconvByISTA(hBasis,Ja,y,alpha,Niter,tau,ISTA_FLAGS);
    [x,c,postfista] = LinearShrinkage(pops,y,alpha,NFISTAit,tau,ISTA_FLAGS);
    toc;
    disp('--------------  solved using FISTA  --------------');
    xfistay = x(:,end);
    for k = 1:NFISTAit
        lerrfista(k) = norm(y - Ja*x(:,k));
        xfistim = reshape(hBasis.Map('S->B',x(:,k)),bx,by);
        xerrfista(k) = norm(xfistim-tgtmuaim );
        perrfista(k) = norm(c(:,k),1);        
    end
  
    figure(2);
    subplot(3,3,5);imagesc(xfistim);colorbar;title('FISTA');
end

%% ADMM - pixels
if ADMM

rho = 1e-2*lambda; % 1e-2*alpha/T; % what is best way to set rho ?
sqrho = sqrt(rho);
tic; 
[xpadm,padm,postpadm] = LinearADMM(pops,y,alpha,NADMMit,tau,rho,LSQRtol,LSQRit,ISTA_FLAGS);
toc;
for k = 1:NADMMit
    lerrpadm(k) = norm(y - Ja*xpadm(:,k));
    xpadmmim = reshape(hBasis.Map('S->B',xpadm(:,k)),bx,by);
    xerrpadm(k) = norm(xpadmmim-tgtmuaim );
    perrpadm(k) = norm(padm(:,k),1);
end
disp('------------- solved using ADMM  -------------- ');

figure(2);
subplot(3,3,6);imagesc(xpadmmim);colorbar;title('ADMM');
end
%% wavelets the hard way...
if ExplicitWM
[lnCT,lnS] = wavedec2(tgtmuaim,Nw,wname);
nwc = length(lnCT); % number of wavelet coefficients;
npix = bx*by;
disp('Building explicit wavelet matrix');
WM = zeros(nwc,nsol);
tic;
for j = 1:nsol
    xp = zeros(nsol,1);
    xp(j) = 1;
    
    WM(:,j) = wavedec2(reshape(hBasis.Map('S->B',xp),bx,by),Nw,wname);
end
toc;
disp('Inverting Wavelet Matrix');
tic;
WMI = pinv(WM);
toc;
else
    if ~strcmp(wname,'haar')
        disp('**********************************************');
        disp(['Wavelet filter ',wname,'-- Need to set ExplicitWM to true']);
        wname = 'haar';
        disp(['switching to ',wname, ' basis']);
        disp('**********************************************');
    end
end
%% Sanity check wavelets

if ExplicitWM
    lnCT = WM*hBasis.Map('B->S',tgtmuaim);
else
    %[lnC,lnS] = wavedec2(tgtmuaim,Nw,wname);
    [lnCT,lnS] = haardec2(tgtmuaim);
end
lnca4 = lnS(1,:);
%T = wf*max(abs(lnCT(lnca4(1)*lnca4(2) + 1:end))); 
nwc = length(lnCT); % number of wavelet coefficients;

if WHT
    disp('HardThresholding');
    lnCr = HardThresh(lnCT,T)   ;
else if WST
        disp('SoftThresholding');
        lnCr = SoftThresh(lnCT,T);
    else
        lnCr = lnCT;
        zc = find(abs(lnCr(lnca4(1)*lnca4(2) + 1:end)) < T);
        cpr = num2str(size(zc,2)*100/(bx*by));
        disp(['setting ',num2str(size(zc,2)),' out of ',num2str(bx*by),' wavelets to zero (',cpr,'%)']);
        lnCr(lnca4(1)*lnca4(2) + zc) = 0;
    end
end
if ExplicitWM
%    tgtrec = reshape(hBasis.Map('S->B',WM\lnCr),bx,by);
     tgtrec = waverec2(lnCr,lnS,wname);

else
    %tgtrec = waverec2(lnCr,lnS,wname);
    tgtrec = haarrec2(lnCr);
end
figure(2);
subplot(3,3,2);imagesc(tgtrec);colorbar;title([wname,'-wavelet compressed']);

%% set up operators for Shrinkage - wavelet case
%wops.hBasis = hBasis; % This is for DOT. Would not be needed in other problems
if ExplicitWM
    wops.W    = @(x) WM*x;
    wops.Winv = @(c)WMI*c;
    wops.WTr  = @(c)WM'*c;
    wops.B    = @(x)x;
    wops.Binv = @(x)x;
    wops.A    = @(x)Ja*x;
    wops.ATr  = @(y) Ja'*y;
    
else
    wops.W    = @(x) haardec2(reshape(x,bx,by));
    wops.Winv = @(c)haarrec2(c);
    wops.WTr  = @(c)haarrec2(c);
    wops.B    = @(x)hBasis.Map('B->S',x);
    wops.Binv = @(x)hBasis.Map('S->B',x);
    wops.A    = @(x)Ja*hBasis.Map('B->S',x);
    wops.ATr  = @(y)reshape(hBasis.Map('S->B',Ja'*y),bx,by);
end

wops.ndat = ndat;
wops.nsol = nsol;
wops.Dims = hBasis.Dims;
wops.ntc  = nwc;
%% Wavelet ISTA
%Niter = 100;
ISTA_FLAGS.pos = true;
ISTA_FLAGS.FISTA = false;
ISTA_FLAGS.Wavelets = true;
ISTA_FLAGS.Iterates = true;
tic;
%[x,c,postawi] = DotWavDeconvByISTA(hBasis,Ja,y,alpha,Niter,tau,ISTA_FLAGS);
[x,c,postawi] = LinearShrinkage(wops,y,alpha,NISTAit,tau,ISTA_FLAGS);
 toc;
disp('-------------- solved using Wavelet-ISTA -------------- ');
xwi = x(:,end);
for k = 1:NISTAit
    lerrawi(k) = norm(y - Ja*x(:,k));
    xwim = reshape(hBasis.Map('S->B',x(:,k)),bx,by);
    xerrawi(k) = norm(xwim-tgtmuaim );
    perrawi(k) = norm(c(:,k),1);
end

figure(2);
subplot(3,3,7);imagesc(xwim);colorbar;title('Wavelet-ISTA');

%% Wavelet FISTA 
%Niter = 100;
ISTA_FLAGS.pos = true;
ISTA_FLAGS.FISTA = true;
ISTA_FLAGS.Wavelets = true;
ISTA_FLAGS.Iterates = true;
tic;
%[x,c,postawfi] = DotWavDeconvByISTA(hBasis,Ja,y,alpha,Niter,tau,ISTA_FLAGS);
[x,c,postawfi] = LinearShrinkage(wops,y,alpha,NFISTAit,tau,ISTA_FLAGS);
toc;

xwfiy = x(:,end);
for k = 1:NFISTAit
    lerrawfi(k) = norm(y - Ja*x(:,k));
    xwfim = reshape(hBasis.Map('S->B',x(:,k)),bx,by);
    xerrawfi(k) = norm(xwfim-tgtmuaim );
    perrawfi(k) = norm(c(:,k),1);
end

disp('-------------- solved using Wavelet-FISTA  -------------- ');

xwfim = reshape(hBasis.Map('S->B',xwfiy(1:nsol)),bx,by);
figure(2);
subplot(3,3,8);imagesc(xwfim);colorbar;title('Wavelet-FISTA');
%% wavelet ADMM

rho = 1e-2*lambda; % 1e-2*alpha/T; % what is best way to set rho ?
sqrho = sqrt(rho);
tic; 
[xwadm,wadm,postwadm] = LinearADMM(wops,y,alpha,NADMMit,tau,rho,LSQRtol,LSQRit,ISTA_FLAGS);
toc;
for k = 1:NADMMit
    lerrwadm(k) = norm(y - Ja*xwadm(:,k));
    xwadmmim = reshape(hBasis.Map('S->B',xwadm(:,k)),bx,by);
    xerrwadm(k) = norm(xwadmmim-tgtmuaim );
    perrwadm(k) = norm(wadm(:,k),1);
end
disp('------------- solved using Wavelet-ADMM  -------------- ');

figure(2);
subplot(3,3,9);imagesc(xwadmmim);colorbar;title('Wavelet-ADMM');

%%
xlin = ceil(bx/2);
ylin = ceil(by/2);
figure(1); clf; hold on; 
plot( tgtmuaim(xlin,:),'k');
plot( x1im(xlin,:),'b');
if SGD
    plot( xsdim(xlin,:),'+b');
end
if ISTA
    plot( xistim(xlin,:),'og');
end
if FISTA
    plot( xfistim(xlin,:),'om');
end
if ADMM
    plot( xpadmmim(xlin,:),'or');
end
plot( xwim(xlin,:),'+g');
plot( xwfim(xlin,:),'+m');
plot( xwadmmim(xlin,:),'+r');

%% error plots
% - posterior 
figure(3); clf; hold on; 
title('posterior - (log scale)')
plot(log(postwadm),'-r+'); 
nawi = length(postawi);
plot(log(postawi([1:ceil(NISTAit/NADMMit):end])),'-g+');
plot(log(postawfi([1:ceil(NFISTAit/NADMMit):end])),'-m+');
if ISTA
    nista = length(postista);
    plot(log(postista([1:ceil(nista/NADMMit):end])),'-go');
end
if FISTA
    nfista = length(postfista);
    plot(log(postfista([1:ceil(nfista/NADMMit):end])),'-mo');
end
if ADMM
plot(log(postpadm),'-ro'); 
end
% - image error
figure(4); clf; hold on; 
title('Image error');
plot(xerrwadm,'-r+');
plot(xerrawi([1:ceil(NISTAit/NADMMit):end]),'-g+');
plot(xerrawfi([1:ceil(NFISTAit/NADMMit):end]),'-m+');
if ISTA
    plot(xerrista([1:ceil(NISTAit/NADMMit):end]),'-go');
end
if ISTA
    plot(xerrfista([1:ceil(NFISTAit/NADMMit):end]),'-mo');
end
if ADMM
    plot(xerrpadm,'-ro');
end

%% clean up
%clear hMesh hBasis hReg;