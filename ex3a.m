%EX3A
%
%   Topology optimization of an MBB beam for minimum compliance using
%   sensitivity or density filtering in combination with MMA.

% Mattias Schevenels
% February 2021

% PROBLEM SETUP
nelx = 120;          % Number of elements in x-direction
nely = 40;           % Number of elements in y-direction
volfrac = 1;       % Maximum volume fraction
penal = 3;           % Penalization power
rmin = 2;          % Filter radius in terms of elements
ft = 2;              % Filter type: 1 for sensitivity, 2 for density

% MATERIAL PROPERTIES
E0 = 210*10^9;
Emin = 10^(-9)*E0;
nu = 0.3;
t = 0.01;

% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2*(nely+1)*(nelx)+2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);

fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

% PREPARE FILTER
[dy,dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1,-ceil(rmin)+1:ceil(rmin)-1);
h = max(0,rmin-sqrt(dx.^2+dy.^2));
Hs = conv2(ones(nely,nelx),h,'same');

% DESIGN VARIABLES
x = ones(nely,nelx);              % Initialization of design variables
xmin = 0;                         % Minimum value design variables
xmax = 1;                         % Maximum value design variables

% OPTIMIZATION LOOP
iter = 1;                         % Iteration counter
change = 1;                       % (Fictitious) initial design change
tol = 0.01;                       % Convergence threshold
mmaparams = [];                   % Internal MMA parameters
while change>tol

  % APPLY FILTER TO OBTAIN PHYSICAL DENSITIES
  if ft == 1
    xPhys = x;
  elseif ft == 2
    xPhys = conv2(x,h,'same')./Hs;
  end

  % FINITE ELEMENT ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);

  % OBJECTIVE FUNCTION AND CONSTRAINT
  Ce = sparse(1:8, edofMat(e, :), ones(1,8), 8, alldofs)
  ue = Ce()*U
  ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);

  ue = edofMat()*U()
  c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));


  v = sum(xPhys(:))/nelx/nely;
  f0 = c/100;
  f = v/volfrac-1;

  % SENSITIVITIES
  if ft == 1
    dcdx = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dvdx = repmat(1/nelx/nely,nely,nelx);
    dcdx = conv2(dcdx.*xPhys,h,'same')./Hs./max(1e-3,xPhys); % Apply sensitivity filter
  elseif ft == 2
    dcdxPhys = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;        % Sensitivity wrt physical densities
    dvdxPhys = repmat(1/nelx/nely,nely,nelx);                % Sensitivity wrt physical densities
    dcdx = conv2(dcdxPhys./Hs,h,'same');                     % Sensitivity wrt design variables
    dvdx = conv2(dvdxPhys./Hs,h,'same');                     % Sensitivity wrt design variables
  end
  df0dx = dcdx(:)'/100;
  dfdx = dvdx(:)'/volfrac;

  % MMA UPDATE
  [xnew,~,~,~,mmaparams,~,change,history] = mma(x(:),xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

  % PLOT CURRENT DESIGN
  colormap(gray);
  imagesc(1-xPhys);
  caxis([0 1]);
  axis('equal');
  axis('off');
  drawnow;

  % CHECK CONVERGENCE CRITERIA AND MOVE ON TO NEXT ITERATION
  if change>tol
    iter = iter+1;
    x(:) = xnew;
  end

end
