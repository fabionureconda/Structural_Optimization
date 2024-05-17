% ELEMENTS 
nelx = 24;
nely = 60;
nel = nely*nelx;

% ELEMENT LENGTHS
l = 800/nely;

% PENALIZATION
penal = 3;

% NODES
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

% BOUNDARY CONDITIONS AND LOADS
F = sparse(2*(nely+1)*(nelx)+2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);

% MATERIALS 
E0 = 210*10^9;
Emin = 10^(-9)*E0;
nu = 0.3;
t = 10;
Ee = ones(nelx, nely);

% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];

% it is depending on the Ee
K_e = {(t*Ee)/((1-nu^2)*24)}*([A11 A12;A12' A11]+[B11 B12;B12' B11]);

nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

% FORMULATION OF THE OPTIMIZATION PROBLEM
% maximum allowable stress
sig_max = 235;
q = 2.5;
r = 10;
Be = 1/(2*l)*[-1, 0, 1, 0, 1, 0, -1, 0,; 0, -1, 0, -1, 0, 1, 0, 1; -1, -1, -1, 1, 1 1 1 -1];
D0 = E0*(1-nu^2)*[1, nu, 0; nu, 1, 0;0, 0, (1-nu)/2];
% what is the location matrix
Ce = locationmatrix;

x = repmat(volfrac,nely,nelx);
xmin = 1e-4;
xmax = 1;

% OPTIMIZATION LOOP
iter = 1;                         % Iteration counter
change = inf;                     % (Fictitious) initial design change
tol = 1e-5;                       % Convergence threshold
mmaparams = [];                   % Internal MMA parameters
while change>tol

  % FINITE ELEMENT ANALYSIS
  sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);

  % OBJECTIVE FUNCTION AND CONSTRAINT            
  
  % FOR CYCLE EVERY ELEMENT

  % elasticity matrix
  ue = Ce*U(freedofs);
  De(:,:) = Ee(:,:)/(1-nu^2)*[1 nu 0;
                              nu 1 0;
                              0  0 (1-nu)/2];
  % stress definitions
  sig_xxe = De(x)*Be(x)*ue(x);
  sig_yye = De(y)*Be(y)*ue(y);
  sig_xye = De(z)*Be(z)*ue(z);
  % Von Mises stress
  sig_vMe = sqrt(sig_xxe^2 + sig_yye^2-sig_xxe*sig_yye+3*sig_xye^2);
  % stresses
  sigma_micr = D0*Be*ue;
  sigma_macr = De*Be*ue;
  % Young
  Ee = (sigma_macr/sigma_micr)*E0;

  rho_e = Ee/E0;

  % SENSITIVITIES
  for 1:nelx && 1:nely
      term11 = (sig_vMe/rho_e^(q-p)*sig_max)^(r);
  end
  term1 = (term11())^(1/r-1);
  term2 = 

  % MMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change] = mma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

  % PLOT CURRENT DESIGN
  SectionsX = Sections;
  SectionsX(:,2) = x;
  plotelemsec(Nodes,Elements,Types,SectionsX,'Numbering','off','GCS','off');
  drawnow;

  % CHECK CONVERGENCE CRITERION AND MOVE ON TO NEXT ITERATION
  if change>tol
    iter = iter+1;
    x = xnew;
  end

end