%EX1D
%
%   Minimum compliance of a ten-bar truss under two point loads using MMA.

% Mattias Schevenels
% February 2021

% ELEMENTS
nelx = 24;
nely = 60;
nel = nely*nelx;

% MATERIALS [iMat E]: E = 1e7 psi = 68.94e9 N/m^2
E0 = 210*10^9;
Emin = 10^(-9)*E0;
nu = 0.3;
t = 10;
Ee = ones(nelx, nely);

% ELEMENT LENGTHS
l = 800/nely;

% BOUNDARY CONDITIONS: 2D analysis; simply supported
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);

% LOADS: 2 point loads of 100 kips = 444822 N
F = sparse(2*(nely+1)*(nelx)+2,1,-1,2*(nely+1)*(nelx+1),1);

% DETERMINE CONTRIBUTION K0j OF EACH ELEMENT TO THE STIFFNESS MATRIX
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
K_e = {(t*E0)/((1-nu^2)*24)}*([A11 A12;A12' A11]+[B11 B12;B12' B11]);


% FORMULATION OF THE OPTIMIZATION PROBLEM
x = repmat(volfrac,nely,nelx);    % Initial design
xmin = 1e-4;                      % Minimum section area
xmax = 0.1;                       % Maximum section area

% OPTIMIZATION LOOP
iter = 1;                         % Iteration counter
change = inf;                     % (Fictitious) initial design change
tol = 1e-5;                       % Convergence threshold
mmaparams = [];                   % Internal MMA parameters
while change>tol

  % FINITE ELEMENT ANALYSIS
  K = zeros(alldofs,alldofs);
  for iel = 1:nel
    K = K + K_e(:,:,iel)*x(iel);
  end
  E_e=(((1-nu^2)*24*K)/t)\([A11 A12;A12' A11]+[B11 B12;B12' B11]);
  U = K\F;

  % OBJECTIVE FUNCTION AND CONSTRAINT
  C = F'*U;                       % Compliance
  V = L'*x;                       % Volume
  f0 = C/1e4;
  f = V/Vmax-1;

  % SENSITIVITIES
  dCdx = zeros(nElem,1);
  for iElem = 1:nElem
    dCdx(iElem) = -U'*K0(:,:,iElem)*U;
  end
  dVdx = L;
  df0dx = dCdx/1e4;
  dfdx = dVdx/Vmax;

  % MMA UPDATE
  [xnew,y,z,lmult,mmaparams,subp,change,history] = mma(x,xmin,xmax,f0,f,df0dx,dfdx,mmaparams);

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

% DETERMINE ELEMENT FORCES
Forces = elemforces(Nodes,Elements,Types,SectionsX,Materials,DOF,U);

% PLOT STRESSES
figure;
plotstress('snorm',Nodes,Elements,Types,SectionsX,Forces);