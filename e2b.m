% PROBLEM SETUP
nelx = 200;          % Number of elements in x-direction
nely = 80;           % Number of elements in y-direction
volfrac = 1;         % Maximum volume fraction
penal = 3;           % Penalization power
rmin = 1;            % Filter radius in terms of elements

% MATERIAL PROPERTIES
E0 = 210e9;          % Young's modulus for the solid material
Emin = 1e-9*E0;      % Young's modulus for the void material
nu = 0.3;            % Poisson's ratio
t = 0.01;            % Thickness
l = 0.01;            % Element length

% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = t/(1-nu^2)/24*([A11 A12; A12' A11] + nu*[B11 B12; B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1, nelx*nely, 1);
edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2*nely+[2 3 0 1] -2 -1], nelx*nely, 1);
iK = reshape(kron(edofMat, ones(8,1))', 64*nelx*nely, 1);
jK = reshape(kron(edofMat, ones(1,8))', 64*nelx*nely, 1);

% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
iF = nelx*(2*nely+2)+[2:2:nely/2]; 
F = sparse(iF,1,-50e3/length(iF),2*(nely+1)*(nelx+1),1);

fixeddofs = 1:2*(nely+1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs, fixeddofs);

% Initialize Values for the Cycle
U = zeros(2*(nely+1)*(nelx+1), 1);
U3 = zeros(2*(nely+1)*(nelx+1), 1);
U4 = zeros(2*(nely+1)*(nelx+1), 1);
Lambda = zeros(2*(nely+1)*(nelx+1), 1);
sig_xxe = zeros(nelx*nely, 1);
sig_yye = zeros(nelx*nely, 1);
sig_xye = zeros(nelx*nely, 1);
qv = zeros(2*(nely+1)*(nelx+1), nelx*nely); 

% STRESS CALCULATION - CONSTRAINT
sig_xxe3 = zeros(nelx*nely, 1);
sig_yye3 = zeros(nelx*nely, 1);
sig_xye3 = zeros(nelx*nely, 1);
sig_vMe3 = zeros(nelx*nely, 1);

% FORMULATION OF THE OPTIMIZATION PROBLEM
s_max = 235e6;
q = 2.5;
r = 2;
delta_r = 0.1;
Be = 1/(2*l)*[-1, 0, 1, 0, 1, 0, -1, 0; 0, -1, 0, -1, 0, 1, 0, 1; -1, -1, -1, 1, 1, 1, 1, -1];
D0 = E0/(1-nu^2)*[1, nu, 0; nu, 1, 0; 0, 0, (1-nu)/2];

% PREPARE FILTER
[dy, dx] = meshgrid(-ceil(rmin)+1:ceil(rmin)-1, -ceil(rmin)+1:ceil(rmin)-1);
h = max(0, rmin-sqrt(dx.^2 + dy.^2));
Hs = conv2(ones(nely, nelx), h, 'same'); 

% DESIGN VARIABLES
x = ones(nely, nelx);   % Initialization of design variables 
xmin = 0;               % Minimum value design variables
xmax = 1;               % Maximum value design variables

% OPTIMIZATION LOOP
iter = 1;               % Iteration counter
change = 1;             % (Fictitious) initial design change
tol = 0.02;             % Convergence threshold
mmaparams = [];         % Internal MMA parameters
move = 0.05;            % Move parameter for MMA

% FINITE DIFFERENCE
h3 = 1e-9;
df3dx = zeros(nelx*nely, 1);

while change > tol

    % APPLY FILTER TO OBTAIN PHYSICAL DENSITIES
    xPhys = conv2(x, h, 'same')./Hs;
    xPhys = xPhys(:);

    % FINITE ELEMENT ANALYSIS
    sK = reshape(KE(:)*(Emin + xPhys(:)'.^penal * (E0 - Emin)), 64*nelx*nely, 1);
    K = sparse(iK, jK, sK);
    K = (K + K') / 2;
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
    
    % OBJECTIVE FUNCTION
    v = sum(xPhys(:)) / nelx / nely;
    f0 = v - 1;

    % STRESS CALCULATION - CONSTRAINT
    
    for el = 1:nelx*nely
        Ue = U(edofMat(el, :));
        sig_xxe(el) = D0(1, :) * Be * Ue;
        sig_yye(el) = D0(2, :) * Be * Ue;
        sig_xye(el) = D0(3, :) * Be * Ue;
    end
    s_vMe = sqrt(sig_xxe.^2 + sig_yye.^2 - sig_xxe .* sig_yye + 3 .* sig_xye.^2);
    f = max(s_vMe./xPhys(:).^(q-penal)/s_max)-1;

    % SENSITIVITIES
    
    for el = 1:nelx*nely
        Ce = sparse(1:8, edofMat(el, :), ones(1,8), 8, 2*(nely+1)*(nelx+1));
        q1 = ((s_vMe(el) / ((xPhys(el)^(q-penal))* s_max))^(r-1)) * (1/(xPhys(el)^(q-penal)*s_max)) * (1/(2*s_vMe(el)));
        qe = q1 * (((2*sig_xxe(el) - sig_yye(el)) * D0(1,:) + (2*sig_yye(el) - sig_xxe(el)) * D0(2,:) + 6*sig_xye(el) * D0(3,:)) * Be * Ce)';
        qv(:, el) = qe;
    end
    qs = sum(qv, 2);
    Lambda(freedofs) = K(freedofs, freedofs) \ qs(freedofs);   

    term22 = sum((Lambda(edofMat)*KE).*U(edofMat),2);
    term22 = (penal*(xPhys.^(penal-1) * (E0 - Emin))).*term22;
    
    term1 = (sum((s_vMe ./ ((xPhys.^(q-penal))* s_max)).^r))^((1/r) -1);
    term21 = ((s_vMe ./ ((xPhys.^(q-penal))* s_max)).^(r-1)) * (penal - q) .* (s_vMe ./ (xPhys.^(q-penal+1) * s_max));
    
    dpdxPhys = term1 * (term21 - term22);
    
    dpdxPhys = reshape(dpdxPhys, nely, nelx);    % Sensitivity wrt physical densities 
    dvdxPhys = repmat(1/nelx/nely, nely, nelx);  % Sensitivity wrt physical densities 
    dpdx = conv2(dpdxPhys ./ Hs, h, 'same');     % Sensitivity wrt design variables
    dvdx = conv2(dvdxPhys ./ Hs, h, 'same');     % Sensitivity wrt design variables
    
    df0dx = dvdx(:);
    dfdx = dpdx(:);
    
    if iter == 10

    % FINITE DIFFERENCE CHECK
    x3 = x + h3;
    xPhys3 = conv2(x3, h, 'same')./Hs;
    xPhys3 = xPhys3(:);
    
    sK3 = reshape(KE(:)*(Emin + xPhys3(:)'.^penal * (E0 - Emin)), 64*nelx*nely, 1);
    K3 = sparse(iK, jK, sK3);
    K3 = (K3 + K3') / 2;
    U3(freedofs) = K3(freedofs, freedofs) \ F(freedofs);
    U3(fixeddofs) = 0;
    
    for el3 = 1:nelx*nely
        Ue3 = U3(edofMat(el3, :));
        sig_xxe3(el3) = D0(1, :) * Be * Ue3;
        sig_yye3(el3) = D0(2, :) * Be * Ue3;
        sig_xye3(el3) = D0(3, :) * Be * Ue3;
        sig_vMe3(el3) = sqrt(sig_xxe3(el3)^2 + sig_yye3(el3)^2 - sig_xxe3(el3)* sig_yye3(el3) + 3* sig_xye3(el3)^2);
        f3 = (sig_vMe3(el3) / ((xPhys3(el3)^(q-penal))* s_max)) - 1;
        f2 = (s_vMe(el3)/ ((xPhys(el3)^(q-penal))* s_max)) - 1;
        df3dx(el3, :) = (f3-f2)/(h3);
    end
    err = max(df3dx - dfdx)
    end

  % MMA UPDATE
  [xnew,~,~,~,mmaparams,subp,change,history] = mma(x(:), xmin, xmax, f0, f, df0dx, dfdx, mmaparams, move);

  % COLLECT CONVERGENCE HISTORY
  history.iter(:,iter) = iter;
  history.x(:,iter) = x(:);
  history.f0(:,iter) = f0;
  history.f(:,iter) = f;
  
  % PLOT CURRENT DESIGN
  xPhys = reshape(xPhys, nely, nelx);
  figure(1)
  colormap(gray);
  imagesc(1-xPhys);
  caxis([0 1]);
  axis('equal');
  axis('off');
  drawnow;
  
  s_vMe = reshape(s_vMe, nely, nelx);
  sig_vMem = s_vMe .* ((xPhys.^penal).*s_max);
  figure(2)
  imagesc(sig_vMem);
  axis('equal');
  axis('off');
  colormap('turbo');
  colorbar;
  s_vMe = s_vMe(:);

  % MOVE ON TO NEXT ITERATION
  iter = iter+1;
  x(:) = xnew;
  x = reshape(x, nely, nelx);
  r = r + delta_r;
  if r > 40
      r = 40;
  end
end

% PLOT CONVERGENCE HISTORY
figure;
plot(history.funevals,history.f0);
xlabel('Function evaluations');
ylabel('Objective function');
figure;
plot(history.funevals,history.f);
xlabel('Function evaluations');
ylabel('Constraint function');