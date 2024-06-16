% PROBLEM SETUP
nelx = 60;          % Number of elements in x-direction
nely = 24;          % Number of elements in y-direction
volfrac = 1;        % Maximum volume fraction
penal = 3;          % Penalization power
rmin = 4.5;           % Filter radius in terms of elements
unitF = 10e3;       % Force divided 5 times (unit force)

% MATERIAL PROPERTIES
E0 = 210e9;         % Young's modulus for the solid material
Emin = 1e-9*E0;     % Young's modulus for the void material
nu = 0.3;           % Poisson's ratio
t = 0.01;           % Thickness
l = 0.01;           % Element length

% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];

KE = 1/(1-nu^2)/24*([A11 A12; A12' A11] + nu*[B11 B12; B12' B11]);

nodenrs = reshape(1:(1+nelx)*(1+nely), 1+nely, 1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1, nelx*nely, 1);
edofMat = repmat(edofVec, 1, 8) + repmat([0 1 2*nely+[2 3 0 1] -2 -1], nelx*nely, 1);
iK = reshape(kron(edofMat, ones(8,1))', 64*nelx*nely, 1);
jK = reshape(kron(edofMat, ones(1,8))', 64*nelx*nely, 1);

% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F1 = sparse(2*((nely+1)*(nelx)+1), 1, -unitF, 2*(nely+1)*(nelx+1), 1);
F2 = sparse(2*((nely+1)*(nelx)+2), 1, -unitF, 2*(nely+1)*(nelx+1), 1);
F3 = sparse(2*((nely+1)*(nelx)+3), 1, -unitF, 2*(nely+1)*(nelx+1), 1);
F4 = sparse(2*((nely+1)*(nelx)+4), 1, -unitF, 2*(nely+1)*(nelx+1), 1);
F5 = sparse(2*((nely+1)*(nelx)+5), 1, -unitF, 2*(nely+1)*(nelx+1), 1);
F = F1+F2+F3+F4+F5;

U = zeros(2*(nely+1)*(nelx+1), 1);
U3 = zeros(2*(nely+1)*(nelx+1), 1);

fixeddofs = 1:2*(nely+1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs, fixeddofs);

% FORMULATION OF THE OPTIMIZATION PROBLEM
sig_max = 235e6;
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
tol = 0.01;             % Convergence threshold
mmaparams = [];         % Internal MMA parameters
move = 0.05;            % Move parameter for MMA

% FINITE DIFFERENCE
h3 = 1e-7;
df3dx = zeros(nelx*nely, 1);

while (change > tol) && (r < 40)

    % APPLY FILTER TO OBTAIN PHYSICAL DENSITIES
    xPhys = conv2(x, h, 'same')./Hs;
    xPhys = xPhys(:);

    % FINITE ELEMENT ANALYSIS
    sK = reshape(KE(:)*(Emin + xPhys(:)'.^penal * (E0 - Emin)), 64*nelx*nely, 1);
    K = sparse(iK, jK, sK);
    K = (K + K') / 2;
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);
    U(fixeddofs) = 0;
    
    % OBJECTIVE FUNCTION
    v = sum(xPhys(:)) / nelx / nely;
    f0 = v - 1;

    % STRESS CALCULATION - CONSTRAINT
    sig_xxe = zeros(nelx*nely, 1);
    sig_yye = zeros(nelx*nely, 1);
    sig_xye = zeros(nelx*nely, 1);

    for el = 1:nelx*nely
        Ue = U(edofMat(el, :));
        sig_xxe(el) = D0(1, :) * Be * Ue;
        sig_yye(el) = D0(2, :) * Be * Ue;
        sig_xye(el) = D0(3, :) * Be * Ue;
    end
    
    sig_vMe = sqrt(sig_xxe.^2 + sig_yye.^2 - sig_xxe .* sig_yye + 3 .* sig_xye.^2);
    
    f = sum((sig_vMe ./ ((xPhys.^(q-penal)).* sig_max)).^r)^(1/r) - 1;

    cs = max(sig_vMe ./ ((xPhys.^(q-penal)).* sig_max)) / sum((sig_vMe ./ ((xPhys.^(q-penal)).* sig_max)).^r)^(1/r);

    % SENSITIVITIES
    qv = zeros(2*(nely+1)*(nelx+1), nelx*nely); 
    for el = 1:nelx*nely
        Ce = sparse(1:8, edofMat(el, :), ones(1,8), 8, 2*(nely+1)*(nelx+1));
        Ce = full(Ce);
        q1 = ((sig_vMe(el) / ((xPhys(el)^(q-penal))* sig_max))^(r-1)) * (1/(xPhys(el)^(q-penal)*sig_max)) * (1/(2*sig_vMe(el)));
        qe = q1 * (((2*sig_xxe(el) - sig_yye(el)) * D0(1,:) + (2*sig_xye(el) - sig_xxe(el)) * D0(2,:) + 6*sig_xye(el) * D0(3,:)) * Be * Ce)';
        qv(:, el) = qe;
    end
    qs = sum(qv, 2);
    Lambda = K \ qs;    
%{
    sK2 = reshape(KE(:)*penal*(xPhys(:)'.^(penal-1) * (E0 - Emin)), 64*nelx*nely, 1);
    K2 = sparse(iK, jK, sK2);
    K2 = (K2 + K2') / 2;
    
    ce = L1 * KE;
    term222 = (penal*(xPhys(:)'.^(penal-1) * (E0 - Emin))).*ce(:);
    %term22 = L1 * K22 * U;
%} 
    term22 = reshape(sum((Lambda(edofMat)*KE).*U(edofMat),2),nely,nelx);
    term22 = (penal*(xPhys(:).^(penal-1) * (E0 - Emin))).*term22(:);
    
    term1 = sum((sig_vMe ./ ((xPhys.^(q-penal)).* sig_max)).^r)^((1/r) -1);
    term21 = ((sig_vMe ./ ((xPhys.^(q-penal)).* sig_max)).^(r-1)) .* (penal - q) .* (sig_vMe ./ (xPhys.^(q-penal+1) .* sig_max));
    
    dpdxPhys = term1 * cs * (term21 - term22);
    
    dpdxPhys = reshape(dpdxPhys, nely, nelx);    % Sensitivity wrt physical densities 
    dvdxPhys = repmat(1/nelx/nely, nely, nelx);  % Sensitivity wrt physical densities 
    dpdx = conv2(dpdxPhys ./ Hs, h, 'same');     % Sensitivity wrt design variables
    dvdx = conv2(dvdxPhys ./ Hs, h, 'same');     % Sensitivity wrt design variables
    
    df0dx = dvdx(:);
    dfdx = dpdx(:);
    
    % FINITE DIFFERENCE CHECK
    x3 = x + h3;
    xPhys3 = conv2(x3, h, 'same')./Hs;
    xPhys3 = xPhys3(:);
    
    sK3 = reshape(KE(:)*(Emin + xPhys3(:)'.^penal * (E0 - Emin)), 64*nelx*nely, 1);
    K3 = sparse(iK, jK, sK3);
    K3 = (K3 + K3') / 2;
    U3(freedofs) = K3(freedofs, freedofs) \ F(freedofs);
    U3(fixeddofs) = 0;
    
    % STRESS CALCULATION - CONSTRAINT
    sig_xxe3 = zeros(nelx*nely, 1);
    sig_yye3 = zeros(nelx*nely, 1);
    sig_xye3 = zeros(nelx*nely, 1);
    
    for el3 = 1:nelx*nely
        Ue3 = U3(edofMat(el3, :));
        sig_xxe3(el3) = D0(1, :) * Be * Ue3;
        sig_yye3(el3) = D0(2, :) * Be * Ue3;
        sig_xye3(el3) = D0(3, :) * Be * Ue3;
    end
    
    sig_vMe3 = sqrt(sig_xxe3.^2 + sig_yye3.^2 - sig_xxe3 .* sig_yye3 + 3 .* sig_xye3.^2);

    for el = 1:nelx*nely
        f3 = (sig_vMe3(el) / ((xPhys3(el)^(q-penal))* sig_max)) - 1;
        f2 = (sig_vMe(el)/ ((xPhys(el)^(q-penal)).* sig_max)) - 1;
        df3dx(el, :) = (f3-f2)/(h3);
    end

  % MMA UPDATE
  [xnew,~,~,~,mmaparams,subp,change,history] = mma(x(:), xmin, xmax, f0, f, df0dx, df3dx, mmaparams, move);

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
  
  sig_vMe = reshape(sig_vMe, nely, nelx);
  
  figure(2)
  imagesc(sig_vMe);
  axis('equal');
  axis('off');
  colormap('turbo');
  colorbar;

  sig_vMe = sig_vMe(:);
  %{
  figure(3)
  plot(dfdx', 'r')
  hold on
  plot(df3dx, 'b')
  hold off
  %}
  % CHECK CONVERGENCE CRITERIA AND MOVE ON TO NEXT ITERATION
  if change>tol
    iter = iter+1;
    x(:) = xnew;
    x = reshape(x, nely, nelx);
    r = r + delta_r;
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