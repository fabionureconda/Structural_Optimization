% PROBLEM SETUP
nelx = 24;          % Number of elements in x-direction
nely = 60;          % Number of elements in y-direction
volfrac = 1;        % Maximum volume fraction
penal = 3;          % Penalization power
rmin = 2;           % Filter radius in terms of elements
unitF = 10e3;

% MATERIAL PROPERTIES
E0 = 210e9;         % Young's modulus for the solid material
Emin = 1e-9*E0;     % Young's modulus for the void material
nu = 0.3;           % Poisson's ratio
t = 0.01;           % Thickness

% ELEMENT LENGTHS
l = 800/nely;       % Element length

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
fixeddofs = 1:2*(nely+1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs, fixeddofs);

% FORMULATION OF THE OPTIMIZATION PROBLEM
sig_max = 235;      % Maximum allowable stress
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
x = ones(nely, nelx);  % Initialization of design variables
xmin = 0;              % Minimum value design variables
xmax = 1;              % Maximum value design variables

% OPTIMIZATION LOOP
iter = 1;              % Iteration counter
change = 1;            % (Fictitious) initial design change
tol = 0.01;            % Convergence threshold
mmaparams = [];        % Internal MMA parameters
while change > tol
    % APPLY FILTER TO OBTAIN PHYSICAL DENSITIES
    xPhys = conv2(x, h, 'same')./Hs;
    xPhys = xPhys(:);

    % FINITE ELEMENT ANALYSIS
    sK = reshape(KE(:)*(Emin + xPhys(:)'.^penal * (E0 - Emin)), 64*nelx*nely, 1);
    K = sparse(iK, jK, sK); K = (K + K') / 2;
    U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

    % OBJECTIVE FUNCTION AND CONSTRAINT
    v = sum(xPhys(:)) / nelx / nely;
    f0 = v / volfrac - 1;

    % STRESS CALCULATION % no p
    term21v = zeros(1, nelx*nely);
    Uv = zeros(nelx * nely, 8);
    term22v = zeros(nelx*nely, 8);
    term11 = 0;
    term22 = 0;
    for el = 1:nelx*nely
        Ue = U(edofMat(el, :));
        sig_xxe = D0(1 ,:) * Be * Ue;
        sig_yye = D0(2, :) * Be * Ue;
        sig_xye = D0(3, :) * Be * Ue;
        sig_vMe = sqrt(sig_xxe^2 + sig_yye^2 - sig_xxe* sig_yye + 3 * sig_xye^2);
        Uv(el, :) = Ue'; % Store displacement vector for element

    % SENSITIVITIES
        % Sums: First term 11, Second term 22
            for ell = 1:nelx*nely
                Ue = U(edofMat(ell, :));
                sig_xxel = D0(1 ,:) * Be * Ue;
                sig_yyel = D0(2, :) * Be * Ue;
                sig_xyel = D0(3, :) * Be * Ue;
                sig_vMel = sqrt(sig_xxel^2 + sig_yye^2 - sig_xxel* sig_yyel + 3 * sig_xyel^2);
                
                term1 = ((sig_vMel / ((xPhys(ell)^(q-penal))* sig_max))^(r));
                term11 = term1+term11;
                f = term1;
                
                Ce = sparse(1:8, edofMat(ell, :), ones(1,8), 8, 2*(nely+1)*(nelx+1));

                term2a = ((sig_vMel / ((xPhys(ell)^(q-penal))* sig_max)))^(r-1);
                term2b = 1/(xPhys(ell)^(q-penal)*sig_max);
                term2c = 1/(2*sig_vMel);
                term2d = (((2*sig_xxel-sig_yyel)*D0(1,:)+(2*sig_xyel-sig_xxel)*D0(2,:)+(6*sig_xyel*D0(3,:)))*Be*Ce(:, ell))';
                term22 = term2a*term2b*term2c*term2d + term22;
            end
            term11 = term11.^(1/r-1);
            term22v(el, : ) = term22;

    % Term 21
    term21 = ((sig_vMe / ((xPhys(el)^(q-penal))* sig_max))^(r-1)) * (penal - q) * (sig_vMe / (xPhys(el)^(q-penal+1) * sig_max));
    term21v(el) = term21;
    end
    
    dpdxPhys = term11 * (term21v - term22v*(Uv'));   % how to link the U of the 8 nodes of every elements with the elements.

    Hs_new = reshape(Hs, 1, 1440);
    dvdxPhys = repmat(1/nelx/nely, nely, nelx); % Sensitivity wrt physical densities
    dpdx = conv2(dpdxPhys ./ Hs_new, h, 'same');    % Sensitivity wrt design variables
    dvdx = conv2(dvdxPhys ./ Hs, h, 'same');    % Sensitivity wrt design variables
    
    df0dx = dvdx;
    dfdx =  dpdx;

    % MMA UPDATE
    xdim = size(x)
    fdim = size(f)
    dfdxdim = size(dfdx)
    df0dxdim = size(df0dx)
    [xnew,~,~,~,mmaparams,~,change,history] = mma(x(:), xmin, xmax, f0, f, df0dx, dfdx, mmaparams);

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
    r = r + delta_r;
  end

end