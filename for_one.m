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
F = F1+F2+F3+F4+F5;                                                          % Load distribution
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
    f = v / volfrac - 1;

    % STRESS CALCULATION  % no p !!!!!!!!!!!!!!!!!!!!!
    term21v = zeros(1, nelx*nely); 
    Uv = zeros(nelx * nely, 8);
    term11 = 0;
    term22v = zeros(nelx*nely, 8);
    term22 = 0;
    for el = 1:nelx*nely
        Ue = U(edofMat(el, :));
        sig_xxe = D0(1 ,:) * Be * Ue;
        sig_yye = D0(2, :) * Be * Ue;
        sig_xye = D0(3, :) * Be * Ue;
        sig_vMe = sqrt(sig_xxe^2 + sig_yye^2 - sig_xxe* sig_yye + 3 * sig_xye^2);
        Uv(el, :) = Ue';                                                               % Store displacement vector for element

        % SENSITIVITIES
        % Sums: First term 11, Second term 22 and the term21 that is not a
        % sum (here we assume that the term22 is just K*LAMDA*U)
                
        term1 = ((sig_vMe / ((xPhys(el)^(q-penal))* sig_max))^(r));
        term11 = sum(term1);
        term11 = term11.^(1/r-1);
                
        Ce = sparse(1:8, edofMat(el, :), ones(1,8), 8, 2*(nely+1)*(nelx+1));            %Location Matrix

        term2a = ((sig_vMe / ((xPhys(el)^(q-penal))* sig_max)))^(r-1);
        term2b = 1/(xPhys(el)^(q-penal)*sig_max);
        term2c = 1/(2*sig_vMe);
        term2d = (((2*sig_xxe-sig_yye)*D0(1,:)+(2*sig_xye-sig_xxe)*D0(2,:)+(6*sig_xye*D0(3,:)))*Be*Ce(:, el))';
        term22 = sum(term2a*term2b*term2c*term2d);
            
        %term22v(el, :) = term22;

        % Term 21 (not a sum)
        term21 = ((sig_vMe / ((xPhys(el)^(q-penal))* sig_max))^(r-1)) * (penal - q) * (sig_vMe / (xPhys(el)^(q-penal+1) * sig_max));
        term21v(el) = term21;
        
    end
    
    dpdxPhys = term11 * (term21v - term22*(Uv));   % how to link the U of the 8 nodes of every elements with the elements.
    
    Hs_new = reshape(Hs, 1, 1440);
    dvdxPhys = repmat(1/nelx/nely, nely, nelx); % Sensitivity wrt physical densities
    dpdx = conv2(dpdxPhys ./ Hs_new, h, 'same');    % Sensitivity wrt design variables
    dvdx = conv2(dvdxPhys ./ Hs, h, 'same');    % Sensitivity wrt design variables
    
    df0dx = dpdx(:)' / 100;
    dfdx = dvdx(:)' / volfrac;
    
    % MMA UPDATE
    x = x(:)';
    xdim = size(x)
    fdim = size(f)
    dfdxdim = size(dfdx)
    df0dxdim = size(df0dx)
    [xnew,~,~,~,mmaparams,~,change,history] = mma(x, xmin, xmax, f, df0dx, dfdx, mmaparams);

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



% meanign of the sum, in the sense that how e is different from j in term
% of pedices? 

% how should we derive the term K*Lambda? what is the meaning of that term
% in the sensitivities with respect to K*Lamda?

% how are we able to put the constrain in the P value? do we need to create
% a if loop in order to check if the P -1 < 0? Or for now we don't need
% that term and it will only be use in the others points of the assignemnt?

% why are we facing a problem in the dimensions of dfdx if they are
% correct? why is matlab still exiting the error? there should be some
% other problems to solve in order to arrive to the solution

% we should increase the value of r at every iteration as said in the
% paper, how can we do that? is it correct to have the implementation we
% already create?

% is correct to change the shape of Hs are reshape it in a vector of 1440
% entries? or should we change the code in order to leave it like it was
% original?

% if we are dealing with all the nodes of the structure in the vector U
% (since it is build with edofMat it has rows = el and for every rows it
% has 8 entries of the nodes) and compare them to the elements? if we are
% building the for loop in the elements we get a total of 1440 elements but
% not 3550 values (ndof)