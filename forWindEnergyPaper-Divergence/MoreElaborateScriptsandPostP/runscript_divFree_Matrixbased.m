clear
clc
% close all

% WARNING: Currently the matrix inversion is done by backslash "\"
% See line 126
% RAM usage is poorly scaled with number of grid points
% Multigrid method recommended for large grids

%% Example usage : input is the divergence of the original wind field

%% See line 42 in "runscript_divFree_FFTbased.m" file

% div_v_og=div_v_og-mean(div_v_og(:)); %source term zero mean

% [x,y,z,u] = solvePoissonEquation3D_direct(Nx,Ny,Nz,Lx,Ly,Lz,div_v_og);
% dx=3.6328;
% dy=10;
% dz=10;

% [gradPhiX, gradPhiY, gradPhiZ] = computeGradient(u, dx, dy, dz); 

%% See line 32 in "runscript_divFree_FFTbased.m" file

% velocity_new(:,:,:,1) = velocity(:,:,:,1) - 1 .* gradPhiX;
% velocity_new(:,:,:,2) = velocity(:,:,:,2) - 1 .* gradPhiY;
% velocity_new(:,:,:,3) = velocity(:,:,:,3) - 1 .* gradPhiZ;

%%
function [x,y,z,u] = solvePoissonEquation3D_direct(Nx,Ny,Nz,Lx,Ly,Lz,f_source)
% solvePoissonEquation3D_direct solves the Poisson equation in 3D
%   u_xx + u_yy + u_zz = f_source
% on a domain with dimensions Lx, Ly, Lz (assumed to be centered at zero)
% with Neumann boundary conditions. To remove the singularity,
% the solution is fixed to u = 0 at one grid point.
%
% INPUTS:
%   Nx, Ny, Nz : number of grid cells in x, y, and z (must be integers > 1)
%   Lx, Ly, Lz : physical dimensions of the domain (positive scalars)
%   f_source   : 3D array of source term values of size [Nx,Ny,Nz]
%
% OUTPUTS:
%   x, y, z : 3D grids (using ndgrid) representing the cell centers
%   u       : 3D solution array (Nx x Ny x Nz)
%
% Example:
%   [x,y,z] = ndgrid( linspace(-1,1,Nx), linspace(-1,1,Ny), linspace(-1,1,Nz) );
%   f_source = -3*pi^2 * cos(pi*x).*cos(pi*y).*cos(pi*z);
%   u = solvePoissonEquation3D_direct(Nx,Ny,Nz,2,2,2,f_source);

arguments
    Nx (1,1) {mustBeInteger, mustBeGreaterThan(Nx,1)}
    Ny (1,1) {mustBeInteger, mustBeGreaterThan(Ny,1)}
    Nz (1,1) {mustBeInteger, mustBeGreaterThan(Nz,1)}
    Lx (1,1) {mustBeNumeric, mustBePositive}
    Ly (1,1) {mustBeNumeric, mustBePositive}
    Lz (1,1) {mustBeNumeric, mustBePositive}
    f_source (:,:,:) {mustBeNumeric}
end

% Optional: check that the size of f_source matches the grid dimensions
assert(isequal(size(f_source),[Nx,Ny,Nz]), 'f_source must be of size [Nx,Ny,Nz]');

% Grid spacing (assumes domain from -L/2 to L/2)
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;

% % Create the cell-centered grid (staggered grid as in the 2D case)
xx = -Lx/2 + dx/2 : dx : Lx/2 - dx/2;
yy = -Ly/2 + dy/2 : dy : Ly/2 - dy/2;
zz = -Lz/2 + dz/2 : dz : Lz/2 - dz/2;
%keyboard
[x,y,z] = ndgrid(xx,yy,zz);

% Assume the center of the domain is at x_0 (for example, x_0 = 1)
%x0 = 1;

% Create the grid centered at x_0 instead of 0
%xx = x0 - Lx/2 + dx/2 : dx : x0 + Lx/2 - dx/2;
%yy = -Ly/2 + dy/2 : dy : Ly/2 - dy/2;  % y-grid still centered at 0
%zz = -Lz/2 + dz/2 : dz : Lz/2 - dz/2;  % z-grid still centered at 0
%[x, y, z] = ndgrid(xx, yy, zz);


%--------------------------------------------------------------------------
% Construct 1D finite-difference matrices for the second derivative 
% with Neumann BC in each direction. For an interior point the stencil is:
%    (u(i-1) - 2u(i) + u(i+1))/h^2.
% At the boundaries the stencil is modified by setting the diagonal to -1.
%--------------------------------------------------------------------------

% x-direction operator
tmp = -2*ones(Nx,1);
tmp([1 end]) = -1;   % adjust for Neumann BC at the boundaries
Ax = diag(tmp) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1);

% y-direction operator
tmp = -2*ones(Ny,1);
tmp([1 end]) = -1;
Ay = diag(tmp) + diag(ones(Ny-1,1),1) + diag(ones(Ny-1,1),-1);

% z-direction operator
tmp = -2*ones(Nz,1);
tmp([1 end]) = -1;
Az = diag(tmp) + diag(ones(Nz-1,1),1) + diag(ones(Nz-1,1),-1);

%--------------------------------------------------------------------------
% Construct the 3D Laplacian using Kronecker products.
% The overall operator A is:
%   A = I_z ⊗ I_y ⊗ Ax/dx^2  +  I_z ⊗ Ay ⊗ I_x/dy^2  +  Az ⊗ I_y ⊗ I_x/dz^2.
%--------------------------------------------------------------------------
Ix = speye(Nx);
Iy = speye(Ny);
Iz = speye(Nz);

A = kron(Iz, kron(Iy, Ax))/(dx^2) + ...
    kron(Iz, kron(Ay, Ix))/(dy^2) + ...
    kron(Az, kron(Iy, Ix))/(dz^2);

%--------------------------------------------------------------------------
% Since the Laplacian with pure Neumann BC is singular, fix u at one point.
% Here we set u(1) = 0.
%--------------------------------------------------------------------------
f = f_source(:);  % convert source term to a vector
A(1,:) = 0;
A(1,1) = 1;
f(1) = 0;

% Solve the linear system
u = A \ f;

% Reshape the solution back to a 3D array
u = reshape(u, [Nx,Ny,Nz]);

end

function divU = computeDivergence(U, V, W, dx, dy, dz)
    [Nx, Ny, Nz] = size(U);
    divU = zeros(Nx, Ny, Nz);
    % Enforce Neumann BC: Copy boundary values from the interior
    U(1,:,:) = U(2,:,:);   U(end,:,:) = U(end-1,:,:);
    V(:,1,:) = V(:,2,:);   V(:,end,:) = V(:,end-1,:);
    W(:,:,1) = W(:,:,2);   W(:,:,end) = W(:,:,end-1,:);

    % Compute divergence in the interior using central differences
    divU(2:Nx-1, :, :) = (U(3:Nx, :, :) - U(1:Nx-2, :, :)) / (2*dx);
    divU(:, 2:Ny-1, :) = divU(:, 2:Ny-1, :) + (V(:, 3:Ny, :) - V(:, 1:Ny-2, :)) / (2*dy);
    divU(:, :, 2:Nz-1) = divU(:, :, 2:Nz-1) + (W(:, :, 3:Nz) - W(:, :, 1:Nz-2)) / (2*dz);

    % Enforce Neumann BCs using forward/backward differences
    divU(1, :, :) = (U(2, :, :) - U(1, :, :)) / dx;  % Left boundary
    divU(Nx, :, :) = (U(Nx, :, :) - U(Nx-1, :, :)) / dx;  % Right boundary

    divU(:, 1, :) = divU(:, 1, :) + (V(:, 2, :) - V(:, 1, :)) / dy;  % Bottom boundary
    divU(:, Ny, :) = divU(:, Ny, :) + (V(:, Ny, :) - V(:, Ny-1, :)) / dy;  % Top boundary

    divU(:, :, 1) = divU(:, :, 1) + (W(:, :, 2) - W(:, :, 1)) / dz;  % Front boundary
    divU(:, :, Nz) = divU(:, :, Nz) + (W(:, :, Nz) - W(:, :, Nz-1)) / dz;  % Back boundary
end

function [gradPhiX, gradPhiY, gradPhiZ] = computeGradient(Phi, dx, dy, dz)
    % Get the size of the scalar field
    [Nx, Ny, Nz] = size(Phi);
    
    % Initialize gradient components
    gradPhiX = zeros(Nx, Ny, Nz);
    gradPhiY = zeros(Nx, Ny, Nz);
    gradPhiZ = zeros(Nx, Ny, Nz);
    
    % Enforce Neumann BCs by copying boundary values
    Phi(1,:,:) = Phi(2,:,:);   Phi(end,:,:) = Phi(end-1,:,:);  % x boundaries
    Phi(:,1,:) = Phi(:,2,:);   Phi(:,end,:) = Phi(:,end-1,:);  % y boundaries
    Phi(:,:,1) = Phi(:,:,2);   Phi(:,:,end) = Phi(:,:,end-1);  % z boundaries

    % Compute the gradient using central differences for interior points
    gradPhiX(2:Nx-1, :, :) = (Phi(3:Nx, :, :) - Phi(1:Nx-2, :, :)) / (2*dx);
    gradPhiY(:, 2:Ny-1, :) = (Phi(:, 3:Ny, :) - Phi(:, 1:Ny-2, :)) / (2*dy);
    gradPhiZ(:, :, 2:Nz-1) = (Phi(:, :, 3:Nz) - Phi(:, :, 1:Nz-2)) / (2*dz);
    
    % Apply one-sided differences at the boundaries
    gradPhiX(1, :, :) = (Phi(2, :, :) - Phi(1, :, :)) / dx;
    gradPhiX(Nx, :, :) = (Phi(Nx, :, :) - Phi(Nx-1, :, :)) / dx;
    
    gradPhiY(:, 1, :) = (Phi(:, 2, :) - Phi(:, 1, :)) / dy;
    gradPhiY(:, Ny, :) = (Phi(:, Ny, :) - Phi(:, Ny-1, :)) / dy;
    
    gradPhiZ(:, :, 1) = (Phi(:, :, 2) - Phi(:, :, 1)) / dz;
    gradPhiZ(:, :, Nz) = (Phi(:, :, Nz) - Phi(:, :, Nz-1)) / dz;
end
