%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 15 2023
% Function used to create the matrix assciated to the semi discretization 
% in space.
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We use sparse matrices and construct a block matrix for the PDE-ODE coupled 
% system.
% To fix things, we assume the BM is lovcated at the top boundary. 

function Mat = construct_matrix(X,Y, dx,dy, params)
  % X,Y: the grid 
  % dx, dy: grid size in x and y-directions
  % params: the parameters of the model (see function parameters to see what it
  % contains).  
  
  width_BM = 

  D_t = 
  D_p =
  kappa_t = 
  kappa_p = 
  beta_t = 
  beta_p = 
  beta_m =  
  beta_d = 
  beta_a =
  beta_tp =
  beta_ta =

  S_t = 
  S_p = 

  kappa_t = 
  kappa_p = 

  alpha_t = 
  alpha_p = 

  rho_0 = 

  k0 = 
  km0 = 
  k1 =
  km1 = 
  k2 = 
  km2 = 
  k3 = 

  gamma = 
  K_M = 
  r_M =
  Mmax =   
  
  % Construct diffusion matrix 
  K = sparse(size(X)); 
  for i = 2:end-1
    K(i,i) = -2/dx^2-2/dy^2;
    K(i,i+Nx) = 1/dy;
    K(i,i-Nx) = 1/dy; 
    K(i,i+1) = 1/dx; 
    K(i,i-1) = 1/dx;  
  endfor
  % BCs 
  i = 1;
  K(i,i) = -2/dx^2-2/dy^2;
  K(i,i+Nx) = 1/dy;
  K(i,i-Nx) = 1/dy; 
  K(i,i+1) = 1/dx; 
  K(i,i-1) = 1/dx; 
 
  
end