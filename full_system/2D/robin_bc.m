%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 08 2023
% 
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_BC] = robin_bc(Nx,dx)
%ROBIN_BC This function computes the robin boundary condition at BM
%    Inputs: 
%        - Nx: number of cell centers along one row. (int)
%        - dx: grid size (float)
%    Outputs: 
%        - R_BC: the matrix with 1 or -1 to connect the two interfaces (sparse)
% 

R_BC = sparse(Nx^2+Nx,Nx^2+Nx); % sparese matrix of the same size as ...
% the diffusion matrix plus the line of size Nx for the BM

for jj = 1:Nx
    % Top row: this side corresponds to the BM. 
    R_BC((Nx^2-(jj-1)),(Nx^2-(jj-1))) = - 1 ; % the cell
    R_BC((Nx^2-(jj-1)),Nx^2+Nx-(jj-1)) = 1;

    R_BC(Nx^2+Nx-(jj-1),Nx^2+Nx-(jj-1)) = -1;
    R_BC(Nx^2+Nx-(jj-1),Nx^2-(jj-1)) = 1;

end

end

