%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 08 2023
% Implementation of transmission coefficients.
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R_BC] = robin_bc(Nx,dx)
%ROBIN_BC This function computes the robin boundary condition at BM
%    BM is assumed to be located at the right boundary!
%    Inputs: 
%        - Nx: number of cell centers along one row. (int)
%        - dx: grid size (float)
%    Outputs: 
%        - R_BC: the matrix with 1 or -1 to connect the two interfaces (sparse)
% 

R_BC = sparse(Nx+1,Nx+1); % sparse matrix of the same size as diffusion matrix
R_BC((Nx),(Nx)) = - 1 ; % the cell
R_BC((Nx),Nx+1) = 1;

R_BC(Nx+1,Nx+1) = -1;
R_BC(Nx+1,Nx) = 1;



end

