%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 08 2023
% Construction of the diffusion matrix
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L] = diffusion_matrix(Nx,dx, BC_type)
    % This function computes the discrete laplacian in matrix form
    % It is assumed that the size of the cntrol volumes are already divided
    % hence /dx^2

    L = sparse(Nx,Nx)
    for ii=2:Nx-1 % interior points
        L(ii,ii) = -2/dx^2;
        L(ii,ii-1) = 1/dx^2;
        L(ii,ii+1) = 1/dx^2;
    end
    
    if BC_type == 1 % periodic
        L(1,1) = -2/dx^2;
        L(1,2) = 1/dx^2;
        L(1,Nx) = 1/dx^2;
        L(Nx,Nx) = -2/dx^2;
        L(Nx,Nx-1) = 1/dx^2;
        L(Nx,1) = 1/dx^2;
    elseif BC_type == 2 % Neumann zero
        L(1,1) = -1/dx^2;
        L(1,2) = 1/dx^2;
        L(Nx,Nx) = -1/dx^2;
        L(Nx,Nx-1) = 1/dx^2;
    elseif BC_type == 3 % dirichlet
        % nothing needed this is taken care in the "mass" matrix
        L(1,1) = 0;
        L(Nx,Nx) = 0;
    else
        exit("Wrong boundary condition asked")
    end
end
