%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Tue Nov 08 2023
% 
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Lx =  diffusion_matrix(Nx,dx)
%
%    This function computes the discrete gradient matrices: one in each direction
%    Inputs: 
%        - Nx: number of cell centers along one row. (int)
%        - dx: grid size (float)
%    Outputs: 
%        - Lx: the laplacian matrix (sparse)
% 
    Lx = sparse(Nx^2,Nx^2);
    % For interior points
    for ii=2:Nx-1
        for jj=2:Nx-1
            Lx((ii-1)*Nx+jj,(ii-1)*Nx+jj) = -2/dx^2 - 2/dx^2;
            Lx((ii-1)*Nx+jj,(ii-1)*Nx+jj-1) = 1/dx^2;
            Lx((ii-1)*Nx+jj,(ii-1)*Nx+jj+1) = 1/dx^2;

            Lx((ii-1)*Nx+jj,(ii)*(Nx)+jj ) = 1/dx^2;
            Lx((ii-1)*Nx+jj,(ii-2)*(Nx)+jj ) = 1/dx^2;
        end
    end
    for jj=1:Nx
        % bottom row: zero flux
        Lx(jj,jj) = -2/dx^2 - 1/dx^2; % the cell in question

        if jj == 1 % if bottom left corner
            Lx(jj,jj) = -1/dx^2 - 1/dx^2; % nothing there
        else
            Lx(jj,jj-1) = 1/dx^2; % left cell 
        end
        if jj == Nx % if rightmost cell of bottom line
            Lx(jj,jj) = -1/dx^2 - 1/dx^2; % nothing there
        else
            Lx(jj,jj+1) = 1/dx^2; % right cell
        end

        Lx(jj,jj+Nx)= 1/dx^2; % cell above
        % Lx[jj,jj-Nx]= 1/dx**2 % cell below ie nothing

        % Top row: zero flux (another matrix will take care of implementing  ...
        % the robin boundary conditioon on this side) 
        Lx((Nx^2-(jj-1)),(Nx^2-(jj-1))) = -2/dx^2 - 1/dx^2; % the cell
        if jj == Nx % top left boundary
            Lx((Nx^2-(jj-1)), Nx^2-(jj-1)) = -1/dx^2 - 1/dx^2; % nothing here 1/dx^2
        else
            Lx((Nx^2-(jj-1)),(Nx^2-1-(jj-1))) = 1/dx^2; % the cell on the left
        end
        if jj == 1 % top right corner
            Lx((Nx^2-(jj-1)),(Nx^2-(jj-1))) = -1/dx^2 - 1/dx^2;
        else
            Lx((Nx^2-(jj-1)),(Nx^2-(jj-1)+1)) = 1/dx^2; % the cell on the right
        end
        % Lx[(Nx**2-1-jj),Nx-1-jj] = 1/dx**2 nothing here
        Lx((Nx^2-(jj-1)),Nx^2-Nx-(jj-1)) = 1/dx^2; % The cell below

        % left boundary (not top or bottom points)
        if jj ~= 1 && jj ~= Nx
            Lx((jj-1)*Nx+1,(jj-1)*Nx+1) =  -1/dx^2 - 2/dx^2; % the cell
            Lx((jj-1)*Nx+1,(jj-1)*Nx+2) = 1/dx^2; % the cell on the right
            % Lx(jj*Nx,jj*Nx+(Nx)) = 1/dx**2 nothing here
            Lx((jj-1)*Nx+1,(jj-1)*Nx+Nx+1) = 1/dx^2; % cell above
            Lx((jj-1)*Nx+1,(jj-1)*Nx-Nx+1) = 1/dx^2; % cell below
        end

        % right boundary  (not top or bottom points)
        if jj ~= 1 && jj~= Nx
            Lx((jj-1)*Nx+(Nx),(jj-1)*Nx+(Nx)) =  -1/dx^2 - 2/dx^2; % the cell
            Lx((jj-1)*Nx+(Nx),(jj-1)*Nx+(Nx)-1) = 1/dx^2; % the cell on the left
            %Lx[jj*Nx+(Nx-1),jj*Nx ] = 1/dx**2; Nothing here
            Lx((jj-1)*Nx+(Nx),(jj-1)*Nx+(Nx)+Nx) = 1/dx^2; % the cell above
            Lx((jj-1)*Nx+(Nx),(jj-1)*Nx+(Nx)-Nx) = 1/dx^2; % the cell below
        end
    end   
    
end