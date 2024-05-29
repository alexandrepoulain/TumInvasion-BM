%% FUNCTION TO EXECUTE PDE+ODE system

function ret_vec = PDE_solve(y,test_case,system,par_case,L,Nx,x,dx,dim)

%%% Parameter values
[gamma,Km,rM,k0,km0,k1,k2,k3,km1,km2,beta_t,beta_ta,beta_m,...
    beta_d,beta_p,beta_a,beta_tp,Mmax,rho0,alpha_m,D_t,D_p,width_BM,...
    kappa_t,kappa_p,sph_t,sph_p,posSF,rt,rp,spread_secretum,alpha_t, alpha_p] = parameters(test_case,L,dim,Nx,x);

% Time parameters
Tf = 1;   % Final time (days)
Tp = 100; % Number of time points to store solution


dx = L/Nx; % the grid size

if dim == 1
    x = linspace(dx/2,L-dx/2,Nx);
    x_inter = linspace(0,L,Nx+1);
elseif dim == 2
    x = linspace(dx/2,L-dx/2,Nx);
    y = linspace(dx/2,L-dx/2,Nx);
    x_inter = linspace(0,L,Nx+1);
    y_inter = linspace(0,L,Nx+1);
    [X,Yx] = meshgrid(x,y);
else
    error("Wrong dimension")
end

% Find position in x array of senescent fibroblast
% xi_SF = find(x==interp1(x,x,pos*L,'nearest'));
xi_SF = posSF*L;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update parameters under SA (make sure entries match Input definition)
if par_case == 1
    if system == 2 % if reduced system cd is an input
        % Block A, for EE
        k1 = y(1);
        km1 = y(2);
        k2 = y(3);
        km2 = y(4);
        k3 = y(5);
        %%% If need to update cd in IC
        beta_t = y(6);
        beta_p = y(7);
        cd = y(8);
    %IC(end-8) = cd;
    elseif system == 1
        k0 = y(1);
        km0 = y(2);
        k1 = y(3);
        km1 = y(4);
        k2 = y(5);
        km2 = y(6);
        k3 = y(7);
        beta_t = y(8);
        beta_ta = y(9);
        beta_m = y(10);
        beta_d = y(11);
        beta_p = y(12);
        beta_a = y(13);
        beta_tp = y(14);

    end

elseif par_case == 2
    if system == 1
        % Block B, for EE
        gamma = y(1) ;
        Km = y(2);
        rM = y(3);
        rho0 = y(4);
        alpha_m = y(5);
        D_t = y(6);
        D_p = y(7);
        kappa_t = y(8);
        kappa_p = y(9);
        sph_t = y(10);
        sph_p = y(11);
        rt = y(12);
        rp = y(13);
        posSF = y(14);
  


    elseif system == 2
        % Block B, for EE
        
        D_t = y(3);
        D_p = y(4);
        kappa_t = y(5);
        kappa_p = y(6);
        sph_t = y(7);
        sph_p = y(8);
        rt = y(9);
        rp = y(10);
        gamma = y(11);
        rM = y(12);
        Km = y(13);
        posSF = y(14);
    end
% 
% elseif par_case == 3 % For FAST and VBSA
%     k2 = y(1);
%     alpha_m = y(2);
%     D_t = y(3);
%     D_p = y(4);
%     rt = y(5);
%     rp = y(6);
%     posSF = y(7);
    
else        
    error(['\n WARNING: specify which block of parameters to study ' ...
            '(A: biochemical reaction rates, B: all others)'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% the source terms \section{Checking numerically the stability results}

% the source terms 
if dim == 1
    
    % Computation of the cell average of SF secretome 
    sf_secretome_t = @(x) (rt./(spread_secretum * sqrt(2*pi))*exp(-(x-posSF).^2/(2*spread_secretum^2) ));
    sf_secretome_p = @(x) (rp./(spread_secretum * sqrt(2*pi))*exp(-(x-posSF).^2/(2*spread_secretum^2) ));
    S_t = zeros(Nx, 1);
    S_p = zeros(Nx, 1);
    for ii=1:Nx
        S_t(ii) = sph_t + 1./dx*integral(sf_secretome_t,x_inter(ii),x_inter(ii+1)); 
        S_p(ii) = sph_p + 1./dx*integral(sf_secretome_p,x_inter(ii),x_inter(ii+1)); 
    end
elseif dim == 2  % pos is a matrix, nb line is the sumber of SFs, size 2 is 2  
    S_t = sph_t*ones(Nx*Nx, 1);
    S_p = sph_p*ones(Nx*Nx, 1);
    for ii = 1:size(posSF,1)
        %xSF =  reshape(X,[Nx^2,1]);
        xSF = X;
        %ySF =  reshape(Yx,[Nx^2,1]);
        ySF = Yx;
        
        sf_secretome_t = @(x,y) rt./(spread_secretum * sqrt(2*pi))*exp(...
                -(x-posSF(ii,1)).^2/(2*spread_secretum^2) -(y-posSF(ii,2)).^2/(2*spread_secretum^2) );
        sf_secretome_p = @(x,y) rp./(spread_secretum * sqrt(2*pi))*exp(...
                -(x-posSF(ii,1)).^2/(2*spread_secretum^2) -(y-posSF(ii,2)).^2/(2*spread_secretum^2) );
        
        

        for jj = 1:Nx % iterate on lines
            for rr = 1:Nx
                xmin = x_inter(rr);
                xmax = x_inter(rr+1);
                ymin = y_inter(jj);
                ymax = y_inter(jj+1);
                S_t((jj-1)*Nx+rr) = S_t((jj-1)*Nx+rr) + 1./(dx^2)*integral2(sf_secretome_t,xmin,xmax,ymin,ymax); 
                S_p((jj-1)*Nx+rr) = S_p((jj-1)*Nx+rr) +1./(dx^2)*integral2(sf_secretome_p,xmin,xmax,ymin,ymax); 
            end
        end
    end
    % recast S_t and S_p in array form
%     S_t_arr = zeros(Nx*Nx,1);
%     S_p_arr = zeros(Nx*Nx,1);
%     
%     for ii = 1:size(S_t,1)
%         S_t_arr((ii-1)*Nx+1:(ii)*Nx) = S_t(ii,:);
%         S_p_arr((ii-1)*Nx+1:(ii)*Nx) = S_p(ii,:);
%     end
%     S_t = S_t_arr;
%     S_p = S_p_arr;
end



%%% Solve the PDE+ODE system
switch system
    case 2 % Solve reduced system
        addpath('../reduced_system/');
        % Initial concentrations 
        IC = initial_conditions(Nx, test_case);
        if par_case == 1
            IC(2*Nx+3) = cd; 
        end
        % simulation of reduced system
        [t_BM,ct_conj, cp_conj, cp_BM,ct_BM,ca_BM,M_BM,cm_BM,cd_BM,c1_BM,c2_BM,c3_BM,ctp_BM,cta_BM]...
            = reduced_system(IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
            rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
            beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_p, D_t, width_BM, ...
            0, 0);

    case 1
        addpath('../full_system/');
        % Initial concentrations 
        IC = initial_conditions(Nx, test_case,dim);
        % simulation of full system
        [t_BM,ct_conj, cp_conj, cp_BM,ct_BM,ca_BM,M_BM,cm_BM,cd_BM,c1_BM,c2_BM,c3_BM,ctp_BM,cta_BM]...
            = full_system(dim,IC, Nx, dx, Tf, Tp, gamma, Km, rM, Mmax, alpha_m, ...
            rho0, k0, km0, beta_m, k1, km1, beta_d, S_p, k2, km2, beta_p, S_t, ...
            beta_t, k3, beta_a, beta_tp, beta_ta, kappa_t, kappa_p, D_p, D_t, ...
            width_BM );

    otherwise
        fprintf('\n WARNING: specify which system must be solved (reduced, full)')
end

% %%% Check convergence to steady state 
% if abs(ca_BM(end)-ca_BM(end-1))>1e-6
%     fprintf('\n WARNING: Solution may not have reached steady state')
% end

%%% Return output
cas = ca_BM(end);
Mas = M_BM(end);



if M_BM(end) < 0 | isnan(M_BM(end))
    'problem'
    M_BM(end)
    quit
else
    'OK'
    M_BM(end)
    % Return full
    ret_vec = Tf/Tp*sum(abs(M_BM))/abs(IC(Nx*2+1));%[cp_BM(end),ct_BM(end),ca_BM(end),M_BM(end),cm_BM(end),...
    %cd_BM(end),c1_BM(end),c2_BM(end),c3_BM(end),ctp_BM(end),cta_BM(end)];
end

%%% IF full system, output cd
% cas = cd_BM(end);

end

