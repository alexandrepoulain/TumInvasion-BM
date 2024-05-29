function IC = initial_conditions(Nx, test,dim)
% This function computes the initial conditions dependding on the test and
% the dimension of Omega
if dim == 1
    M = 50.0e3; % BM density nM
    cm = 0.0; % MT1-MMP monomer concentration
    if test == 1
        cd = 1e-4;
    elseif  test == 2 | test == 3
        cd = 1e-4; % in nmol/mm^3
    else
        error('Wrong test');
    end
    ct = 9.6; % TIMP-2
    cp = 18.0; % ProMMP-2
    cta = 17; % TIMP-2/active MMP-2
    ctp = 6.72; % TIMP-2/proMMP-2
    ca = 0.44; % active MMP-2: ALT: 0.44; % in nmol/mm^3 
    
    % IC: ct cp in conjunctive tissue
    IC = [ct*ones(Nx,1); cp*ones(Nx,1)];
    % IC (BM): M cm cd cp  ct  ca  c1 c2 c3 ctp cta
    IC = [IC; M; cm; cd; cp; ct; ca; 0; 0; 0; ctp; cta];

elseif dim == 2
    ct_BM = 9.6*ones(Nx,1); %TIMP-2
    cp_BM = 18*ones(Nx,1); % proMMP-2
    ca_BM = 0.44*ones(Nx,1); % active MMP-2
    cd_BM = 0.14*ones(Nx,1); % dimeric MT1-MMP
    cm_BM = 0.14*ones(Nx,1); % monomeric MT1-MMP
    c1_BM = zeros(Nx,1); % MT1-MMP/TIMP-2
    c2_BM = zeros(Nx,1);% MT1-MMP/TIMP-2/TIMP-2
    c3_BM = zeros(Nx,1);% MT1-MMP/TIMP-2/proMMP-2
    ctp_BM = 6.72*zeros(Nx,1); % TIMP-2/proMMP-2
    cta_BM = 17*ones(Nx,1); % TIMP-2/active MMP-2
    M = (50.e3)*ones(Nx,1); % BM density
    
    ct_conj = max(ct_BM) * ones(Nx^2,1); % in conjunctive tissue
    cp_conj = max(cp_BM) * ones(Nx^2,1);
    
    IC = [ct_conj; cp_conj; M; cm_BM; cd_BM; cp_BM; ct_BM; ca_BM; c1_BM; c2_BM; c3_BM; ctp_BM; cta_BM];
else
    error("Wrong dimension")
end

end
