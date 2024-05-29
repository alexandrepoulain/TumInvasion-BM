function IC = initial_conditions(Nx, test)
%INITIAL_CONDITIONS Tests the conditions that are derived in the article 

M = 50.e3; % BM density

cm = 0.0; % MT1-MMP monomer concentration
if test == 1
    cd = 0;
elseif  test == 2 | test == 3
    cd = 1e-5; % in nmol/mm^3
else
    error('Wrong test');
end
ct = 9.6;  
cp = 18.0;
cta = 17;
ca = 0.44; % ALT: 0.44; % in nmol/mm^3

% IC: ct cp in conjunctive tissue
IC = [ct*ones(Nx,1); cp*ones(Nx,1)];
% IC (BM): M cm cd cp  ct  ca  c1 c2 c3 ctp cta
IC = [IC; M; cm; cd; cp; ct; ca; 0; 0; 0; 0; cta];

end