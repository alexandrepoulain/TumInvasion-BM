%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created on Thu Oct 23 2023
% Function that tests the stability of the steady states of the reduced system.
% @author: Alexandre Poulain, Chiara Villa
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = test_stability(kappa_p,kappa_t,k1,k2,k3,km1,km2,Km,gamma,rM...
    ,Mmax,ct,cp,ca,cd, S_t , S_p, rho0, alpha_t,alpha_p,width_BM)
%TEST_STABILITY Tests the conditions that are derived in the article 
    kappa_p = kappa_p/width_BM;
    kappa_t = kappa_t/width_BM;
    
    S_t = kappa_p*ct + rho0*alpha_t;
    S_p = kappa_p*cp + rho0*alpha_p;
    K_R = k2+ k1*k2*k3*cd/(km1*(km2+k3));
    

    % Computing the equilibrium values
    cp = sqrt( ( (S_t+S_p)/(2*kappa_p) + kappa_t/(2*K_R) )^2 - S_t*S_p/kappa_p^2 ) - (S_t-S_p)/(2*kappa_p) - kappa_t/(2*K_R);
    ca = k1*k2*cd/(km1*(km2+k3))*cp; % in nmol/mm^3
    ct = (S_t-S_p)/kappa_t + kappa_p/kappa_t*cp;  

    M_crit = Mmax/2;
    if (ca > (Km+M_crit)/(gamma * M_crit)*rM*(1-M_crit/Mmax))
        "Degradation"
    else 
        "No degradation"
    end
    T1 = kappa_p + kappa_t + K_R*(ct+cp) + k2*ct;
    T2 = kappa_p*kappa_t + K_R *(kappa_p*cp+kappa_t*ct) +k2*ct*(kappa_p +...
        kappa_t)+  K_R^2*cp*ct + K_R*k2*ct^2;
    T3 = k2*ct*kappa_t*kappa_p + K_R * k2*ct^2*kappa_t +K_R*k2*ct*cp*kappa_p;
    T4 = T1*T2/6 + ((T2/3 -T1^2/9)^3 + (T1^3/27 - T1*T2/6 +T3/2)^2 )^(1/2)...
        - T1^3/27 - T3/2;
    
    if (T4^(1/3)-T1/3-(T2/3-T1^2/9)/(T4^(1/3)) <0 && (T2/3-T1^2/9)/(2*T4^(1/3)) -T1/3 - T4^(1/3)/2 < 0)
        "stability"
    elseif ((T2/3-T1^2/9)/(2*T4^(1/3)) -T1/3 - T4^(1/3)/2 == 0) && (T4^(1/3) + (T2/3 - T1^2/9)/T4^(1/3) ~= 0)
        "oscillate"
    else
        "unstable system"
    end
end

