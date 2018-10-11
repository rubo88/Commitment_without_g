function [erg]=ErgDistUpwind(Pgrid,Pdot,S2) 
% This function calculates the ergodic distribution using upwind
% approximation.
% It requires the grid of the state space (1dim), the drift (Pdot) and the
% variace(S2)

%% Read input size
	N = numel(Pgrid)-1;    
    dPs=1/N;
% Define Delta t acordint to CTL    
    dt=(dPs).^2/max(S2);
% Compute discrete probabilities (upwind diff)
	pi_u=dt.*(S2./2+max(dPs.*Pdot,0))./((dPs).^2);
	pi_d=dt.*(S2./2+max(-dPs.*Pdot,0))./((dPs).^2);
	pi_m=1-pi_u-pi_d;
% Create transition matrix in Pgrid space     
    PI=zeros(N+1,N+1);
    PI(logical(eye(N+1)))=pi_m;
    for j=1:N+1
        for k=1:N+1
            if j==k+1
                PI(j,k)=pi_d(j);
            elseif k==j+1
                PI(j,k)=pi_u(j);
            end
        end
    end
% Solve for ergodic distribution     
    PI3=[PI'-eye(N+1) ; ones(1,N+1)];
    b=[zeros(1,N+1) 1]';
    erg=PI3\b./dPs;
end