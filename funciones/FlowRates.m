    function [fup, fdn, flv] = FlowRates(a,s,dP)
        flv = s.^2 ./ dP^2;         % Flow rate of leaving current grid point
        fup = 0.5*(flv + a/dP);     % Flow rate of going one grid point up.
        fdn = flv - fup;            % Flow rate of going one grid point down.

        fdn(1  ) = 0;               % Make boundaries absorbing: No Brownian
        fup(1  ) = max( a(1  )/dP, 0);  % shocks here, only drift may take us 
        flv(1  ) = fup(1  );
        fup(end) = 0;               % into grid.
        fdn(end) = min(-a(end)/dP, 0);
        flv(end) = fdn(end);
    end   