classdef AxCaliberSMT1 < handle
    properties (Constant = true, Access = protected)
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
        b;
        Delta;
        delta;
        q;
        D0;
    end
    
    properties (GetAccess = private, SetAccess = protected)
        
    end
    
    methods (Access = public)
        function this = AxCaliberSMT1(b, delta, Delta, D0)
%AXCALIBERSMT1 Axon size estimation using SMT model of IAS
% smt = AxCaliberSMT1(bval, delta, Delta, D0)
%       output:
%           - smt: object of a fitting class
%
%       input:
%           - bval: b-value [ms/um2]
%           - delta: gradient duration [ms]
%           - Delta: gradient seperation [ms]
%           - D0: intra-cellular intrinsic diffusivity [um2/ms]
%
%       usage:
%           smt = AxCaliberSMT1(bval, delta, Delta, D0);
%           model = 'Neuman';
%           [r, beta]  = smt.getAxonRadius(S, model);
%           Sfit = smt.AxonDiameterFWD([r, beta], model);
%
%  Authors: 
%  Hong-Hsi Lee (hlee84@mgh.harvard.edu)
%  Copyright (c) 2022 Massachusetts General Hospital
%
%  Adapted from the code of
%  Jelle Veraart (jelle.veraart@nyulangone.org)
%  Copyright (c) 2019 New York University
            
            this.b = b;
            this.Delta = Delta;
            this.delta = delta;
            % q is the gradient strength in unit 1/um/ms
            this.q = sqrt(b./delta.^2./(Delta-delta/3));
            this.D0 = D0;
        end
        
        function [r, beta]  = mcmc(this, y, model, N)
            r = zeros(N,1);
            beta = zeros(N,1);
            
            for i = 1:N
                [r(i), beta(i)] = this.getAxonRadius(y, model);
            end
            idx = kmeans([r, beta], 2);
            I = idx == mode(idx);
            r = mean(r(I));
            beta = mean(beta(I));
            
        end
        
        function [r, beta]  = getAxonRadius(this, y, model)
%GETAXONRADIUS Estimate the axon radius from multi-shell data.
% [r, beta]  = getAxonRadius(y, model)
%       output:
%           - r: effective MR radius
%           - beta: sqrt(pi/4/Da), Da = pi/4/beta^2
%
%       input:
%           - y:      powder-averaged diffusion-weighted signal,
%                     normalized to y(b=0) = 1. If the dot comparment
%                     cannot be ignored, it need to be subtracted from y prior to fitting.  
%           - model:  VanGelderen or Neuman (the latter only applied in long pulse regimes, but typically applies)
%
%  Authors: 
%  Maxina Sheft (msheft@mit.edu)
%  Hong-Hsi Lee (hlee84@mgh.harvard.edu)
%  Copyright (c) 2022 Massachusetts General Hospital
%
%  Adapted from the code of
%  Jelle Veraart (jelle.veraart@nyulangone.org)
%  Copyright (c) 2019 New York University

            if ~exist('model', 'var')
                model = 'Neuman';
            end

            options = optimset('lsqnonlin'); 
            options = optimset(options,'Jacobian','on','TolFun',1e-12,'TolX',1e-12,'MaxIter',10000,'Display','off');
            
            % Fitting parameters: [r, beta]
            start = [1.5+2*rand, sqrt(4*pi)*rand]; % initial values
            lb = [0, 0];            % lower bound
            ub = [5, sqrt(4*pi)];   % upper bound
            
            pars = lsqnonlin(@(x)this.residuals(x, y, model),...
                start,lb,ub,options);

            r    = pars(1);
            beta = pars(2);
        end
        
        function [E, J] = residuals(this, pars, y, model)
            [shat, dshat] = this.AxonDiameterFWD(pars, model); 
            E = shat(:) - y(:);
            J = dshat;
        end
        
        function [s, ds] = AxonDiameterFWD(this, pars, model)
            r    = pars(1);
            beta = pars(2);
            
            % Forward model
            switch model
                case 'Neuman'
                    [C, dC] = this.neuman(this.delta, this.q, r, this.D0);
                case 'VanGelderen'
                    [C, dC] = this.vg(this.delta, this.Delta, this.q, r, this.D0);
            end
            S = exp(-C);
            dS = S.*(-dC);
            
            erf_a = erf(sqrt(this.b*pi/4)./beta);
            s = beta./sqrt(this.b) .* S .* erf_a;
            ds_dbeta = 1./sqrt(this.b) .* S  .* erf_a + ...
                beta./sqrt(this.b) .* S .* exp(-this.b*pi/4./beta.^2).*sqrt(this.b)./(-beta.^2);
            ds_dr = beta./sqrt(this.b) .* dS .* erf_a;
            
            ds = [ds_dr, ds_dbeta];
        end
        
    end
    
    methods(Static)
        function [s, ds] = vg(delta, Delta, q, r, D0)
            td=r^2/D0;
            bardelta=delta/td; dbardelta = -2*delta*D0 / r^3;
            barDelta=Delta/td; dbarDelta = -2*Delta*D0 / r^3;

            N=15; 
            b = [1.8412    5.3314    8.5363   11.7060   14.8636   18.0155   21.1644 24.3113   27.4571   30.6019 ...
                 33.7462   36.8900   40.0334   43.1766   46.3196 49.4624   52.6050   55.7476   58.8900   62.0323];

            s = 0; ds = 0;
            for k=1:N
               s = s + (2/(b(k)^6*(b(k)^2-1)))*(-2 ...
                                + 2*b(k)^2*bardelta ...
                                + 2*exp(-b(k)^2*bardelta) ...
                                + 2*exp(-b(k)^2*barDelta) ...
                                - exp(-b(k)^2*(barDelta+bardelta))...
                                - exp(-b(k)^2*(barDelta-bardelta))); 
               ds = ds +  (2/(b(k)^6*(b(k)^2-1)))*( ...
                                2*b(k)^2*dbardelta ...
                                + 2*exp(-b(k)^2*bardelta) *(- b(k)^2*dbardelta) ...
                                + 2*exp(-b(k)^2*barDelta) *(- b(k)^2*dbarDelta) ...
                                - exp(-b(k)^2*(barDelta+bardelta)) *(- b(k)^2*(dbarDelta + dbardelta)) ... 
                                - exp(-b(k)^2*(barDelta-bardelta)) *(- b(k)^2*(dbarDelta - dbardelta)));           
            end
            ds = ds.*D0.*q.^2.*td^3 + 6*s.*q.^2.*r^5 / D0^2;
            s = s.*D0.*q.^2.*td^3;
        end
        
        function [s, ds] = neuman(delta, q, r, D0)
            s = (7/48)*q.^2.*delta*r^4/D0;
            ds = (7/12)*q.^2.*delta*r^3/D0;
        end
        
    end
end
