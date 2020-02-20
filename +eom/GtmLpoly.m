classdef GtmLpoly < eom.GtmLong & aerootools.pkg.PolyApproximations
% Polynomial longitudinal EOM for the Generic Transport Model.
%
% Trigonometric and other non-polynomial functions are replaced by Taylor
% approximations (see aerotools.pkg.POLYAPPROXIMATIONS).
%
%% About
%
% * Author:     Torbjoern Cunis
% * Email:      <mailto:tcunis@umich.edu>
% * Created:    2020-02-11
% * Changed:    2020-02-11
%
% This file is part of GTMpw -- Piecewise polynomial model of the GTM
% published under the GNU General Public License v3.
%
%% Variables, constants, and their units
%
% * |alpha|    :  angle of attack,                              rad
% * |gamma|    :  flight-path angle,                            rad
% * |gammadot| :  change in flight-path angle,                  rad/s
% * |eta|      :  elevator deflection,                          rad
% * |rho|      :  air density,                                  kg/m^3
% * |b|        :  reference aerodynamic span,                   m
% * |c|        :  reference (mean) aerodynamic coord,           m
% * |Cdrag|    :  aerodynamic drag coefficient,                 -
% * |Clift|    :  aerodynamic lift coefficient,                 -
% * |Cm|       :  aerodynamic coefficient moment body y-axis,   -
% * |Cx|       :  aerodynamic coefficient force body x-axis,    -
% * |Cz|       :  aerodynamic coefficient force body z-axis,    -
% * |g|        :  gravitational constant,                       m/s^2
% * |Iy|       :  inertia y-axis,                               kg-m^2
% * |m|        :  aircraft mass,                                kg
% * |M|        :  pitch moment, body y-axis,                    N-m
% * |q|        :  pitch rate, body y-axis,                      rad/s
% * |qhat|     :  normalized pitch rate, body y-axis            rad
% * |S|        :  reference wing aera,                          m^2
% * |T|        :  thrust,                                       N
% * |V|        :  airspeed,                                     m/s
% * |Vdot|     :  change in airspeed,                           m/s^2
%%
%
% See also aerotools.EOM3, eom.GTMLONG
%
%%

properties (Access=protected)
    %AC     -- inherited from GtmLong
    %AM     -- inherited from GtmLong
    
    X0;
end

methods 
    function obj = GtmLpoly(x0, varargin)
        % Instance of longitudinal equations of motion
        % with reference trim condition |x0|
        % with aircraft model |AM|, aerodynamic coefficients |AC|, and
        % optional gravitational constant |g|.
        obj@eom.GtmLong(varargin{:});

        obj.X0 = eom.GtmLong.X(x0);
    end
end

methods (Access=protected)
    function iv = Vinv(obj, X, varargin)
        % inverse air speed
        iv = obj.inv(V(X),V(obj.X0));
    end
    
    function CD = Cdrag(obj, X, varargin)
        % polynomial aerodynamic drag
        if ~isa(obj.AM, 'aerootools.coeffs.AeroXZModel')
            CD = Cdrag@eom.GtmLong(obj,X,varargin{:});
            return
        end
        
        % else:
        CD = - obj.Cz(X, varargin{:}).*obj.sin(alpha(X)) ...
             - obj.Cx(X, varargin{:}).*obj.cos(alpha(X));
    end
    
    function CL = Clift(obj, X, varargin)
        % polynomial aerodynmic lift
        if ~isa(obj.AM, 'aerootools.coeffs.AeroXZModel')
            CL = Clift@eom.GtmLong(obj,X,varargin{:});
            return
        end
        
        % else:
        CL = - obj.Cz(X, varargin{:}).*obj.cos(alpha(X)) ...
             + obj.Cx(X, varargin{:}).*obj.sin(alpha(X));
    end
end

end
