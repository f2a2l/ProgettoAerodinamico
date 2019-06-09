classdef solverVHS
    %solverVHS Inviscid + viscid aerodynamic solver.
    %   Uses a Hess-Smith method for external flow; boundary layer is
    %   accounted for with an integral equation with Thwaites' closure
    %   (laminar flow) or Head's (turbulent).
    
    properties

        nArfls
        alpha
        
        inviscidCL
        inviscidCD
        inviscidCFs
        
        x
        y
        xmax
        ymax

        Cp
        maxdCp

        panels
        metaPan

        SOL

    end
    
    methods
        function obj = solverVHS(npoint, aname, alpha, varargin)
            %solverVHS Constructor: calls inviscid solver for the specified
            %geometry and saves data. Keep in mind that external viscous
            %solution depends only on the geometry, not on Re.

            if length(alpha) == 1 && isempty(varargin)
                [Cl, Cd, xmax, ymax, Cp, maxdCp, x, y, p, SOL, metaPan, nairfoils] = solverHS(npoint, aname, alpha);
            elseif length(alpha) > 1 && length(varargin) == 2
                dist = varargin{1};
                crel = varargin{2};
                [Cl, Cd, xmax, ymax, Cp, maxdCp, x, y, p, SOL, metaPan, nairfoils] = solverHS(npoint, aname, alpha, dist, crel);
            else
                error('wrong input; please check documentation.')
            end
                        
            obj.nArfls = nairfoils;
            obj.alpha = alpha;

            obj.inviscidCFs{2} = Cd;
            obj.inviscidCFs{1} = Cl;
            obj.inviscidCL = sum(Cl);
            obj.inviscidCD = sum(Cd);

            obj.x = x;
            obj.y = y;
            obj.xmax = xmax;
            obj.ymax = ymax;

            obj.Cp = Cp;
            obj.maxdCp = maxdCp;

            obj.panels = p;
            obj.metaPan = metaPan;

            obj.SOL = SOL;

        end


        
        function plotStreamlines(obj,varargin)
            %plotStreamlines Plots streamlines of inviscid solution.

            if isempty(varargin)
                showNodes = true;
            elseif length(varargin) == 1
                showNodes = varargin{1};
                if ~islogical(showNodes)
                    error('input must be true or false.')
                end
            else
                error('too many inputs; just pass true or false to plot nodes or not.')
            end

            % load input from class
            Cl = obj.inviscidCFs{1};
            Cd = obj.inviscidCFs{2};

            % call plot function
            PlotStream(obj.nArfls, obj.panels, obj.x, obj.y, obj.metaPan, obj.alpha, obj.SOL, Cl, Cd, showNodes);

        end




    end
end

