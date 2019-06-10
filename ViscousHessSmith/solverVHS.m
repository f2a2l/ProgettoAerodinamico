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
        xpanc
        ypanc
        xmax
        ymax

        Cp
        ue
        maxdCp

        panels
        metaPan

        SOL

        stagnIdx
        bl

    end
    
    methods



        function obj = solverVHS(npoint, aname, alpha, varargin)
            %solverVHS Constructor: calls inviscid solver for the specified
            %geometry and saves data. Keep in mind that external viscous
            %solution depends only on the geometry, not on Re.

            if length(alpha) == 1 && isempty(varargin)
                [Cl, Cd, xmax, ymax, Cp, v, maxdCp, x, y, p, SOL, metaPan, nairfoils] = solverHS(npoint, aname, alpha);
            elseif length(alpha) > 1 && length(varargin) == 2
                dist = varargin{1};
                crel = varargin{2};
                [Cl, Cd, xmax, ymax, Cp, v, maxdCp, x, y, p, SOL, metaPan, nairfoils] = solverHS(npoint, aname, alpha, dist, crel);
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
            obj.ue = v;
            obj.maxdCp = maxdCp;

            obj.panels = p;
            obj.metaPan = metaPan;

            obj.SOL = SOL;

            obj.stagnIdx = zeros(1,nairfoils);

            tic

            % get panel info FIXME: this won't be necessary if you change the way panels are handled
            for k = nairfoils:-1:1
                xe = x{k};
                ye = y{k};
                xc = zeros(size(xe));
                xc(end) = [];
                yc = xc;
                for ii = 1:length(xc)
                    xc(ii) = (xe(ii) + xe(ii+1))/2;
                    yc(ii) = (ye(ii) + ye(ii+1))/2;
                end
                xpanc{k} = xc;
                ypanc{k} = yc;
            end
            obj.xpanc = xpanc;
            obj.ypanc = ypanc;

            % preprocess BL info
            bldata = cell(nairfoils, 2);
            for k = nairfoils:-1:1
                cptemp = Cp{k};
                if k == 1
                    cap = 15;
                else
                    cap = 5;
                end
                for ii = 1:cap % temporarily set trailing edge Cp values to 0
                    cptemp(ii) = 0;
                    cptemp(end+1-ii) = 0;
                end
                [~, idx] = max(cptemp);
                obj.stagnIdx(k) = idx;

                u = v{k};
                xp = xpanc{k};
                yp = ypanc{k};

                vbot = u(1:idx);
                xbot = xp(1:idx);
                ybot = yp(1:idx);
                % pbot = cptemp(1:idx); % DEBUG
                vbot = flip(vbot);
                xbot = flip(xbot);
                ybot = flip(ybot);
                % pbot = flip(pbot); % DEBUG

                vtop = u(idx:end);
                xtop = xp(idx:end);
                ytop = yp(idx:end);
                % ptop = cptemp(idx:end); % DEBUG
                
                % DEBUG
                % figure
                % plot(xtop, ptop);
                % figure
                % plot(xbot, pbot);
                % figure
                % hold on
                % plot(xtop,ytop)
                % plot(xbot, ybot)
                % legend({'top' 'bottom'})
                % axis equal
                % xtop
                % ytop
                % xbot

                bldata{k,1} = [xtop'; ytop'; vtop'];
                bldata{k,2} = [xbot'; ybot'; vbot'];
                
            end

            obj.bl = bldata;

            disp(['Time elapsed for BL preprocessing: ' num2str(toc) ' seconds.'])

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



        function plotCp(obj)

            for k = 1:obj.nArfls
                figName = ['Cp - airfoil ' int2str(k)];
                figTitle = ['Airfoil ' int2str(k)];
                figure('Name', figName)
                plot(obj.xpanc{k}, obj.Cp{k}, '--', 'LineWidth',1.5, 'Color',[0, 0, 0.1])
                xlabel('x/c')
                ylabel('Cp')
                set(gca,'Ydir','reverse')
                title(figTitle)
                grid on
            end

        end



        function plotUe(obj)

            for k = 1:obj.nArfls
                figName = ['ue - airfoil ' int2str(k)];
                figTitle = ['Airfoil ' int2str(k)];
                figure('Name', figName)
                plot(obj.xpanc{k}, obj.ue{k}, '--', 'LineWidth',1.5, 'Color',[0, 0, 0.1])
                xlabel('x/c')
                ylabel('ue/U')
                title(figTitle)
                grid on
            end

        end


        function getBL(obj, Re, ht)

            for k=1:obj.nArfls

                % top
                blt = obj.bl{k,1};
                xt = blt(1,:);
                yt = blt(2,:);
                ut = blt(3,:);
                [warnOutT, x_transitionT, CfT] = solverBL(Re, xt, yt, ut, ht, 0, false);
                disp(['Top transition at x = ' num2str(x_transitionT)])

                % bottom
                blb = obj.bl{k,2};
                xb = blb(1,:);
                yb = blb(2,:);
                ub = blb(3,:);
                [warnOutB, x_transitionB, CfB] = solverBL(Re, xb, yb, ub, ht, 0, false);
                disp(['Bottom transition at x = ' num2str(x_transitionB)])

                if isempty(warnOutB) && isempty(warnOutT)
                    CfB = flip(CfB);
                    obj.Cf{k} = [CfB CfT];

                    % plotting
                    figName = ['Cf - airfoil ' int2str(k)];
                    figure('Name', figName)
                    hold on
                    grid on
                    plot(xb, CfB)
                    plot(xt, CfT)
                else
                    for ii = 1:length(warnOutT)
                        disp(warnOutT{ii})
                    end
                    for ii = 1:length(warnOutB)
                        disp(warnOutB{ii})
                    end
                end
            
            end

        end

    end
end

