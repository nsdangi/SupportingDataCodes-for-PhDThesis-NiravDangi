%% DIVERGENCE
close all
Umean=[11.5];

custom=gray(10);
subtitlelabeldiv={'Original','Divergence- free'};
for b = 1:length(basePaths)
    figure(Visible="on")
    for c=1:2%length(cases)
        datahere=windData.(cases{c}).(basePaths{b}){1};

        [Nt, ~, Ny, Nz] = size(datahere);
        U_mean = Umean(b);   % Mean advection speed [m/s] (adjust as needed)
        dx = U_mean * dt(b);       % Spatial step in x
        Nx = Nt;                % Let the number of x points equal the number of timesteps
        x = (0:Nx-1) * dx;      % x coordinate vector

        % 3. Grid spacings in y and z
        dy = 10;
        dz = 10;
        
        
        velocity = buildVelocityField(datahere, Nx, Ny, Nz);

        div_v = computeDivergencefourth(velocity, Nx, Ny, Nz, dx, dy, dz);

        subplot(1,2,c)
        histogram(div_v(:),Normalization="probability",BinWidth=0.01)
        xlim([-0.15 0.15])
        ylim([0 0.6])
        xticks([-0.15:0.05:0.15])
        ylabel('$Probability [-]$','Interpreter','Latex')
        
        box on
        grid on
        grid minor
        xlabel('$\nabla \cdot V$ [1/s]','Interpreter','Latex')
        
        box on
        grid on
        grid minor

        ax = gca;
        ax.FontSize = 16;
        ax.TickLabelInterpreter = "latex";
        ax.LineWidth = 3;
       
        
        subtitle(subtitlelabeldiv{c},Interpreter="latex")
        

        
    end
    exportgraphics(gcf,['Images\DIVU_',cases{1},'_',basePaths{b},'.png'])
end
%% DIVERGENCE CONTOUR
close all
Umean=[11.5];


custom=flipud(gray);
subtitlelabeldiv={'Original','Divergence- free'};
for b = 1:length(basePaths)
    figure(Visible="on",Position=[488.0000  319.4000  705.8000  438.6000])
    for c=1:length(cases)
        datahere=windData.(cases{c}).(basePaths{b}){1};

        [Nt, ~, Ny, Nz] = size(datahere);
        U_mean = Umean(b);   % Mean advection speed [m/s] (adjust as needed)
        dx = U_mean * dt(b);       % Spatial step in x
        Nx = Nt;                % Let the number of x points equal the number of timesteps
        x = (0:Nx-1) * dx;      % x coordinate vector
        yyy=-150:10:150;
        zzz=25:10:325;
        [X, Y, Z] = meshgrid(yyy, zzz, x(500:600));
        % 3. Grid spacings in y and z
        dy = 10;
        dz = 10;
        
        
        velocity = buildVelocityField(datahere, Nx, Ny, Nz);

        div_v = computeDivergencefourth(velocity, Nx, Ny, Nz, dx, dy, dz);
        div_v_use=div_v(500:600,3:33,3:33,:); % to prevent RAM overload
        subplot(2,1,c)
        scatter3(Z(:)./280, X(:)./280, Y(:)./280, 20, abs(div_v_use(:)), 'filled')
        colormap(custom)
        cccc=colorbar;
        cccc.Ticks=[0:0.02:0.06];
        cccc.Label.String = '$\mid\nabla \cdot V\mid$ [1/s]'; 
        cccc.Label.Interpreter = 'latex';
        cccc.TickLabelInterpreter= 'latex';
        zticks([0:0.25:1])
        zticklabels({'0','0.25','0.5','0.75','1'})
        yticks([-1:0.5:1])
        yticklabels({'-1','-0.5','0','0.5','1'})
        xticks([0:0.25:10])
        if c==2
            
            zlabel('$z/D [-]$', Interpreter='latex')
            
            ylabel('$y/D [-]$', Interpreter='latex')
            
            xlabel('$x/D [-]$', Interpreter='latex')
            
        end
        clim([0 0.05])
        box on
        grid on
        grid minor

        ax = gca;
        ax.FontSize = 16;
        ax.TickLabelInterpreter = "latex";
        ax.LineWidth = 3;
        axis tight
        
        subtitle(subtitlelabeldiv{c},Interpreter="latex")
        

        
    end
   exportgraphics(gcf,['Images\DIVUCONT_',cases{1},'_',basePaths{b},'.png'])
end

function div_v = computeDivergencefourth(velocity, Nx, Ny, Nz, dx, dy, dz)
    div_v = zeros(Nx, Ny, Nz);
    for k = 1:Nz
        for j = 1:Ny
            for i = 1:Nx
                % du/dx (Fourth-order centered difference)
                if i <= 2 || i >= Nx-1
                    du_dx = 0; % Boundary condition for first and last points
                else
                    du_dx = (-velocity(i+2, j, k, 1) + 8*velocity(i+1, j, k, 1) - 8*velocity(i-1, j, k, 1) + velocity(i-2, j, k, 1)) / (12 * dx);
                end

                % dv/dy (Fourth-order periodic difference)
                jp = mod(j, Ny) + 1;
                jm = mod(j-2, Ny) + 1;
                jp2 = mod(j+1, Ny) + 1;
                jm2 = mod(j-3, Ny) + 1;
                dv_dy = (-velocity(i, jp2, k, 2) + 8*velocity(i, jp, k, 2) - 8*velocity(i, jm, k, 2) + velocity(i, jm2, k, 2)) / (12 * dy);

                % dw/dz (Fourth-order with no-flow BC)
                if k <= 2 || k >= Nz-1
                    dw_dz = 0;
                else
                    dw_dz = (-velocity(i, j, k+2, 3) + 8*velocity(i, j, k+1, 3) - 8*velocity(i, j, k-1, 3) + velocity(i, j, k-2, 3)) / (12 * dz);
                end

                div_v(i, j, k) = du_dx + dv_dy + dw_dz;
            end
        end
    end
end
function velocity = buildVelocityField(velocity2D, Nx, Ny, Nz)
velocity = zeros(Nx, Ny, Nz, 3);
velocity(:,:,:,1) = velocity2D(:,1,:,:);  % u component (x-direction)
velocity(:,:,:,2) = velocity2D(:,2,:,:);  % v component (y-direction)
velocity(:,:,:,3) = velocity2D(:,3,:,:);  % w component (z-direction)
end
