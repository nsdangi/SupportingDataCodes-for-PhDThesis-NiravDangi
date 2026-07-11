clear 
clc
load('..\Data\3Dvelarray.mat_3Dup.mat')
addpath("spod_matlab-1.1.2")
%%

close all
for jj=1:length(INP_common.subfolders)
    subfolder=INP_common.subfolders{jj};
    disp(subfolder)
    

    X1_u=OPT_3dVel.(subfolder).vel_arr_U(:,:)';
    X1_v=OPT_3dVel.(subfolder).vel_arr_V(:,:)';
    X1_w=OPT_3dVel.(subfolder).vel_arr_W(:,:)';

    if jj == 4 || jj == 5 || jj == 6
        X1_u=X1_u(:,1:2100,:);
        X1_v=X1_v(:,1:2100,:);
        X1_w=X1_w(:,1:2100,:);
    end
    X1_u_avg = mean (X1_u ,2) ;

    X1_v_avg = mean (X1_v ,2) ;
   
    X1_w_avg = mean (X1_w ,2) ;
    if jj == 4 || jj == 5 || jj == 6
        dt=0.33;
        whatsizetimedim=2100;
        WINDOW=300;
        
        NOVERLAP=round(WINDOW/2);
    else
        dt=0.25;
        whatsizetimedim=2800;
        WINDOW=400;
        
        NOVERLAP=round(WINDOW/2);
    end
    
    [L_u,P_u,F_u,~,A_u] = spod(transpose(X1_u-X1_u_avg),WINDOW,[],NOVERLAP,dt); 
    [L_v,P_v,F_v,~,A_v] = spod(transpose(X1_v-X1_v_avg),WINDOW,[],NOVERLAP,dt);
    [L_w,P_w,F_w,~,A_w] = spod(transpose(X1_w-X1_w_avg),WINDOW,[],NOVERLAP,dt); 
    
    totalenergy=sum(L_u(:));
    
   
    
    
    
    
    dy=max(unique(diff(points.Y)));
    dz=max(unique(diff(points.Z)));
    x=points.Y(sorted_indices); % Cartesian X (east)     co-ordinate
    y=points.Z(sorted_indices); % Cartesian Y (north)    co-ordinate
    [Xgrid, Ygrid]=meshgrid(x,y); % For contour
    X_cont=Xgrid(1,:)';
    Y_cont=Ygrid(:,1);
    onearr=ones(height(X_cont),size(P_u,3));
   
    % figure
    count = 1;
    for fi = [ 2 3 4 14 ]
        for mi = [1 ]
            figure
            contourf(reshape(X_cont,35,[])-1505,reshape(Y_cont,35,[]),reshape(real(squeeze(P_u(fi,:,mi))),35,[]),50,'edgecolor','none')
            hold on
            contour(reshape(X_cont,35,[])-1505,reshape(Y_cont,35,[]),reshape(real(squeeze(P_u(fi,:,mi))),35,[]), 5, 'k','LineWidth',1.5) % 'k' for black color lines
            axis equal
            colorbar
            % colormap(flipud(gray(10)))
            colormap((gray(10)))
            centerX = 1505-1505;
            centerY = 170;
            radius = 142;
            turbine_center = [0, 170]; % Center of the turbine (x=0, y=170)
            hub_radius = 8; % Radius of the turbine hub
            blade_length = 140; % Length of each blade
            blade_width = 2; % Width of each blade

            draw_turbine = @(center) drawTurbine2D(center, hub_radius, blade_length, blade_width);
            draw_turbine(turbine_center);
            rectangle('Position',[centerX-radius, centerY-radius, 2*radius, 2*radius], 'Curvature', [1,1], 'EdgeColor', 'k','LineWidth',2);
            clim([-0.1 0.1])
            count = count + 1;
            xlabel('$y/D (-)$','Interpreter','latex'), ylabel('$z/D (-)$','Interpreter','latex') 
            xticks([-1:0.5:1].*280)
            xticklabels({'-1','-0.5','0','0.5','1'})
            yticks([0:0.25:1].*280)
            yticklabels({'0','0.25','0.5','0.75','1'})
            title(['$',strrep(subfolder,'_',' '),'$'],'Interpreter','latex')
            axis tight

            if fi==1
                subtitle(['$f < 0.01\,Hz, Energy=',num2str(100.*sum(L_u(fi,mi))./totalenergy,'%.2f'),'\%$'],'Interpreter','latex')%, mode ' num2str(mi) ', 
            else
                subtitle(['$f=', num2str(F_u(fi),'%.2f'),'\,Hz, Energy=',num2str(100.*sum(L_u(fi,:))./totalenergy,'%.2f'),'\%$'],'Interpreter','latex')%, mode ' num2str(mi) ', 
            end
            
            
            set(gca, 'FontSize', 14);
            ax = gca;
            ax.LineWidth=3;
            ax.TickLabelInterpreter="latex";
            box on
            axis equal
                
            if ~exist("..\Images", 'dir')
                mkdir("..\Images");
            end
            exportgraphics(gcf, ['..\Images\', subfolder, '_F_', num2str(fi), '.png']);
            
            close all
        end

    end
    
    clearvars -except INP_common OPT_3dVel points sorted_indices timesteps
    keyboard
end
%%
function drawTurbine2D(center, hub_radius, blade_length, blade_width)
    % Draw the hub
    rectangle('Position', [center(1)-hub_radius, center(2)-hub_radius, 2*hub_radius, 2*hub_radius], ...
              'Curvature', [1 1], 'FaceColor', 'k', 'EdgeColor', 'none');
    
    % Draw the blades
    angles = [0, 120, 240]; % Angles for the three blades
    for angle = angles
        % Calculate blade endpoints
        x_end = center(1) + blade_length * cosd(angle);
        y_end = center(2) + blade_length * sind(angle);
        % Draw the blade as a thick line
        plot([center(1), x_end], [center(2), y_end], 'k', 'LineWidth', blade_width);
    end
end
%%

