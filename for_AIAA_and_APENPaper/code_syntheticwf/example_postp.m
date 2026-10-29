clear
clc
addpath("Outputs\")
[velocity, y, z, nz, ny, dz, dy, dt, zHub, z1, SummVars] = func_readfile_WND('Outputs/example_owncoh.wnd');
timearr=0:dt:(size(velocity,1)-1)*dt;
%%

[YQ,ZQ]=meshgrid(y,z);
uatts100= transpose(squeeze(velocity(100,1,:,:)));
meanu = transpose(squeeze(mean(velocity(:,1,:,:),1)));


figure(1)
subplot(1,2,1)
contourf(YQ,ZQ,uatts100,100,EdgeColor="none")
title('Instantaneous u at 100th timestep')
xlabel('Y (m)')
ylabel('Z (m)')

subplot(1,2,2)
contourf(YQ,ZQ,meanu,100,EdgeColor="none")
title('Mean u ')
xlabel('Y (m)')
ylabel('Z (m)')

figure(2)
plot(timearr,squeeze(velocity(:,1,10,10)))
hold on
plot(timearr,squeeze(velocity(:,2,10,10)))
plot(timearr,squeeze(velocity(:,3,10,10)))
title('Instantaneous u,v,w at arbitrary location')
xlabel('Time (s)')
ylabel('u,v,w (m/s)')
