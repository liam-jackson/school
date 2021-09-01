%Numerical Modeling of the Male Testosterone Cycle

%Daniel Carbonero, Liam Jackson, Noshin Nawar, & Kathryn Regan
close all; 
clear all; 
clc;

%% Healthy Case
%Parameter Initializations
history = ...       %Initial Hormone Concentrations:
    [0;             %LHRH
    10;             %LH
    12];            %T

tspan = [0:1440];    %1440 mins = 24 hrs interval. 

%Set Time Delays
Thp = 3; 
Tpt = 5; 
Tth = 5; 
Tph = 5; 
T0 = 25;
lags = [Tph Tth Thp Tpt+T0];

%Hormone System Simulation 
sol = dde23(@ddefun,lags,history,tspan);

%Interpolation for Phase Plots
t_res = 1; %use this value to control the phase plot resolutions
t_grid = min(tspan):t_res:max(tspan);

for i = 1:length(sol.y(:,1))
    solY_grid(i,:) = interp1(sol.x, sol.y(i,:), t_grid); 
    solYp_grid(i,:) = interp1(sol.x, sol.yp(i,:), t_grid); 
end

%Data Plotting:
%All Hormones vs time
fig1 = figure(1);
    sgtitle('Hormone Levels During Cyclic Testosterone Production');
sub1fig1 = subplot(2,1,1);
    plot(sol.x, sol.y, 'Linewidth',1)
    grid minor; 
    legend('[LHRH]','[LH]','[T]', 'FontSize',12); 
    xlabel('Time (minutes)'); 
    ylabel('Concentration (ng/mL)');

%LHRH vs time     
sub2fig1 = subplot(2,1,2);
    nano2pico = sol.y(1,:).*1000;
    plot(sol.x, nano2pico, 'Linewidth',1) 
    grid minor; 
    legend('[LHRH]', 'FontSize',12)
    xlabel('Time (minutes)'); 
    ylabel('Concentration (pg/mL)');

%3D plot of data
%Note: a box with 2 sliders will pop up, which controls the dynamic label
%positioning of the 3D figure, aligning the labels with the respective axes
fig2 = figure(2);    
    ax = axes('parent',fig2);
    colormap(jet);
    c = linspace(1,length(t_grid),length(t_grid));
    datYfig2 = scatter3(ax,solY_grid(1,:),solY_grid(2,:),solY_grid(3,:),[],c,'*'); 
    hold on;
    datYICfig2 = scatter3(ax,solY_grid(1,1),solY_grid(2,1),solY_grid(3,1),100,c(1,1),'filled');
    grid on; 
    title(ax,'Phase Plot of [C] vs. [C]','FontSize',18);
    legend([datYfig2 datYICfig2],'[C]','[C]|t=0');
    set(gca, 'dataaspectratio',...
       [max(max(max(solY_grid(1,:)))),...
        max(max(max(solY_grid(2,:)))),...
        max(max(max(solY_grid(3,:))))]);
    xlabel('Concentration LHRH (ng/mL)','FontSize',15); 
    ylabel('Concentration LH (ng/mL)','FontSize',15);
    zlabel('Concentration T (ng/mL)','FontSize',15);
    h = rotate3d; 
    set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)') 
    set(gcf, 'ResizeFcn', @align_axislabel) 
    align_axislabel([], gca) 
    axislabel_translation_slider;
    c5 = colorbar;
    title(c5,'Time (min)')
    hold off;

%Projections of Figure 5
fig3 = figure(3);
    sgtitle('Phase Plots of Hormone Concentrations','FontSize',18);
sub1fig3 = subplot(1,3,1); 
    colormap(jet);
    c = linspace(1,length(t_grid),length(t_grid));
    datYsub1fig3 = scatter3(solY_grid(1,:),...
        solY_grid(2,:),solY_grid(3,:),[],c,'*'); 
    hold on;
    datYICsub1fig3 = scatter3(solY_grid(1,1),...
        solY_grid(2,1),solY_grid(3,1),100,c(1,1),'filled');
    grid on; 
    legend([datYsub1fig3 datYICsub1fig3],'[C]','[C]|t=0','Location','northeast');
    xlabel('Concentration LHRH (ng/mL)','FontSize',15); 
    ylabel('Concentration LH (ng/mL)','FontSize',15);
    zlabel('Concentration T (ng/mL)','FontSize',15);
    hold off;

sub2fig3 = subplot(1,3,2); 
    datYsub2fig3 = scatter3(solY_grid(1,:),...
        solY_grid(2,:),solY_grid(3,:),[],c,'*'); 
    hold on;
    datYICsub2fig3 = scatter3(solY_grid(1,1),...
        solY_grid(2,1),solY_grid(3,1),100,c(1,1),'filled');
    grid on; 
    legend([datYsub2fig3 datYICsub2fig3],'[C]','[C]|t=0','Location','northeast');
    xlabel('Concentration LHRH (ng/mL)','FontSize',15); 
    ylabel('Concentration LH (ng/mL)','FontSize',15);
    zlabel('Concentration T (ng/mL)','FontSize',15);
    hold off;
    
sub3fig3 = subplot(1,3,3); 
    datYsub3fig3 = scatter3(solY_grid(1,:),...
        solY_grid(2,:),solY_grid(3,:),[],c,'*'); 
    hold on;
    datYICsub3fig3 = scatter3(solY_grid(1,1),...
        solY_grid(2,1),solY_grid(3,1),100,c(1,1),'filled');
    grid on; 
    legend([datYsub3fig3 datYICsub3fig3],'[C]','[C]|t=0','Location','northeast');
    xlabel('Concentration LHRH (ng/mL)','FontSize',15); 
    ylabel('Concentration LH (ng/mL)','FontSize',15);
    zlabel('Concentration T (ng/mL)','FontSize',15);
    c6 = colorbar;
    title(c6,'Time (min)')
    hold off;

%Change Camera Views to show Planar Projections
    view(sub1fig3,0,90); 
    view(sub2fig3,90,0);    
    view(sub3fig3,0,0);         

%% Castrated Case
clear sol;

%Parameter Initializations
history = ...       %Initial Hormone Concentrations:
    [0;             %LHRH
    10;             %LH
    12];            %T

tspan = [0:1440];    %1440 mins = 24 hrs interval. 

%Set Time Delays
Thp = 3; 
Tpt = 5; 
Tth = 5; 
Tph = 5; 
T0 = 25;
lags = [Tph Tth Thp Tpt+T0];

%Hormone System Simulation 
sol = dde23(@ddefun_cast,lags,history,tspan);

% %Interpolation for Phase Plots
% t_res = 1; %use this value to control the phase plot resolutions
% t_grid = min(tspan):t_res:max(tspan);

for i = 1:length(sol.y(:,1))
    solY_grid(i,:) = interp1(sol.x, sol.y(i,:), t_grid); 
    solYp_grid(i,:) = interp1(sol.x, sol.yp(i,:), t_grid); 
end

%Data Plotting:
%All Hormones vs time
fig4 = figure(4);
    sgtitle({'Hormone Levels During Cyclic Testosterone Production','(Castration)'});
sub1fig4 = subplot(2,1,1);
    plot(sol.x, sol.y, 'Linewidth',1)
    grid minor; 
    legend('[LHRH]','[LH]','[T]', 'FontSize',12); 
    xlabel('Time (minutes)'); 
    ylabel('Concentration (ng/mL)');

%LHRH vs time     
sub2fig4 = subplot(2,1,2);
    nano2pico = sol.y(1,:).*1000;
    plot(sol.x, nano2pico, 'Linewidth',1) 
    grid minor; 
    legend('[LHRH]', 'FontSize',12)
    xlabel('Time (minutes)'); 
    ylabel('Concentration (pg/mL)');

%3D plot of data
fig5 = figure(5);    
    ax = axes('parent',fig5);
    colormap(jet);
    c = linspace(1,length(t_grid),length(t_grid));
    datYfig5 = scatter3(ax,solY_grid(1,:),solY_grid(2,:),solY_grid(3,:),[],c,'*'); 
    hold on;
    datYICfig5 = scatter3(ax,solY_grid(1,1),solY_grid(2,1),solY_grid(3,1),100,c(1,1),'filled');
    grid on; 
    title(ax,{'Phase Plot of [C] vs. [C]','(Castration)'}, 'FontSize',18);
    legend([datYfig5 datYICfig5],'[C]','[C]|t=0');
    set(gca, 'dataaspectratio',...
       [max(max(max(solY_grid(1,:)))),...
        max(max(max(solY_grid(2,:)))),...
        max(max(max(solY_grid(3,:))))]);
    xlabel('Concentration LHRH (ng/mL)','FontSize',15); 
    ylabel('Concentration LH (ng/mL)','FontSize',15);
    zlabel('Concentration T (ng/mL)','FontSize',15);
    h = rotate3d; 
    set(h, 'ActionPreCallback', 'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)') 
    set(gcf, 'ResizeFcn', @align_axislabel) 
    align_axislabel([], gca) 
    axislabel_translation_slider;
    c5 = colorbar;
    title(c5,'Time (min)')
    hold off;

%Projections of Figure 5
fig6 = figure(6);
sgtitle({'Phase Plots of Hormone Concentrations','(Castration)'}, 'FontSize',18);
sub1fig6 = subplot(1,3,1); 
    colormap(jet);
    c = linspace(1,length(t_grid),length(t_grid));
    datYsub1fig6 = scatter3(solY_grid(1,:),...
        solY_grid(2,:),solY_grid(3,:),[],c,'*'); 
    hold on;
    datYICsub1fig6 = scatter3(solY_grid(1,1),...
        solY_grid(2,1),solY_grid(3,1),100,c(1,1),'filled');
    grid on; 
    legend([datYsub1fig6 datYICsub1fig6],'[C]','[C]|t=0','Location','northeast');
    xlabel('Concentration LHRH (ng/mL)','FontSize',15); 
    ylabel('Concentration LH (ng/mL)','FontSize',15);
    zlabel('Concentration T (ng/mL)','FontSize',15);
    hold off;

sub2fig6 = subplot(1,3,2); 
    datYsub2fig6 = scatter3(solY_grid(1,:),...
        solY_grid(2,:),solY_grid(3,:),[],c,'*'); 
    hold on;
    datYICsub2fig6 = scatter3(solY_grid(1,1),...
        solY_grid(2,1),solY_grid(3,1),100,c(1,1),'filled');
    grid on; 
    legend([datYsub2fig6 datYICsub2fig6],'[C]','[C]|t=0','Location','northeast');
    xlabel('Concentration LHRH (ng/mL)','FontSize',15); 
    ylabel('Concentration LH (ng/mL)','FontSize',15);
    zlabel('Concentration T (ng/mL)','FontSize',15);
    hold off;
    
sub3fig6 = subplot(1,3,3); 
    datYsub3fig6 = scatter3(solY_grid(1,:),...
        solY_grid(2,:),solY_grid(3,:),[],c,'*'); 
    hold on;
    datYICsub3fig6 = scatter3(solY_grid(1,1),...
        solY_grid(2,1),solY_grid(3,1),100,c(1,1),'filled');
    grid on; 
    legend([datYsub3fig6 datYICsub3fig6],'[C]','[C]|t=0','Location','northeast');
    xlabel('Concentration LHRH (ng/mL)','FontSize',15); 
    ylabel('Concentration LH (ng/mL)','FontSize',15);
    zlabel('Concentration T (ng/mL)','FontSize',15);
    c6 = colorbar;
    title(c6,'Time (min)')
    hold off;

view(sub1fig6,0,90); 
view(sub2fig6,90,0);    
view(sub3fig6,0,0);         

%%
function dydt = ddefun(t,y,Z)
%Set constants
rR =0.1; rL =5; rT =0.008; %0.01;       
dR =0.1; dL =0.015;  dT =0.04;  %1/min

ylag1 = Z(:,1); ylag2 = Z(:,2); ylag3 = Z(:,3); ylag4 = Z(:,4);

%Define Heaviside Function
L =30; T =8;
Hx = (2 - (ylag1(2)/L) - (ylag2(3)/T));
if Hx > 0
    H = 1;
elseif Hx == 0
    H = 0.5;
else
    H = 0;
end

%set differential equations (order: [R];[L];[T])
dydt = [  -dR*y(1) + rR*H         
          -dL*y(2) + rL*ylag3(1)                                    
          -dT*y(3) + rT*ylag4(2)      ];
end 

function dydt = ddefun_cast(t,y,Z)
%Set constants
rR =0.1; rL =5; rT =0.008; %0.01;       
dR =0.1; dL =0.015;  dT =0.04;  %1/min

ylag1 = Z(:,1); ylag2 = Z(:,2); ylag3 = Z(:,3); ylag4 = Z(:,4);

%Define Heaviside Function
L =30; T =8;
if t < 600
    Hx = (2 - (ylag1(2)/L) - (ylag2(3)/T));
    castration=1;
else
    Hx = (2 - (ylag1(2)/L)); 
    castration=0;
end

if Hx > 0
    H = 1;
elseif Hx == 0;
    H = 0.5;
else
    H = 0;
end

%set differential equations (order: [R];[L];[T])
dydt = [  -dR*y(1) + rR*H      
          -dL*y(2) + rL*ylag3(1)                                    
          -dT*y(3) + rT*ylag4(2)*castration];
end 

%The following is a function that I included here because it was not native
%to Matlab, but was used to align the axes labels in the 3D plots.
%Definitely not worth scrolling through. 

function align_axislabel(~,ax)
% Note: in order to maintain graphical neatness, I found this function on
% the mathworks website. 

% This function first rotates x, y and z labels to the direction of their 
% corresponding axes when a rotate operation has finished. So that labels 
% will be parallel to axes. Then, moves the axis labels to a proper 
% distance from the axes by calling 'axislabel_translation'.
% The function is used as an 'ActionPostCallback', which is a callback 
% function of a rotate3d object.
% It still works when the projection mode is perspective or when the data 
% aspect ratio is not [1 1 1].
%
% Example 1: % set ActionPostCallback to align_axislabel
%  h = rotate3d;
%  set(h,'ActionPostCallback',@align_axislabel);
%
% Example 2: % set WindowButtonMotionFcn to align_axislabel
%  h = rotate3d;
%  set(h,'ActionPreCallback',...
%      'set(gcf,''windowbuttonmotionfcn'',@align_axislabel)')
%  set(h,'ActionPostCallback',...
%      'set(gcf,''windowbuttonmotionfcn'','''')')
%
% Jul/13/2015, Adds a trans_mode option
% Feb/06/2015
% Apr/10/2016

% trans_mode controls the way we define trans_vec_x and trans_vec_y, which
% is an input argument passed to function axislabel_translation.
% Take x-axis for example. If trans_mode == 2, the translation direction is
% parallel to y-axis, and if trans_moede == 1, the translation direction is
% perpendicular to x-axis.
trans_mode = 1;

if isa(ax,'struct')
    h_a = ax.Axes;
elseif isa(ax,'matlab.graphics.axis.Axes')
    h_a = ax;
else
    h_a = gca;
end

proj_mode = get(h_a,'projection');
if strcmpi(proj_mode,'orthographic')
    proj_mode = 1;
else
    proj_mode = 2;
end
data_aspect_ratio = get(h_a,'dataaspectratio');
data_aspect_ratio = data_aspect_ratio(:);

set(h_a, 'DataAspectRatio', data_aspect_ratio); % optional, you can comment this line out

x_lim = get(h_a,'xlim');
y_lim = get(h_a,'ylim');
z_lim = get(h_a,'zlim');
cam_pos = get(h_a,'cameraposition');
cam_tar = get(h_a,'cameratarget');
cam_vec = cam_pos-cam_tar;
N = norm(cam_vec,2); % Near field
[az,el] = view; % get current azimuth and elevation angles
if abs(el) == 90
    az = atan2d(sind(az)*data_aspect_ratio(2),cosd(az)*data_aspect_ratio(1));
end % Adjust az when el == 90
az = az/180*pi;
el = pi/2-el/180*pi;
% Rotation Matrices
R_az = [cos(-az),-sin(-az),0,0;sin(-az),cos(-az),0,0;0,0,1,0;0,0,0,1];
R_el = [1,0,0,0;0,cos(-el),-sin(-el),0;0,sin(-el),cos(-el),0;0,0,0,1];
% Translation Matrices
T_tar = [1,0,0,-cam_tar(1);0,1,0,-cam_tar(2);0,0,1,-cam_tar(3);0,0,0,1];
T_cam = [1,0,0,0;0,1,0,0;0,0,1,-N;0,0,0,1];
pnt = [x_lim(1),y_lim(1),z_lim(1),1;...
    x_lim(2),y_lim(1),z_lim(1),1;...
    x_lim(2),y_lim(2),z_lim(1),1;...
    x_lim(1),y_lim(2),z_lim(1),1;...
    x_lim(1),y_lim(1),z_lim(2),1;...
    x_lim(2),y_lim(1),z_lim(2),1;...
    x_lim(2),y_lim(2),z_lim(2),1;...
    x_lim(1),y_lim(2),z_lim(2),1]';
pnt_cam = zeros(4,8);
if proj_mode == 1
    F = N+100; % Far field, the value would not influence the projection
    a = -(F+N)/(F-N);
    b = -2*F*N/(F-N);
    M_orth = [N,0,0,0;0,N,0,0;0,0,a,b;0,0,-1,0];
    pnt_cam = zeros(4,4);
    for kk = 1:8
        p = pnt(:,kk);
        p = T_tar*p;
        p(1:3) = p(1:3)./data_aspect_ratio;
        p = R_el*R_az*p;
        p = T_cam*p;
        pnt_cam(:,kk) = M_orth*p;
    end
    x_vec = pnt_cam(:,2)-pnt_cam(:,1);
    y_vec = pnt_cam(:,4)-pnt_cam(:,1);
    z_vec = pnt_cam(:,5)-pnt_cam(:,1);
    % find the most-left point
    [~,ix] = min(pnt_cam(1,:));
    % find the lowest point
    [~,iy] = min(pnt_cam(2,:));
    switch ix % use ix to determine z axis, then use iy to determine x and y axis
        case {1,5}
            trans_pnt_z = [x_lim(1),y_lim(1),(z_lim(1)+z_lim(2))/2]';
            if iy == 2
                trans_vec_x_edge = pnt_cam(1:2,2)-pnt_cam(1:2,3);
                trans_vec_y_edge = pnt_cam(1:2,2)-pnt_cam(1:2,1);
                trans_pnt_x = (pnt(1:3,2)+pnt(1:3,1))/2;
                trans_pnt_y = (pnt(1:3,3)+pnt(1:3,2))/2;
            else % iy == 4
                trans_vec_x_edge = pnt_cam(1:2,4)-pnt_cam(1:2,1);
                trans_vec_y_edge = pnt_cam(1:2,4)-pnt_cam(1:2,3);
                trans_pnt_y = (pnt(1:3,4)+pnt(1:3,1))/2;
                trans_pnt_x = (pnt(1:3,3)+pnt(1:3,4))/2;
            end
        case {2,6}
            trans_pnt_z = [x_lim(2),y_lim(1),(z_lim(1)+z_lim(2))/2]';
            if iy == 3
                trans_vec_x_edge = pnt_cam(1:2,3)-pnt_cam(1:2,2);
                trans_vec_y_edge = pnt_cam(1:2,3)-pnt_cam(1:2,4);
                trans_pnt_y = (pnt(1:3,3)+pnt(1:3,2))/2;
                trans_pnt_x = (pnt(1:3,4)+pnt(1:3,3))/2;
            else % iy == 1
                trans_vec_x_edge = pnt_cam(1:2,1)-pnt_cam(1:2,4);
                trans_vec_y_edge = pnt_cam(1:2,1)-pnt_cam(1:2,2);
                trans_pnt_x = (pnt(1:3,1)+pnt(1:3,2))/2;
                trans_pnt_y = (pnt(1:3,4)+pnt(1:3,1))/2;
            end
        case {3,7}
            trans_pnt_z = [x_lim(2),y_lim(2),(z_lim(1)+z_lim(2))/2]';
            if iy == 4
                trans_vec_x_edge = pnt_cam(1:2,4)-pnt_cam(1:2,1);
                trans_vec_y_edge = pnt_cam(1:2,4)-pnt_cam(1:2,3);
                trans_pnt_x = (pnt(1:3,4)+pnt(1:3,3))/2;
                trans_pnt_y = (pnt(1:3,1)+pnt(1:3,4))/2;
            else % iy == 2
                trans_vec_x_edge = pnt_cam(1:2,2)-pnt_cam(1:2,3);
                trans_vec_y_edge = pnt_cam(1:2,2)-pnt_cam(1:2,1);
                trans_pnt_y = (pnt(1:3,2)+pnt(1:3,3))/2;
                trans_pnt_x = (pnt(1:3,1)+pnt(1:3,2))/2;
            end
        case {4,8}
            trans_pnt_z = [x_lim(1),y_lim(2),(z_lim(1)+z_lim(2))/2]';
            if iy == 1
                trans_vec_x_edge = pnt_cam(1:2,1)-pnt_cam(1:2,4);
                trans_vec_y_edge = pnt_cam(1:2,1)-pnt_cam(1:2,2);
                trans_pnt_y = (pnt(1:3,1)+pnt(1:3,4))/2;
                trans_pnt_x = (pnt(1:3,2)+pnt(1:3,1))/2;
            else % iy == 3
                trans_vec_x_edge = pnt_cam(1:2,3)-pnt_cam(1:2,2);
                trans_vec_y_edge = pnt_cam(1:2,3)-pnt_cam(1:2,4);
                trans_pnt_x = (pnt(1:3,3)+pnt(1:3,4))/2;
                trans_pnt_y = (pnt(1:3,2)+pnt(1:3,3))/2;
            end
    end
else
    for kk = 1:8
        p = pnt(:,kk);
        p = T_tar*p;
        p(1:3) = p(1:3)./data_aspect_ratio;
        p = R_el*R_az*p;
        p = T_cam*p;
        z = p(3);
        M_pers = [-N/z,0,0,0;0,-N/z,0,0;0,0,-N,0;0,0,0,1];
        pnt_cam(:,kk) = M_pers*p;
    end
    % find the most-left point
    [~,ix] = min(pnt_cam(1,:));
    % find the lowest point
    [~,iy] = min(pnt_cam(2,:));
    % find the proper axis to label
    switch ix % use ix to determine z axis, then use iy to determine x and y axis
        case {1,5}
            z_vec = pnt_cam(:,5)-pnt_cam(:,1);
            trans_pnt_z = [x_lim(1),y_lim(1),(z_lim(1)+z_lim(2))/2]';
            if iy == 2
                x_vec = pnt_cam(:,2)-pnt_cam(:,1);
                y_vec = pnt_cam(:,3)-pnt_cam(:,2);
                trans_vec_x_edge = pnt_cam(1:2,2)-pnt_cam(1:2,3);
                trans_vec_y_edge = pnt_cam(1:2,2)-pnt_cam(1:2,1);
                trans_pnt_x = (pnt(1:3,2)+pnt(1:3,1))/2;
                trans_pnt_y = (pnt(1:3,3)+pnt(1:3,2))/2;
            else % iy == 4
                y_vec = pnt_cam(:,4)-pnt_cam(:,1);
                x_vec = pnt_cam(:,3)-pnt_cam(:,4);
                trans_vec_x_edge = pnt_cam(1:2,4)-pnt_cam(1:2,1);
                trans_vec_y_edge = pnt_cam(1:2,4)-pnt_cam(1:2,3);
                trans_pnt_y = (pnt(1:3,4)+pnt(1:3,1))/2;
                trans_pnt_x = (pnt(1:3,3)+pnt(1:3,4))/2;
            end
        case {2,6}
            z_vec = pnt_cam(:,6)-pnt_cam(:,2);
            trans_pnt_z = [x_lim(2),y_lim(1),(z_lim(1)+z_lim(2))/2]';
            if iy == 3
                y_vec = pnt_cam(:,3)-pnt_cam(:,2);
                x_vec = pnt_cam(:,4)-pnt_cam(:,3);
                trans_vec_x_edge = pnt_cam(1:2,3)-pnt_cam(1:2,2);
                trans_vec_y_edge = pnt_cam(1:2,3)-pnt_cam(1:2,4);
                trans_pnt_y = (pnt(1:3,3)+pnt(1:3,2))/2;
                trans_pnt_x = (pnt(1:3,4)+pnt(1:3,3))/2;
            else % iy == 1
                x_vec = pnt_cam(:,1)-pnt_cam(:,2);
                y_vec = pnt_cam(:,4)-pnt_cam(:,1);
                trans_vec_x_edge = pnt_cam(1:2,1)-pnt_cam(1:2,4);
                trans_vec_y_edge = pnt_cam(1:2,1)-pnt_cam(1:2,2);
                trans_pnt_x = (pnt(1:3,1)+pnt(1:3,2))/2;
                trans_pnt_y = (pnt(1:3,4)+pnt(1:3,1))/2;
            end
        case {3,7}
            z_vec = pnt_cam(:,7)-pnt_cam(:,3);
            trans_pnt_z = [x_lim(2),y_lim(2),(z_lim(1)+z_lim(2))/2]';
            if iy == 4
                x_vec = pnt_cam(:,4)-pnt_cam(:,3);
                y_vec = pnt_cam(:,1)-pnt_cam(:,4);
                trans_vec_x_edge = pnt_cam(1:2,4)-pnt_cam(1:2,1);
                trans_vec_y_edge = pnt_cam(1:2,4)-pnt_cam(1:2,3);
                trans_pnt_x = (pnt(1:3,4)+pnt(1:3,3))/2;
                trans_pnt_y = (pnt(1:3,1)+pnt(1:3,4))/2;
            else % iy == 2
                y_vec = pnt_cam(:,2)-pnt_cam(:,3);
                x_vec = pnt_cam(:,1)-pnt_cam(:,2);
                trans_vec_x_edge = pnt_cam(1:2,2)-pnt_cam(1:2,3);
                trans_vec_y_edge = pnt_cam(1:2,2)-pnt_cam(1:2,1);
                trans_pnt_y = (pnt(1:3,2)+pnt(1:3,3))/2;
                trans_pnt_x = (pnt(1:3,1)+pnt(1:3,2))/2;
            end
        case {4,8}
            z_vec = pnt_cam(:,8)-pnt_cam(:,4);
            trans_pnt_z = [x_lim(1),y_lim(2),(z_lim(1)+z_lim(2))/2]';
            if iy == 1
                y_vec = pnt_cam(:,1)-pnt_cam(:,4);
                x_vec = pnt_cam(:,2)-pnt_cam(:,1);
                trans_vec_x_edge = pnt_cam(1:2,1)-pnt_cam(1:2,4);
                trans_vec_y_edge = pnt_cam(1:2,1)-pnt_cam(1:2,2);
                trans_pnt_y = (pnt(1:3,1)+pnt(1:3,4))/2;
                trans_pnt_x = (pnt(1:3,2)+pnt(1:3,1))/2;
            else % iy == 3
                x_vec = pnt_cam(:,3)-pnt_cam(:,4);
                y_vec = pnt_cam(:,2)-pnt_cam(:,3);
                trans_vec_x_edge = pnt_cam(1:2,3)-pnt_cam(1:2,2);
                trans_vec_y_edge = pnt_cam(1:2,3)-pnt_cam(1:2,4);
                trans_pnt_x = (pnt(1:3,3)+pnt(1:3,4))/2;
                trans_pnt_y = (pnt(1:3,2)+pnt(1:3,3))/2;
            end
    end
end

if trans_mode == 1
    trans_vec_x = [x_vec(2),-x_vec(1)]';
    trans_vec_y = [y_vec(2),-y_vec(1)]';
    trans_vec_z = [-z_vec(2),z_vec(1)]';
    if trans_vec_x'*trans_vec_x_edge < 0
        trans_vec_x = -trans_vec_x;
    end
    if trans_vec_y'*trans_vec_y_edge < 0
        trans_vec_y = -trans_vec_y;
    end
elseif trans_mode == 2
    trans_vec_x = trans_vec_x_edge;
    trans_vec_y = trans_vec_y_edge;
    trans_vec_z = [-z_vec(2),z_vec(1)]';
end

% Normalize translation vectors
trans_vec_x = trans_vec_x/norm(trans_vec_x,2);
trans_vec_y = trans_vec_y/norm(trans_vec_y,2);
trans_vec_z = trans_vec_z/norm(trans_vec_z,2);

% Compute rotation angles
theta_x = atan2d(x_vec(2),x_vec(1));
theta_y = atan2d(y_vec(2),y_vec(1));
theta_z = atan2d(z_vec(2),z_vec(1));
if abs(theta_x) >= 90
    theta_x = theta_x+180;
end
if abs(theta_y) >= 90
    theta_y = theta_y+180;
end
% if abs(theta_z) >= 90
%     theta_z = theta_z+180;
% end
set(get(h_a,'xlabel'),'rotation',theta_x);
set(get(h_a,'ylabel'),'rotation',theta_y);
set(get(h_a,'zlabel'),'rotation',theta_z);

%% Axis label translation
axislabel_translation(h_a,...
    [trans_pnt_x,trans_pnt_y,trans_pnt_z],...
    [trans_vec_x,trans_vec_y,trans_vec_z],...
    'character')
end

function axislabel_translation(h_a,trans_pnt,trans_vec,trans_units)
% This function moves axis labels to a proper distance from the axes. It is
% called by function 'align_axislabels', which first, rotates the x, y and 
% z labels.
% If you just want to move the labels but not rotate them, you have to 
% specify the position where to put the labels and the direction along 
% which you move labels away from the axes.
%
% Input arguments: 
%  trans_units: 'pixels', move labels away from axes in pixels;
%  'characters', move labels away from axes in characters;
%  trans_vec: translation vector, a 2-by-3 matrix, each column defines the 
%  direction along which you move labels away from the axes, the vector 
%  should be normalized;
%  trans_pnt: translation point, a 3-by-3 matrix, each column defines the 
%  position (in the original 3D space rather than in the 2D canvas space) 
%  where you put the axis labels;
%  h_a: handle of the axes.
%
% Apr/02/2015

global AXISALIGN_TRANS_A AXISALIGN_TRANS_B
if isempty(AXISALIGN_TRANS_A), AXISALIGN_TRANS_A = 1; end
if isempty(AXISALIGN_TRANS_B), AXISALIGN_TRANS_B = 1; end

if nargin < 4
    trans_units = 'pixels';
end
vers = version();
vers = str2double(vers(1:3));
h_xlb = get(h_a,'xlabel');
h_ylb = get(h_a,'ylabel');
h_zlb = get(h_a,'zlabel');
if strcmpi(trans_units,'pixels')
    trans_units = 1;
    % Modify the method by which you move labels here
    %======================================%
    char_pixel_len = get(h_xlb,'FontSize')*0.55; % 0.67
    base_pixel_len = get(h_xlb,'FontSize')*2.0; % 2.0
    %======================================%
elseif strcmpi(trans_units,'character')
    trans_units = 2;
    % Modify the method by which you move labels here
    %======================================%
    % pos = pos+trans_vec'.*aniso_ratio*(A*char_len+B)
    A = AXISALIGN_TRANS_A;
    B = AXISALIGN_TRANS_B;
    aniso_ratio = [2, 1];
    %======================================%
end
%% Move XLabel
set(h_xlb,'HorizontalAlignment','center','VerticalAlignment','Middle','Units','data');
set(h_xlb,'Position',trans_pnt(:,1))
char_len_x = label_max_len('x');
if trans_units == 2
    set(h_xlb,'units','characters');
    pos_xlb = get(h_xlb,'Position');
    pos_xlb(1:2) = pos_xlb(1:2)+trans_vec(:,1)'.*aniso_ratio*(A*char_len_x+B);
else
    set(h_xlb,'Units','Pixels');
    pos_xlb = get(h_xlb,'Position');
    pos_xlb(1:2) = pos_xlb(1:2)+trans_vec(:,1)'*(base_pixel_len+char_len_x*char_pixel_len);
end
set(h_xlb,'Position',pos_xlb)
%% Move YLabel
set(h_ylb,'HorizontalAlignment','center','VerticalAlignment','Middle','Units','data');
set(h_ylb,'Position',trans_pnt(:,2))
char_len_y = label_max_len('y');
if trans_units == 2
    set(h_ylb,'units','characters');
    pos_ylb = get(h_ylb,'Position');
    pos_ylb(1:2) = pos_ylb(1:2)+trans_vec(:,2)'.*aniso_ratio*(A*char_len_y+B);
else
    set(h_ylb,'Units','Pixels');
    pos_ylb = get(h_ylb,'Position');
    pos_ylb(1:2) = pos_ylb(1:2)+trans_vec(:,2)'*(base_pixel_len+char_len_y*char_pixel_len);
end
set(h_ylb,'Position',pos_ylb)
%% Move ZLabel
set(h_zlb,'HorizontalAlignment','center','VerticalAlignment','Middle','Units','data');
set(h_zlb,'Position',trans_pnt(:,3))
char_len_z = label_max_len('z');
if trans_units == 2
    set(h_zlb,'units','characters');
    pos_zlb = get(h_zlb,'Position');
    pos_zlb(1:2) = pos_zlb(1:2)+trans_vec(:,3)'.*aniso_ratio*(A*char_len_z+B);
else
    set(h_zlb,'Units','Pixels');
    pos_zlb = get(h_zlb,'Position');
    pos_zlb(1:2) = pos_zlb(1:2)+trans_vec(:,3)'*(base_pixel_len+char_len_z*char_pixel_len);
end
set(h_zlb,'Position',pos_zlb)

    function maxlen = label_max_len(ax)
        ticklabel = get(h_a, [ax, 'TickLabel']);
        if vers < 8.4
            maxlen = size(ticklabel,2);
        else
            maxlen = 0;
            for k_f = 1:length(ticklabel)
                if length(ticklabel{k_f}) > maxlen
                    maxlen = length(ticklabel{k_f});
                end
            end
        end
    end
end
function axislabel_translation_slider
% Create UI figure window and components

global AXISALIGN_TRANS_A AXISALIGN_TRANS_B
if isempty(AXISALIGN_TRANS_A), AXISALIGN_TRANS_A = 1; end
if isempty(AXISALIGN_TRANS_B), AXISALIGN_TRANS_B = 1; end

h_f = uifigure('Position',[100 100 350 275]);

% cg = uigauge(fig,'Position',[100 100 120 120]);

sld1 = uislider(h_f,...
    'Position',[50 200 250 3],...
    'Limits',[0 5],...
    'ValueChangingFcn',@(sld,event) sliderMoving1(event));
sld1.Value = 1;
sld2 = uislider(h_f,...
    'Position',[50 75 250 3],...
    'Limits',[0 5],...
    'ValueChangingFcn',@(sld,event) sliderMoving2(event));
sld2.Value = 1;
end

% Create ValueChangingFcn callback
function sliderMoving1(event)
global AXISALIGN_TRANS_A
AXISALIGN_TRANS_A = event.Value;
align_axislabel([],gca);
end

function sliderMoving2(event)
global AXISALIGN_TRANS_B
AXISALIGN_TRANS_B = event.Value;
align_axislabel([],gca);
end
