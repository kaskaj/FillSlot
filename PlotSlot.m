function [] = PlotSlot(x,y,R,c)

figure;

%% Slot

plot(x,y,'k');
hold on; axis equal; grid on;

%% Conductors

c_x = c(:,1) + R*cos(0:0.2:2*pi);
c_y = c(:,2) + R*sin(0:0.2:2*pi);
patch('XData',c_x','YData',c_y','FaceColor','None');
hold off;

end