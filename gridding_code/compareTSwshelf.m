function [outputArg1,outputArg2] = compareTS(grid,name,shelf)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


objmapT = grid.T.objmap;
errmapT = grid.T.errmap;
residualsT = grid.T.residuals;
objmapS = grid.S.objmap;
errmapS = grid.S.errmap;
residualsS = grid.S.residuals;
xg = grid.T.xg;
yg = grid.T.yg;
dataT = grid.T.data;
dataS = grid.S.data;
xdata = grid.T.xdata;
ydata = grid.T.ydata;

shelfminCT = shelf.CTmin_z;
shelfminSA = shelf.SAmin_z;

%% Variables needed for plotting TS plots
%add Gade line. This requires defining the Atlantic water mass...
% T_SMW = -87.4;
% AW_CT = AW5.CT(9); %2019 calculed AW average properties
% AW_SA = AW5.SA(9);
%%%%%%Find line slope and create data
% m = (-1*T_SMW+AW_CT)/(AW_SA);
% xmelt = linspace(32.7,AW_SA,100);
% ymelt = m*xmelt+T_SMW;

%%%%A more custom T-S plot
%First we need to range of temperature and salinity values that we expect
%our TS plot to have.
sx = [20:.01:36]; % range of absolute salinity values 
ty = [-2:.01:6]; % range of conservative temp values

[S,T] = meshgrid(sx,ty); %creates a grid we can calculate density on

%%%%Make freezing line
CT_freezing = gsw_CT_freezing(sx,zeros(size(sx)));

%%%Add runoff line
%%%%%%Find line slope and create data
% xrunoff = linspace(27,AW_SA,100);
% yrunoff = ones(size(xrunoff)).*AW_CT;

%Now that we have grids of salinity and temperature, we can make a grid of
%potential density anomaly from that
PDEN = gsw_rho(S,T,0)-1000;

%Now we want just the PDEN grid for only where GMW is to shade in on the
%plot
[ii, jj] = size(PDEN);


% figure;
% set(gcf,'position',[10,10,800,800])
% %Plot all contour lines
% [C,h] = contour(S,T,PDEN,[20:29],'color',[.58 .58 .65],'ShowText','on'); % contour plot of density
% hold on
% clabel(C,h,'FontSize',12,'color',[.58 .58 .65])
% %plot mixing lines
% %h1 = plot(xmelt,ymelt,'k--','LineWidth',1);
% h2 = plot(sx,CT_freezing,'r--','LineWidth',1);
% %h3 = plot(xrunoff,yrunoff,'k--','LineWidth',1);
% xlim([28 36])
% ylim([-2 5.5])
% xlabel('S_{A} [g/kg]')
% ylabel('\Theta (^oC)')
% %title('SMW %')
% 
% h4 = scatter(objmapS(:),objmapT(:),20,'filled','MarkerFaceColor',[.58 .58 .65],'DisplayName','ObjMap');
% h5 = scatter(dataS(:),dataT(:),20,'filled','b','DisplayName','Data');
% 
% set(gca,'FontSize',30)
% 
% file_name = strcat(name,'_TScomp','.png');
% 
% % Save the plot as a PNG file in the specified folder
% saveas(gcf, file_name);

%Plot as profile lines instead of scatter
figure;
set(gcf,'position',[10,10,800,800])
%Plot all contour lines
[C,h] = contour(S,T,PDEN,[20:29],'color',[.58 .58 .65],'ShowText','on'); % contour plot of density
hold on
clabel(C,h,'FontSize',12,'color',[.58 .58 .65])
%plot mixing lines
%h1 = plot(xmelt,ymelt,'k--','LineWidth',1);
h2 = plot(sx,CT_freezing,'k--','LineWidth',1);
%h3 = plot(xrunoff,yrunoff,'k--','LineWidth',1);
xlim([28 36])
ylim([-2 5.5])
xlabel('S_{A} [g/kg]')
ylabel('\Theta (^oC)')
%title('SMW %')
[ll,kk] = size(objmapT);
[jj,mm] = size(dataT);

h3 = plot(shelfminSA, shelfminCT,'Color','r','LineWidth',1.5)

for i = 1:kk
    plot(objmapS(:,i),objmapT(:,i),'Color',[.58 .58 .65],'LineWidth',1.5);
end
for j = 1:mm
    plot(dataS(:,j),dataT(:,j),'Color','b','LineWidth',1.5);
end

set(gca,'FontSize',30)

file_name = strcat(name,'_TScomp3_wshelf','.png');

% Save the plot as a PNG file in the specified folder
saveas(gcf, file_name);

end

