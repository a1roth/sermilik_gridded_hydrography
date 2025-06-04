function [outputArg1,outputArg2] = makeSectionPlots(grid,xbathy,ybathy,name)
%UNTITLED2 Summary of this function goes here
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

rawT = grid.T.data_raw;
rawS = grid.S.data_raw;

%Get x and y data coordinate pairs for only where there is data of input
%map data
ii=find(~isnan(residualsT));
x = xdata(ii);
y=ydata(ii);

%Get x and y data coordinate pairs for only where there is data from raw
%data
iiraw=find(~isnan(rawT));
xraw = xdata(iiraw);
yraw=ydata(iiraw);
Traw = rawT(iiraw);

%Calculate density for isopycnals (requires Gibbs Seawater Toolbox TEOS-10)

sigma0 = gsw_sigma0(objmapS,objmapT);
%smooth ispycnals
sigma0  = smoothdata(sigma0,2,'movmean',5,'omitnan');
sigma0  = smoothdata(sigma0,2,'movmean',10,'omitnan');
sigma0  = smoothdata(sigma0,1,'movmean',4,'omitnan');

%List of isopycnals to show and label on plot
isopyc_lev = [26,27,27.4,27.5]; 


%% Make a nice along fjord section
%Start with TEMP
%figure1 = figure;
f = figure;
%sp1 = subplot(1,3,1,'Parent',f);
%axes1 = axes('Parent',figure1);
set(gcf,'position',[10,10,950,370])
hold on;
%hold(axes1,'on');

%h1 = contourf(xg,yg.*1000,objmap',30,'linecolor','none');
h1 = pcolor(xg,yg.*1000,objmapT);
shading flat


%Create colormap from cbrewer
colors = cbrewer('div','RdYlBu',150,'cubic');
colors(find(colors>=1))=1;
colors(find(colors<=0))=0;
colors = flipud(colors(3:end,:)); %Gets rid of dark dark red that I don't like

%Set colormap
%colormap(axes1,colors);
colormap(colors)
%Draw colorbar
axes1 = gca;
bar = colorbar(axes1,'FontSize',16);
drawnow;

%k1 = get(gca,'position');
%k=get(bar, 'Position');
%k(3) = k(3)*1.3;
%set(bar, 'Position',k);
%set(gca,'position',k1);

%Label colorbar and adjust position
ylabel(bar,'\Theta (°C)','FontSize',16,'Rotation',360);
bar.Label.Position(1) = 1.9;
bar.Label.Position(2) = 5.4;

%Plot locations of the CTD stations
h2 = scatter(xdata(1,:),20*ones(1,length(xdata(1,:))),65,'v','MarkerEdgeColor','none','MarkerFaceColor',[128/255 133/255 133/255]);
%h3 = scatter(x(:),y(:).*1000,'.','MarkerEdgeColor',[88/255 93/255 93/255]);
h3 = scatter(x(:),y(:).*1000,2,'filled','MarkerFaceColor',[218/255 255/255 255/255]);%plotting bottom filled profiles in white 
h4 = scatter(xraw(:),yraw(:).*1000,20,Traw(:),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',0.3);%plotting raw profiles - not bottom filled - in grey


%h3 = scatter(xx,depth.*1000,10,data+bgrnd_ctd,'filled');

%%Create text labels for CTD stations. Only apply text to 2 - 19 and then
%32. Other stations are too dense for labels
%a = [5,1,1,1,6,1,2,2,2]'; 
%b= num2str(a);
%c = cellstr(b);
%dx = -2.4; dy = 30;
%T = text(latpoints,20.*ones(1,length(a))+dy,c,'FontSize',12);
%T.Color = [128/255 133/255 133/255];


%Add title to plot and adjust position so that it sits above CTD station
%text labels
%title(strcat('2019 - Objective Map - ','xcora= ', num2str(xcora), ' ycora= ', num2str(ycora),' error= ', num2str(err)),'FontSize',16)
set(gca, 'Units', 'normalized')
titleHandle = get(gca, 'Title');
pos = get(titleHandle, 'position');
newpos = pos + [0 35 0];
set(titleHandle,'position', newpos);

%List of isopycnals to show and label on plot
%isopyc_lev = [26,27,27.4,27.5,27.6,27.7]; 
%Plot isopycnals overtop temperature field
h4 = contour(xg,yg.*1000,sigma0 ,isopyc_lev,'linecolor','k','ShowText','on');

%Plot bathymetry over everything. 
basevalue = -1000;
h5 = area(xbathy,ybathy,basevalue);
h5(1).FaceColor = [0.7 0.7 0.7]; %Creates light grey bathymetry
h5(1).EdgeColor = [0.7 0.7 0.7];
h5(1).LineWidth = 1;

T1 = text(2,-980,'Glacier','FontSize',16);
T1.Color = [128/255 133/255 133/255];
set(T1,'Rotation',90);
T2 = text(108,-980,'Shelf','FontSize',16);
set(T2,'Rotation',90);
T2.Color = [128/255 133/255 133/255];


ylim([-1000,25])
xlim([0,110])
xlabel('Distance from glacier (km)', 'FontSize',18)
ylabel('Depth (m)','FontSize',18)
axes1 = gca;
box(axes1,'off');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','CLim',[-1 5],'FontSize',14,'Layer','top');


hold off

file_name = strcat(name,'_Tsection_wRaw','.png');

% Save the plot as a PNG file in the specified folder
saveas(gcf, file_name);


%% Make a nice along fjord section with raw data on top 
%Start with TEMP
%figure1 = figure;
f = figure;
%sp1 = subplot(1,3,1,'Parent',f);
%axes1 = axes('Parent',figure1);
set(gcf,'position',[10,10,950,370])
hold on;
%hold(axes1,'on');

%h1 = contourf(xg,yg.*1000,objmap',30,'linecolor','none');
h1 = pcolor(xg,yg.*1000,objmapT);
shading flat


%Create colormap from cbrewer
colors = cbrewer('div','RdYlBu',150,'cubic');
colors(find(colors>=1))=1;
colors(find(colors<=0))=0;
colors = flipud(colors(3:end,:)); %Gets rid of dark dark red that I don't like

%Set colormap
%colormap(axes1,colors);
colormap(colors)
%Draw colorbar
axes1 = gca;
bar = colorbar(axes1,'FontSize',16);
drawnow;

%k1 = get(gca,'position');
%k=get(bar, 'Position');
%k(3) = k(3)*1.3;
%set(bar, 'Position',k);
%set(gca,'position',k1);

%Label colorbar and adjust position
ylabel(bar,'\Theta (°C)','FontSize',16,'Rotation',360);
bar.Label.Position(1) = 1.9;
bar.Label.Position(2) = 5.4;

%Plot locations of the CTD stations
h2 = scatter(xdata(1,:),20*ones(1,length(xdata(1,:))),65,'v','MarkerEdgeColor','none','MarkerFaceColor',[128/255 133/255 133/255]);
%h3 = scatter(x(:),y(:).*1000,'.','MarkerEdgeColor',[88/255 93/255 93/255]);
h3 = scatter(x(:),y(:).*1000,2,'filled','MarkerFaceColor',[218/255 255/255 255/255]);%plotting bottom filled profiles in white 
%h4 = scatter(xraw(:),yraw(:).*1000,2,'filled','MarkerFaceColor',[88/255 93/255 93/255]);%plotting raw profiles - not bottom filled - in grey
h4 = scatter(xraw(:),yraw(:).*1000,2,'filled','MarkerFaceColor',[88/255 93/255 93/255]);%plotting raw profiles - not bottom filled - in grey



%h3 = scatter(xx,depth.*1000,10,data+bgrnd_ctd,'filled');

%%Create text labels for CTD stations. Only apply text to 2 - 19 and then
%32. Other stations are too dense for labels
%a = [5,1,1,1,6,1,2,2,2]'; 
%b= num2str(a);
%c = cellstr(b);
%dx = -2.4; dy = 30;
%T = text(latpoints,20.*ones(1,length(a))+dy,c,'FontSize',12);
%T.Color = [128/255 133/255 133/255];


%Add title to plot and adjust position so that it sits above CTD station
%text labels
%title(strcat('2019 - Objective Map - ','xcora= ', num2str(xcora), ' ycora= ', num2str(ycora),' error= ', num2str(err)),'FontSize',16)
set(gca, 'Units', 'normalized')
titleHandle = get(gca, 'Title');
pos = get(titleHandle, 'position');
newpos = pos + [0 35 0];
set(titleHandle,'position', newpos);

%List of isopycnals to show and label on plot
%isopyc_lev = [26,27,27.4,27.5,27.6,27.7]; 
%Plot isopycnals overtop temperature field
h4 = contour(xg,yg.*1000,sigma0 ,isopyc_lev,'linecolor','k','ShowText','on');

%Plot bathymetry over everything. 
basevalue = -1000;
h5 = area(xbathy,ybathy,basevalue);
h5(1).FaceColor = [0.7 0.7 0.7]; %Creates light grey bathymetry
h5(1).EdgeColor = [0.7 0.7 0.7];
h5(1).LineWidth = 1;

T1 = text(2,-980,'Glacier','FontSize',16);
T1.Color = [128/255 133/255 133/255];
set(T1,'Rotation',90);
T2 = text(108,-980,'Shelf','FontSize',16);
set(T2,'Rotation',90);
T2.Color = [128/255 133/255 133/255];


ylim([-1000,25])
xlim([0,110])
xlabel('Distance from glacier (km)', 'FontSize',18)
ylabel('Depth (m)','FontSize',18)
axes1 = gca;
box(axes1,'off');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','CLim',[-1 5],'FontSize',14,'Layer','top');


hold off

file_name = strcat(name,'_Tsection','.png');

% Save the plot as a PNG file in the specified folder
saveas(gcf, file_name);






%% Plot Error

f = figure;
%sp1 = subplot(1,3,1,'Parent',f);
%axes1 = axes('Parent',figure1);
set(gcf,'position',[10,10,950,370])
hold on;
%hold(axes1,'on');

%h1 = contourf(xg,yg.*1000,objmap',30,'linecolor','none');
h1 = pcolor(xg,yg.*1000,errmapT);
shading flat
colorbar



%Plot locations of the CTD stations
h2 = scatter(xdata(1,:),20*ones(1,length(xdata(1,:))),65,'v','MarkerEdgeColor','none','MarkerFaceColor',[128/255 133/255 133/255]);
h3 = scatter(x(:),y(:).*1000,1,'filled','MarkerFaceColor','k');%[88/255 93/255 93/255]);
h3 = scatter(x(:),y(:).*1000,2,'filled','MarkerFaceColor',[218/255 255/255 255/255]);%plotting bottom filled profiles in white 
h4 = scatter(xraw(:),yraw(:).*1000,2,'filled','MarkerFaceColor',[88/255 93/255 93/255]);%plotting raw profiles - not bottom filled - in grey


%%Create text labels for CTD stations. Only apply text to 2 - 19 and then
%32. Other stations are too dense for labels
%a = [5,1,1,1,6,1,2,2,2]'; 
%b= num2str(a);
%c = cellstr(b);
%dx = -2.4; dy = 30;
%T = text(latpoints,20.*ones(1,length(a))+dy,c,'FontSize',12);
%T.Color = [128/255 133/255 133/255];


%Add title to plot and adjust position so that it sits above CTD station
%text labels
%title(strcat('2019 - Uncertainty - ','xcora= ', num2str(xcora), ' ycora= ', num2str(ycora),' error= ', num2str(err)),'FontSize',16)
set(gca, 'Units', 'normalized')
titleHandle = get(gca, 'Title');
pos = get(titleHandle, 'position');
newpos = pos + [0 35 0];
set(titleHandle,'position', newpos);

%List of isopycnals to show and label on plot
%isopyc_lev = [26,27,27.4,27.55,27.7,27.8]; 
%Plot isopycnals overtop temperature field
% h4 = contour(xg,yg.*1000,sigma0 ,isopyc_lev,'linecolor','k','ShowText','on');

%Plot bathymetry over everything. 
basevalue = -1000;
h5 = area(xbathy,ybathy,basevalue);
h5(1).FaceColor = [0.7 0.7 0.7]; %Creates light grey bathymetry
h5(1).EdgeColor = [0.7 0.7 0.7];
h5(1).LineWidth = 1;


T1 = text(2,-980,'Glacier','FontSize',16);
T1.Color = [128/255 133/255 133/255];
set(T1,'Rotation',90);
T2 = text(108,-980,'Shelf','FontSize',16);
set(T2,'Rotation',90);
T2.Color = [128/255 133/255 133/255];


ylim([-1000,25])
xlim([0,110])
xlabel('Distance from glacier (km)', 'FontSize',18)
ylabel('Depth (m)','FontSize',18)
axes1 = gca;
box(axes1,'off');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','CLim',[0 0.35],'FontSize',14,'Layer','top');

hold off

file_name = strcat(name,'_Tsection_error','.png');

% Save the plot as a PNG file in the specified folder
saveas(gcf, file_name);


%% Plot residuals

f = figure;
%sp1 = subplot(1,3,1,'Parent',f);
%axes1 = axes('Parent',figure1);
set(gcf,'position',[10,10,950,370])
hold on;
%hold(axes1,'on');

%h1 = contourf(xg,yg.*1000,objmap',30,'linecolor','none');
%h1 = pcolor(xg,yg.*1000,errmap');
h1 = scatter(xdata(:),ydata(:).*1000,20,residualsT(:),'filled');
%shading flat



%Create colormap from cbrewer
colors = cbrewer('div','PRGn',50,'cubic');
colors(find(colors>=1))=1;
colors(find(colors<=0))=0;
colors = flipud(colors(1:end,:)); %Gets rid of dark dark red that I don't like

%Set colormap
%colormap(axes1,colors);
colormap(colors)
%Draw colorbar
axes1 = gca;
bar = colorbar(axes1,'FontSize',16);
drawnow;

%k1 = get(gca,'position');
%k=get(bar, 'Position');
%k(3) = k(3)*1.3;
%set(bar, 'Position',k);
%set(gca,'position',k1);

%Label colorbar and adjust position
ylabel(bar,'\Theta (°C)','FontSize',16,'Rotation',360);
bar.Label.Position(1) = 1.9;
bar.Label.Position(2) = 5.4;

%Plot locations of the CTD stations


%Plot locations of the CTD stations
h2 = scatter(xdata(1,:),20*ones(1,length(xdata(1,:))),65,'v','MarkerEdgeColor','none','MarkerFaceColor',[128/255 133/255 133/255]);
%h3 = scatter(xdata(:),ydata(:).*1000,'.','MarkerEdgeColor',[88/255 93/255 93/255]);

%%Create text labels for CTD stations. Only apply text to 2 - 19 and then
%32. Other stations are too dense for labels
%a = [5,1,1,1,6,1,2,2,2]'; 
%b= num2str(a);
%c = cellstr(b);
%dx = -2.4; dy = 30;
%T = text(latpoints,20.*ones(1,length(a))+dy,c,'FontSize',12);
%T.Color = [128/255 133/255 133/255];


%Add title to plot and adjust position so that it sits above CTD station
%text labels
%title(strcat('2019 - Uncertainty - ','xcora= ', num2str(xcora), ' ycora= ', num2str(ycora),' error= ', num2str(err)),'FontSize',16)
set(gca, 'Units', 'normalized')
titleHandle = get(gca, 'Title');
pos = get(titleHandle, 'position');
newpos = pos + [0 35 0];
set(titleHandle,'position', newpos);

%List of isopycnals to show and label on plot
%isopyc_lev = [26,27,27.4,27.55,27.7,27.8]; 
%Plot isopycnals overtop temperature field
%h4 = contour(xg,yg.*1000,objmap_rho' ,isopyc_lev,'linecolor','k','ShowText','on');

%Plot bathymetry over everything. 
basevalue = -1000;
h5 = area(xbathy,ybathy,basevalue);
h5(1).FaceColor = [0.7 0.7 0.7]; %Creates light grey bathymetry
h5(1).EdgeColor = [0.7 0.7 0.7];
h5(1).LineWidth = 1;


T1 = text(2,-980,'Glacier','FontSize',16);
T1.Color = [128/255 133/255 133/255];
set(T1,'Rotation',90);
T2 = text(108,-980,'Shelf','FontSize',16);
set(T2,'Rotation',90);
T2.Color = [128/255 133/255 133/255];


ylim([-1000,25])
xlim([0,110])
xlabel('Distance from glacier (km)', 'FontSize',18)
ylabel('Depth (m)','FontSize',18)
axes1 = gca;
box(axes1,'off');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','CLim',[-0.6 0.6],'FontSize',14,'Layer','top');

hold off

file_name = strcat(name,'_Tsection_residuals','.png');

% Save the plot as a PNG file in the specified folder
saveas(gcf, file_name);

%% Make a nice along fjord section
%SALINIRY
%figure1 = figure;
f = figure;
%sp1 = subplot(1,3,1,'Parent',f);
%axes1 = axes('Parent',figure1);
set(gcf,'position',[10,10,950,370])
hold on;
%hold(axes1,'on');

%h1 = contourf(xg,yg.*1000,objmap',30,'linecolor','none');
h1 = pcolor(xg,yg.*1000,objmapS);
shading flat


% %Create colormap from cbrewer
% colors = cbrewer('seq','YlGnBu',150,'cubic');
% colors(find(colors>=1))=1;
% colors(find(colors<=0))=0;
% colors = flipud(colors(1:135,:)); %Gets rid of dark dark red that I don't like
% 
% %Set colormap
% %colormap(axes1,colors);
% colormap(colors)
% %cmap = cmocean('haline',150);
% %colormap(cmap)
% %Draw colorbar
% axes1 = gca;
% bar = colorbar(axes1,'FontSize',16);
% drawnow;


% %Create colormap from cbrewer
colors = cbrewer('div','Spectral',170,'cubic');
colors(find(colors>=1))=1;
colors(find(colors<=0))=0;
colors = flipud(colors(50:end,:)); %Gets rid of dark dark red that I don't like

%Set colormap
%colormap(axes1,colors);
colormap(colors)
%cmap = cmocean('haline',150);
%colormap(cmap)
%Draw colorbar
axes1 = gca;
bar = colorbar(axes1,'FontSize',16);
drawnow;



%Create colormap from cbrewer
% colors = cbrewer('div','Spectral',28,'cubic');
% colors(find(colors>=1))=1;
% colors(find(colors<=0))=0;
% colors = flipud(colors(6:end,:));
% 
% 
% %Adjust colormap to make uneven intervals - use for salinity
% cmap = [ repmat(colors(1,:), 40, 1) %%32.50 - 32.90
%          repmat(colors(2,:), 45, 1) %%32.90 - 33.30
%          repmat(colors(3,:), 50, 1) %%33.30 - 33.80
%          repmat(colors(4,:), 45, 1) %%33.80 - 34.20
%          repmat(colors(5,:), 30, 1) %%34.20 - 34.50
%          repmat(colors(6,:), 10, 1) %%34.50 - 34.60
%          repmat(colors(7,:), 10, 1) %%34.60 - 34.70
%          repmat(colors(8,:), 5, 1) %%34.70 - 34.75
%          repmat(colors(9,:), 5, 1) %%34.75 - 34.80
%          repmat(colors(10,:), 5, 1) %%34.80 - 34.85
%          repmat(colors(11,:), 5, 1) %%34.85 - 34.90
%          repmat(colors(12,:), 5, 1) %%34.90 - 34.95
%          repmat(colors(13,:), 5, 1) %%34.95 - 35.00
%          repmat(colors(14,:), 2, 1) %%35.00 - 35.02
%          repmat(colors(15,:), 3, 1) %%35.02 - 35.05
%          repmat(colors(16,:), 2, 1) %%35.05 - 35.07
%          repmat(colors(17,:), 3, 1) %%35.07 - 35.10
%          repmat(colors(18,:), 5, 1) %%35.10 - 35.15
%          repmat(colors(19,:), 5, 1) %%35.15 - 35.20
%          repmat(colors(20,:), 5, 1) %%35.20 - 35.25
%          repmat(colors(21,:), 5, 1) %%35.25 - 35.3
%          repmat(colors(22,:), 5, 1) %%>35.3      
%          ]; %
%      
% colormap(cmap)
% 
% %Draw colorbar with only certain ticks labeled
% bar = colorbar;%(axes1,'FontSize',16);%,'Ticks',[32.50,32.90,33.30,33.80,34.20,34.50,34.60,34.70,34.80,34.90,35.00,35.10,35.20,35.3]);%,'TickLabels',{32.50,32.90,33.30,33.80,34.20,34.50,34.60,34.70,34.80,34.90,35.00,35.10,35.20,35.3});
% bar.Ticks = [32.50,32.90,33.30,33.80,34.20,34.50,34.60,34.70,34.80,34.90,35.00,35.10,35.20,35.3];
% bar.TickLabels = [32.50,32.90,33.30,33.80,34.20,34.50,34.60,34.70,34.80,34.90,35.00,35.10,35.20,35.3];
% drawnow;


%k1 = get(gca,'position');
%k=get(bar, 'Position');
%k(3) = k(3)*1.3;
%set(bar, 'Position',k);
%set(gca,'position',k1);

%Label colorbar and adjust position
ylabel(bar,'S_{A} (g/kg)','FontSize',16,'Rotation',360);
bar.Label.Position(1) = 1.9;
bar.Label.Position(2) = 5.4;

%Plot locations of the CTD stations
h2 = scatter(xdata(1,:),20*ones(1,length(xdata(1,:))),65,'v','MarkerEdgeColor','none','MarkerFaceColor',[128/255 133/255 133/255]);
h3 = scatter(x(:),y(:).*1000,2,'filled','MarkerFaceColor',[218/255 255/255 255/255]);%plotting bottom filled profiles in white 
h4 = scatter(xraw(:),yraw(:).*1000,2,'filled','MarkerFaceColor',[88/255 93/255 93/255]);%plotting raw profiles - not bottom filled - in grey


%%Create text labels for CTD stations. Only apply text to 2 - 19 and then
%32. Other stations are too dense for labels
%a = [5,1,1,1,6,1,2,2,2]'; 
%b= num2str(a);
%c = cellstr(b);
%dx = -2.4; dy = 30;
%T = text(latpoints,20.*ones(1,length(a))+dy,c,'FontSize',12);
%T.Color = [128/255 133/255 133/255];


%Add title to plot and adjust position so that it sits above CTD station
%text labels
%title(strcat('2019 - Objective Map - ','xcora= ', num2str(xcora), ' ycora= ', num2str(ycora),' error= ', num2str(err)),'FontSize',16)
set(gca, 'Units', 'normalized')
titleHandle = get(gca, 'Title');
pos = get(titleHandle, 'position');
newpos = pos + [0 35 0];
set(titleHandle,'position', newpos);

%List of isopycnals to show and label on plot
%isopyc_lev = [26,27,27.4,27.5,27.6,27.7]; 
%Plot isopycnals overtop temperature field
h4 = contour(xg,yg.*1000,sigma0 ,isopyc_lev,'linecolor','k','ShowText','on');

%Plot bathymetry over everything. 
basevalue = -1000;
h5 = area(xbathy,ybathy,basevalue);
h5(1).FaceColor = [0.7 0.7 0.7]; %Creates light grey bathymetry
h5(1).EdgeColor = [0.7 0.7 0.7];
h5(1).LineWidth = 1;

T1 = text(2,-980,'Glacier','FontSize',16);
T1.Color = [128/255 133/255 133/255];
set(T1,'Rotation',90);
T2 = text(108,-980,'Shelf','FontSize',16);
set(T2,'Rotation',90);
T2.Color = [128/255 133/255 133/255];


ylim([-1000,25])
xlim([0,110])
xlabel('Distance from glacier (km)', 'FontSize',18)
ylabel('Depth (m)','FontSize',18)
axes1 = gca;
box(axes1,'off');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','CLim',[27 35],'FontSize',14,'Layer','top');


hold off

file_name = strcat(name,'_Ssection','.png');

% Save the plot as a PNG file in the specified folder
saveas(gcf, file_name);





%% Plot Error

f = figure;
%sp1 = subplot(1,3,1,'Parent',f);
%axes1 = axes('Parent',figure1);
set(gcf,'position',[10,10,950,370])
hold on;
%hold(axes1,'on');

%h1 = contourf(xg,yg.*1000,objmap',30,'linecolor','none');
h1 = pcolor(xg,yg.*1000,errmapS);
shading flat
colorbar



%Plot locations of the CTD stations
h2 = scatter(xdata(1,:),20*ones(1,length(xdata(1,:))),65,'v','MarkerEdgeColor','none','MarkerFaceColor',[128/255 133/255 133/255]);
h3 = scatter(x(:),y(:).*1000,1,'filled','MarkerFaceColor','k');%[88/255 93/255 93/255]);
h3 = scatter(x(:),y(:).*1000,2,'filled','MarkerFaceColor',[218/255 255/255 255/255]);%plotting bottom filled profiles in white 
h4 = scatter(xraw(:),yraw(:).*1000,2,'filled','MarkerFaceColor',[88/255 93/255 93/255]);%plotting raw profiles - not bottom filled - in grey


%%Create text labels for CTD stations. Only apply text to 2 - 19 and then
%32. Other stations are too dense for labels
%a = [5,1,1,1,6,1,2,2,2]'; 
%b= num2str(a);
%c = cellstr(b);
%dx = -2.4; dy = 30;
%T = text(latpoints,20.*ones(1,length(a))+dy,c,'FontSize',12);
%T.Color = [128/255 133/255 133/255];


%Add title to plot and adjust position so that it sits above CTD station
%text labels
%title(strcat('2019 - Uncertainty - ','xcora= ', num2str(xcora), ' ycora= ', num2str(ycora),' error= ', num2str(err)),'FontSize',16)
set(gca, 'Units', 'normalized')
titleHandle = get(gca, 'Title');
pos = get(titleHandle, 'position');
newpos = pos + [0 35 0];
set(titleHandle,'position', newpos);

%List of isopycnals to show and label on plot
%isopyc_lev = [26,27,27.4,27.55,27.7,27.8]; 
%Plot isopycnals overtop temperature field
% h4 = contour(xg,yg.*1000,sigma0 ,isopyc_lev,'linecolor','k','ShowText','on');

%Plot bathymetry over everything. 
basevalue = -1000;
h5 = area(xbathy,ybathy,basevalue);
h5(1).FaceColor = [0.7 0.7 0.7]; %Creates light grey bathymetry
h5(1).EdgeColor = [0.7 0.7 0.7];
h5(1).LineWidth = 1;


T1 = text(2,-980,'Glacier','FontSize',16);
T1.Color = [128/255 133/255 133/255];
set(T1,'Rotation',90);
T2 = text(108,-980,'Shelf','FontSize',16);
set(T2,'Rotation',90);
T2.Color = [128/255 133/255 133/255];


ylim([-1000,25])
xlim([0,110])
xlabel('Distance from glacier (km)', 'FontSize',18)
ylabel('Depth (m)','FontSize',18)
axes1 = gca;
box(axes1,'off');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','CLim',[0 0.35],'FontSize',14,'Layer','top');

hold off

file_name = strcat(name,'_Ssection_error','.png');

% Save the plot as a PNG file in the specified folder
saveas(gcf, file_name);


%% Plot residuals

f = figure;
%sp1 = subplot(1,3,1,'Parent',f);
%axes1 = axes('Parent',figure1);
set(gcf,'position',[10,10,950,370])
hold on;
%hold(axes1,'on');

%h1 = contourf(xg,yg.*1000,objmap',30,'linecolor','none');
%h1 = pcolor(xg,yg.*1000,errmap');
h1 = scatter(xdata(:),ydata(:).*1000,20,residualsS(:),'filled');
%shading flat



%Create colormap from cbrewer
colors = cbrewer('div','PRGn',50,'cubic');
colors(find(colors>=1))=1;
colors(find(colors<=0))=0;
colors = flipud(colors(1:end,:)); %Gets rid of dark dark red that I don't like

%Set colormap
%colormap(axes1,colors);
colormap(colors)
%Draw colorbar
axes1 = gca;
bar = colorbar(axes1,'FontSize',16);
drawnow;

%k1 = get(gca,'position');
%k=get(bar, 'Position');
%k(3) = k(3)*1.3;
%set(bar, 'Position',k);
%set(gca,'position',k1);

%Label colorbar and adjust position
ylabel(bar,'\Theta (°C)','FontSize',16,'Rotation',360);
bar.Label.Position(1) = 1.9;
bar.Label.Position(2) = 5.4;

%Plot locations of the CTD stations


%Plot locations of the CTD stations
h2 = scatter(xdata(1,:),20*ones(1,length(xdata(1,:))),65,'v','MarkerEdgeColor','none','MarkerFaceColor',[128/255 133/255 133/255]);
%h3 = scatter(xdata(:),ydata(:).*1000,'.','MarkerEdgeColor',[88/255 93/255 93/255]);

%%Create text labels for CTD stations. Only apply text to 2 - 19 and then
%32. Other stations are too dense for labels
%a = [5,1,1,1,6,1,2,2,2]'; 
%b= num2str(a);
%c = cellstr(b);
%dx = -2.4; dy = 30;
%T = text(latpoints,20.*ones(1,length(a))+dy,c,'FontSize',12);
%T.Color = [128/255 133/255 133/255];


%Add title to plot and adjust position so that it sits above CTD station
%text labels
%title(strcat('2019 - Uncertainty - ','xcora= ', num2str(xcora), ' ycora= ', num2str(ycora),' error= ', num2str(err)),'FontSize',16)
set(gca, 'Units', 'normalized')
titleHandle = get(gca, 'Title');
pos = get(titleHandle, 'position');
newpos = pos + [0 35 0];
set(titleHandle,'position', newpos);

%List of isopycnals to show and label on plot
%isopyc_lev = [26,27,27.4,27.55,27.7,27.8]; 
%Plot isopycnals overtop temperature field
%h4 = contour(xg,yg.*1000,objmap_rho' ,isopyc_lev,'linecolor','k','ShowText','on');

%Plot bathymetry over everything. 
basevalue = -1000;
h5 = area(xbathy,ybathy,basevalue);
h5(1).FaceColor = [0.7 0.7 0.7]; %Creates light grey bathymetry
h5(1).EdgeColor = [0.7 0.7 0.7];
h5(1).LineWidth = 1;


T1 = text(2,-980,'Glacier','FontSize',16);
T1.Color = [128/255 133/255 133/255];
set(T1,'Rotation',90);
T2 = text(108,-980,'Shelf','FontSize',16);
set(T2,'Rotation',90);
T2.Color = [128/255 133/255 133/255];


ylim([-1000,25])
xlim([0,110])
xlabel('Distance from glacier (km)', 'FontSize',18)
ylabel('Depth (m)','FontSize',18)
axes1 = gca;
box(axes1,'off');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'BoxStyle','full','CLim',[-0.6 0.6],'FontSize',14,'Layer','top');

hold off

file_name = strcat(name,'_Ssection_residuals','.png');

% Save the plot as a PNG file in the specified folder
saveas(gcf, file_name);





end

