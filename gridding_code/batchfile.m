%batchfile
%Run this in the folder you would like all your results in. 
%Use this script to grid multiple profile sets in one go and plot
%everything up too! 
%This is where you can specify parameters. 

%List parameters for this run
xcor1 = 50;
ycor1 = 0.100;
a1 = 0.6;
xcor2 = 15;
ycor2 = 0.010;
a2 = 0.4;
noise = 0.5;
xres = 2; %km
yres = 5; %m
xstart = 0;%km
xend = 104;%km
ystart = 0;%m
yend = 900;%m
notes = 'This version uses a slightly different way of projecting the original profile coordinates onto the thalweg section lat/lon. This version also interpolates to a standard grid';


%Or load netcdf file if that's better
%background_CT = ncread('backgrounds.nc','background_CT')
%background_SA = ncread('backgrounds.nc','background_SA')
%xback = ncread('backgrounds.nc','dist')
%yback = ncread('backgrounds.nc','depth')

%Create folder for everything to be saved in
foldername = 'ObjMap_allyears_2km5mres_paramset_1';
mkdir(foldername)
%Add folder to path
addpath(foldername)
%Move into that folder
cd(foldername)



%% RUN ALL YEARS

objmap_2009(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2010(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2011(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2012(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2013(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2015(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2016(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2017(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)
 
objmap_2018(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2019(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2021(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2022(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2023(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

objmap_2023xctd(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)



%%

names = ['2009','2010','2011','2012','2013','2015','2016','2017','2018','2019','2021','2022','2023','2023xctd'];

%load thalweg section
dist_bed = ncread('SermilikThalwegBathy.nc','dist'); %distance along thalweg section from Helheim 2019 terminus position
latZbed = ncread('SermilikThalwegBathy.nc','lat'); %latitude of points along thalweg section
lonZbed = ncread('SermilikThalwegBathy.nc','lon'); %longitude of points along thalweg section
Zsmooth = ncread('SermilikThalwegBathy.nc','depth'); %smoothed bathymetry depth values at Lat/Lon locations

makeSectionPlots(grids2009,dist_bed,Zsmooth,'2009')
makeSectionPlots(grids2010,dist_bed,Zsmooth,'2010')
makeSectionPlots(grids2011,dist_bed,Zsmooth,'2011')
makeSectionPlots(grids2012,dist_bed,Zsmooth,'2012')
makeSectionPlots(grids2013,dist_bed,Zsmooth,'2013')
makeSectionPlots(grids2015,dist_bed,Zsmooth,'2015')
makeSectionPlots(grids2016,dist_bed,Zsmooth,'2016')
makeSectionPlots(grids2017,dist_bed,Zsmooth,'2017')
makeSectionPlots(grids2018,dist_bed,Zsmooth,'2018')
makeSectionPlots(grids2019,dist_bed,Zsmooth,'2019')
makeSectionPlots(grids2021,dist_bed,Zsmooth,'2021')
makeSectionPlots(grids2022,dist_bed,Zsmooth,'2022')
makeSectionPlots(grids2023,dist_bed,Zsmooth,'2023')
makeSectionPlots(grids2023xctd,dist_bed,Zsmooth,'2023xctd')

% compareTSwshelf(grids2009,'2009',shelf2009)
% compareTS(grids2010,'2010')
% compareTSwshelf(grids2011,'2011',shelf2011)
% compareTS(grids2012,'2012')
% compareTSwshelf(grids2013,'2013',shelf2013)
% compareTSwshelf(grids2015,'2015',shelf2015)
% compareTS(grids2016,'2016')
% compareTSwshelf(grids2017,'2017',shelf2017)
% compareTSwshelf(grids2018,'2018',shelf2018)
% compareTSwshelf(grids2019,'2019',shelf2019)
% compareTSwshelf(grids2021,'2021',shelf2021)
% compareTS(grids2022,'2022')
% compareTSwshelf(grids2023,'2023',shelf2023)
% compareTS(grids2023xctd,'2023xctd')

compareTS(grids2009,'2009')
compareTS(grids2010,'2010')
compareTS(grids2011,'2011')
compareTS(grids2012,'2012')
compareTS(grids2013,'2013')
compareTS(grids2015,'2015')
compareTS(grids2016,'2016')
compareTS(grids2017,'2017')
compareTS(grids2018,'2018')
compareTS(grids2019,'2019')
compareTS(grids2021,'2021')
compareTS(grids2022,'2022')
compareTS(grids2023,'2023')
compareTS(grids2023xctd,'2023xctd')

