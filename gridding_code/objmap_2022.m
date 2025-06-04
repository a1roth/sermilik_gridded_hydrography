function objmap_2022(xcor1,ycor1,a1,xcor2,ycor2,a2,noise,xres,yres,xstart,xend,ystart,yend,notes)

%objmap_2022.m -

%   Input Variables Definitions
% variable: A string that will be part of saved names of output figures and
% logfile. Example '2019Temp'
% xx: x-coordinates (horizontal) of profile data. These can be in lat/lon
% or distance along a track
% yy: y-coordinates (vertical) of profile data. These are usually expressed
% as depth. 
% data: The profile data you wish to objectively map
% background: This is the background field that you wish to use. The mean
% of the field can be used if appropriate for your dynamics. 
% xcor1: X horizontal correlation length scale for first gaussian
% ycor1: Y vertical correlation length scale for first gaussian
% a1: First guassian amplitude 
% xcor2: X horizontal correlation length scale for second gaussian
% ycor2: Y vertical correlation length scale for second gaussian
% a2: Second guassian amplitude 
% noise: Noise parameter inherent to any observation. Larger noise the
% objective map will not fit to the data as closely. Smaller noise and the
% objective map will prioritize fitting the data closely over the
% smoothness of the domain.
% xres: resolution in x direction of final gridded data
% yres: resolution in y direction of final gridded data
% xstart: Final gridded data x axis starting point. Defining the limits of
% the gridded field domain.
% xend: Final gridded data x axis ending point.Defining the limits of
% the gridded field domain.
% ystart: Same as above but for the y-direction
% yend: Same as above but for the y-direction
% 
% notes: A text string that will appear in the output file information.
% Leave as '' if no notes. 

%% Load data

%The data you are loading should be manually selected profiles that create
%the best synoptic along-fjord section possible. These profiles should
%already have been processed (corrected for anomalies and spikes, binned,
%and smoothed!). They need to have lat and lon coordinates AND also
%along-fjord distances as "distance away from glacier" in km. This depends
%on your specific site and can be calculated with the function
%calc_distance_from_glacier.m if needed.

%Fields for the matfile should be
%lat - latitude on thalweg section in your distance coordinates
%lon - longitude on thalweg section in your distance coordinates
%CT - conservative temperature profiles
%SA - absolutely salinity profiles
%lat_raw - original latitude and longitude of profiles
%lon_raw - original longitude and longitude of profiles
%distance - distance in the along-fjord direction from the prescribed
%glacier terminus
%time - time of profiles
%depth - Depth coordinates for profiles. Should all be the same. 

%Import data for the given year
%as Matfile.

load('sf2022_input_profiles.mat')
CT = CT'; % transpose so that profiles are rows and depth is columns ie 15 profiles would be 15 x 901
SA  = SA';
dist_plot = dist; %try manually inputting
station  = stationID;

%or as netcdf
% lat = ncread('sf2022_input_profiles.nc','lat');
% lon = ncread('sf2022_input_profiles.nc','lon');
% CT = ncread('sf2022_input_profiles.nc','CT');
% SA = ncread('sf2022_input_profiles.nc','SA');
% station = ncread('sf2022_input_profiles.nc','stationID');
% dist_plot = ncread('sf2022_input_profiles.nc','distance');
% time = ncread('sf2022_input_profiles.nc','time');

%Also load bathymetry data of thalweg section that shows along-fjord
%thalweg coordinates that everything is plotted on. For Sermilik the file
%to use is "SermilikThalwegBathy.nc"
dist_bed = ncread('SermilikThalwegBathy.nc','dist'); %distance along thalweg section from Helheim 2019 terminus position
latZbed = ncread('SermilikThalwegBathy.nc','lat'); %latitude of points along thalweg section
lonZbed = ncread('SermilikThalwegBathy.nc','lon'); %longitude of points along thalweg section
Zsmooth = ncread('SermilikThalwegBathy.nc','depth'); %smoothed bathymetry depth values at Lat/Lon locations

%plot the thalweg section lat lons and the lat/lons of the profiles (on the
%thalweg section) that are being used for this objective map.
figure;
geoscatter(latZbed,lonZbed,5,'filled','k')
hold on
geoscatter(lat_raw,lon_raw,'b','filled') %shows original profile locations before put into alongfjord coordinates
geoscatter(lat,lon,'r','filled') %profile locations in distance coordinates on the thalweg section
geobasemap colorterrain

%% Fill in bottoms to be constant from deepest measurement. Identify which casts are deeper than 550m that can reasonable be extended to help grid the deeper parts of the fjord.
%If cast is shallow (ends above 550m) then don't fill in bottom, but if
%there is just a couple hundred meters to infill from 550-bottom then we can assume that the
%bottom is fairly constant.

whichcasts = [1,2,3,4,8,9]; %%which rows are deep enough to be extended. CHANGE THIS AS NEEDED FOR YOUR SPECIFIC SITE/YEAR. 

%% Extend data of selected casts

disp(strcat('Number of casts to be bottom-filled:  ',num2str(length(whichcasts))))

endindex = nan(size(whichcasts));
avg_bottom_value = nan(size(whichcasts));

CT_plot_bottomfill = CT;
SA_plot_bottomfill = SA;

%Find where NaNs start at each cast
for i = 1:length(whichcasts)
    c = whichcasts(i);
    nn = find(isnan(CT(c,10:end)),1,'first'); %have to start at index 3 because there are nans.
    endindex(i) = nn + 8;
    
    avg_bottom_valueCT(i) = mean(CT(c,(nn-7:nn + 2)));
    avg_bottom_valueSA(i) = mean(SA(c,(nn-7:nn + 2)));
    
    gg = 901 - endindex(i) + 1;
    
    CT_plot_bottomfill(c,endindex(i):end) = ones(1,gg).*avg_bottom_valueCT(i);
    SA_plot_bottomfill(c,endindex(i):end) = ones(1,gg).*avg_bottom_valueSA(i);
end

%endindex now shows which column the Nans start. 

%Take a quick look to make sure it looks okay
figure;
imagesc(CT_plot_bottomfill')


%%
%Vertially smooth data with moving mean to make gridding smoother. Can change these
%parameters if needed. Currently smoothing by 5m. 
CT_plot_smooth = smoothdata(CT_plot_bottomfill,2,'movmean',5,'omitnan'); %includenan or omitnan
SA_plot_smooth = smoothdata(SA_plot_bottomfill,2,'movmean',5,'omitnan');

%Add NaNs back in so that smoothdata isn't extending data beyond where we
%actually have measurements...
CT_plot_smooth(isnan(CT_plot_bottomfill))= NaN;
SA_plot_smooth(isnan(SA_plot_bottomfill))= NaN;



%% Objective map

%Adapted from Rudnick code offxy.m
%Objective Map of hydrographic profiles for temperature and salinity

%Figure out domain constants to set up grids. We need xg and yg which
%provide the horizontal and vertical gridded domain steps.

%Figure out the deepest actual measurement (not bottom filled). This
%determines the vertical extent that you can objectively map.
depthtest = nansum(CT,1); %This collapses all the profiles to find deepest extent
inds = find(depthtest);
yex=inds(end); %Maximum depth of this years objective map grid!

ny = floor((yex/yres)); %Calculates how many depth steps are needed at specified resolution 
ad1=-0.000; ad2=-.001*(yex-1); %depth range of objective map in kilometeres

%Figuring out how many horizontal grid steps are needed
numkm = floor(dist_plot(end)-dist_plot(1));
numkm = round(numkm/xres);

xg=linspace(dist_plot(1),dist_plot(end),numkm); %objective map along distance grid. Ensures a xres value spacing between first and last CTD cast. 
nx = numkm;

yg=linspace(ad1,ad2,ny); %objective map depth grid

%% SPECIFIC INPUT DATA TO BE OBJECTIVELY MAPPED and the individual coordinate pairs (distance, depth) for every data point in every profile. 

%Now specify input data to be objectively mapped

TempData = CT_plot_smooth'; %check if needed to fliplr or transpose..??
SalData = SA_plot_smooth';
[cc,vv] = size(TempData); 

%Now we need horizontal distance and vertical distance coordinate pairs for
%every data points. These need to be the same size as TempData
depth2=repmat(depth,1,vv)./(-1000); %Now in km.  Now in km and negative. %Needs to be same dimensions as TempData
xx2 = repmat(dist_plot,cc,1); %Needs to be same dimensions as TempData

%%
%Load background data for objective map
load('background_input_temp.mat')
background_CT = backgroundbin_CT_2km;
xback = xback_2km;
yback = yback;

load('background_input_sal.mat') %must be same dimensions as background_CT
background_SA = backgroundbin_SA_2km;

%Or load netcdf file if that's better
%background_CT = ncread('backgrounds.nc','background_CT')
%background_SA = ncread('backgrounds.nc','background_SA')
%xback = ncread('backgrounds.nc','dist')
%yback = ncread('backgrounds.nc','depth')

%% CREATE OBJECTIVE MAP FOR 2022 TEMPERATURE

%Make objective map
[objmap,errmap,residuals_mat,XG,YG] = objmap_from_profiles('2022T',xx2,depth2,TempData,xg,yg,background_CT,xback,yback,xcor1,ycor1,a1,xcor2,ycor2,a2,noise);

%% save objective map outputs into structure
grids2022.T.objmap = objmap;
grids2022.T.errmap = errmap;
grids2022.T.residuals = residuals_mat;
grids2022.T.xg = XG;
grids2022.T.yg = YG;

%save input data into structure
grids2022.T.data = TempData;
grids2022.T.xdata = xx2;
grids2022.T.ydata = depth2;
grids2022.T.latdata = lat;
grids2022.T.londata = lon;
grids2022.T.dist = dist_plot;

%Want some way to save non-bottomfilled points to data CT_plot
grids2022.T.data_raw = CT';
grids2022.T.lat_raw = lat_raw;
grids2022.T.lon_raw = lon_raw;
grids2022.T.info = 'data_raw is original profiles with no bottom-filled points. Only actually measured points. data is profiles with bottoms-filled in and what goes into the mapping algorith. objmap is the objective map of the profiles to the specified grid. errmap is the calculated uncertainty of the grid. residuals is the difference between input data (data) and objmap at the points of the input data.';

%% RUN OBJECTIVE MAP FOR SALINITY

[objmap,errmap,residuals_mat,XG,YG] = objmap_from_profiles('2022S',xx2,depth2,SalData,xg,yg,background_SA,xback,yback,xcor1,ycor1,a1,xcor2,ycor2,a2,noise);

%% Save salinity outputs into structure
grids2022.S.objmap = objmap;
grids2022.S.errmap = errmap;
grids2022.S.residuals = residuals_mat;
grids2022.S.xg = XG;
grids2022.S.yg = YG;
%save input data into structure
grids2022.S.data = SalData;
grids2022.S.xdata = xx2;
grids2022.S.ydata = depth2;
grids2022.S.latdata = lat;
grids2022.S.londata = lon;
grids2022.S.dist = dist_plot;


%Want some way to save non-bottomfilled points to data SA_plot
grids2022.S.data_raw = SA';
grids2022.S.lat_raw = lat_raw;
grids2022.S.lon_raw = lon_raw;
grids2022.S.info = 'data_raw is original profiles with no bottom-filled points. Only actually measured points. data is profiles with bottoms-filled in and what goes into the mapping algorith. objmap is the objective map of the profiles to the specified grid. errmap is the calculated uncertainty of the grid. residuals is the difference between input data (data) and objmap at the points of the input data.';


%% APPLY DATA CORRECTIONS

%CORREECT FOR TEMERATURE VALUES BELOW FREEZING POINT AND RESAVE
[objmapT_corr,residualsT_corr] = correctFreezing(grids2022.T.objmap,grids2022.S.objmap,grids2022.T.xg,grids2022.T.yg,grids2022.T.data,grids2022.T.xdata,grids2022.T.ydata,lat(1));

%save objective map outputs into structure
grids2022.T.objmap = objmapT_corr;
grids2022.T.residuals = residualsT_corr;
%% Interpolate onto standardized grid for all grids
%convert objective maps and error to standardized grid across whole domain
xq = xstart:xres:xend; %These are inputs into this function! 
yq = ystart:-1.*yres/1000:-1.*yend/1000; %in kilometers

[Xq,Yq] = meshgrid(xq,yq);

grids2022.T.errmap_stdgrid = interp2(grids2022.T.xg,grids2022.T.yg,grids2022.T.errmap,Xq,Yq);
grids2022.S.errmap_stdgrid = interp2(grids2022.S.xg,grids2022.S.yg,grids2022.S.errmap,Xq,Yq);

grids2022.T.objmap_stdgrid = interp2(grids2022.T.xg,grids2022.T.yg,grids2022.T.objmap,Xq,Yq);
grids2022.S.objmap_stdgrid = interp2(grids2022.S.xg,grids2022.S.yg,grids2022.S.objmap,Xq,Yq);

grids2022.T.Xq_stdgrid = Xq;
grids2022.S.Xq_stdgrid = Xq;

grids2022.T.Yq_stdgrid = Yq;
grids2022.S.Yq_stdgrid = Yq;

grids2022.notes = notes;
% 

%Save parameters for this run also
parameters = [{'xcor1'},xcor1;{'ycor1'},ycor1;{'a1'},a1;{'xcor2'},xcor2;{'ycor2'},ycor2;{'a2'},a2;{'noise'},noise;{'xres km'},xres;{'yres m'},yres];

grids2022.parameteres = parameters;

grids2022.raw_profile_stationID = station;
grids2022.raw_profile_date = time;
grids2022.lat_raw = lat_raw;
grids2022.lon_raw = lon_raw;

save('grids2022.mat','grids2022')
disp('grids2022 saved')
end
