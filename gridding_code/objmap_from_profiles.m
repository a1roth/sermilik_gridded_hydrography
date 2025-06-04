function [objmap,errmap,residuals_mat,XG,YG] = objmap_from_profiles(id,xx_mat,yy_mat,data,xg,yg,background,xb,yb, xcor1,ycor1,a1,xcor2,ycor2,a2,noise)
%objmap_from_profiles.m - Takes in CTD profiles and creates a gridded,
%smoothed field using objective mapping (optimal interpolation, optimum
%analysis) methods. 
%Author: Aurora Roth
%Version 1 February 14 2024

%Code was intially adapted from Dan Rudnick's SIOXXX
%class code offxy.m Further code modification with the help of Matthew
%Mazzloff (SIO). 

%

%   Input Variables Definitions
% id: A string that will be part of saved names of output figures and
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
% smoothness of the map.

disp(strcat('Starting objective map for: ',id,'...'));

xx = xx_mat(:); %vector format of horizontal locations
yy = yy_mat(:); %vector format of vertical locations


tic
xcor1sq=xcor1*xcor1;
ycor1sq=ycor1*ycor1;

xcor2sq=xcor2*xcor2;
ycor2sq=ycor2*ycor2;
phi = 0; %no rotation in this grid
cosphi=cos(phi/180*pi); %rotation won't matter here
sinphi=sin(phi/180*pi); %same
err = noise;

%now have to make background fit to the model grid.
[xxback,yyback] = meshgrid(xb,yb);
background_grid = interp2(xxback,yyback,background,xg, yg');

%Now get background at CTD profile data locations
bgrnd_ctd = interp2(xg,yg,background_grid,xx,yy);

%data needs to be in vector form and subtracting background field
data = data(:)-bgrnd_ctd;  

%Put data into vector form and get rid of NaN values
ii=find(~isnan(data(:,1)));
d=data(ii); %vector format of data
x = xx(ii); %get rid of nan values
y=yy(ii); %get rid of nan values

nfact=length(d); %get length of vector of data
%At this point, the data, horizontal coordinates, and vertical coordinates, are all
%in 3 diferent long single column vectors, with all nan data points
%removed.

normdata = std(d);
meandata = mean(d);
d = (d-meandata)/normdata; %nondimensionalized d (demeaned and normalized)

%Calculate the data-data covariance matrix E.
[X1,X2]=meshgrid(x);
[Y1,Y2]=meshgrid(y);
dx=X1-X2;
dy=Y1-Y2;
dxr=dx*cosphi-dy*sinphi; %doesn't matter because of no rotation
dyr=dx*sinphi+dy*cosphi; %doesn't matter because of no rotation

%Set up gaussians from correlation scales
E1=a1.*exp(-dxr.^2/xcor1sq-dyr.^2/ycor1sq)+a2.*exp(-dxr.^2/xcor2sq-dyr.^2/ycor2sq); %Can specify multiple guassians here
E=E1+(err*err)*eye(nfact);
toc %7 seconds
%Calculate the data-model points covariance matrix C.
C=zeros(nfact,length(xg)*length(yg));

[Y,X]=meshgrid(yg,xg);
for n=1:nfact
   dx=X(:)'-x(n);
   dy=Y(:)'-y(n);
   dxr=dx*cosphi-dy*sinphi;
   dyr=dx*sinphi+dy*cosphi;
   
   C(n,:)=a1.*exp(-dxr.^2/xcor1sq-dyr.^2/ycor1sq)+a2.*exp(-dxr.^2/xcor2sq-dyr.^2/ycor2sq);%+ampb.*exp(-dxr.^2/xcor2-dyr.^2/ycor2); %specify multiple guassians
end
toc %27 seconds
%Look at condition number. Range that you're dividing by. 
eigval = eig(E);
condnum = min(eigval)/max(eigval);
disp(strcat('Condition Number: ', num2str(condnum)));
toc %50 seconds
wtddat = E\d; %project data into data-data matrix Same as inv(E)*d

dhat = C'*wtddat; %projecting data information into model land 

%Formal mapping error - but matrices are too large for matlab
%M = diag((C'*(E^-1)*C)); %formal mapping error, based on locations of data and model grid locations
%M = diag((C'*(E\C)));
toc %64 seconds
EC = E\C; %this step takes a little bit of time depending on size of target grid. Same as inv(E)*c
toc %305
for i = 1:length(dhat)
    u(i) = sum(EC(:,i).*C(:,i)); %How to interpret this number??
end
toc %368
dhat = (dhat*normdata)+meandata; %degrees C units 
u = sqrt((a1+a2) - u)*normdata; %degrees C units &have sum of amps here sqrt((ampa + ampb)-u)
%Amplitudes might have to squared in above code. Or all amplitudes have to
%sum to 1. Then don't worry about squaring. 

objmap = reshape(dhat,length(xg),length(yg)) + background_grid';

%Get objective map values at CTD profile locations
objmap_ctdlocations = interp2(xg,yg,objmap',xx,yy);

residuals = (data(:)+bgrnd_ctd)-objmap_ctdlocations;  %data needs to be in vector form and subtracting background field

errmap = reshape(u,length(xg),length(yg));
toc
%% Plots to show objective map and original data
%Plot just objective map. 
figure;
set(gcf,'position',[10,10,700,1200])
subplot(3,1,1)
    pcolor(xg,yg.*1000,objmap')
    shading flat
    colorbar; 
    %caxis([-1 5])
    hold on
    scatter(xx,yy.*1000,30,data+bgrnd_ctd,'filled')
    xlabel('distance [km]')
    ylabel('depth [m]')
    set(gca,'fontsize',14)
    title([ 'xcor_a=' num2str(xcor1) '; ycor_a= ' num2str(ycor1) '; amp_a= ' num2str(a1) ';xcor_b=' num2str(xcor2) '; ycor_b= ' num2str(ycor2) '; amp_b= ' num2str(a2) '; error= ' num2str(noise) ])

%figure;
subplot(3,1,2)
    pcolor(xg,yg.*1000,errmap')
    shading flat
    colorbar; 
    hold on
    scatter(x,y.*1000,'k.')
    xlabel('distance [km]')
    ylabel('depth [m]')
    set(gca,'fontsize',14)
    title([ 'xcor_a=' num2str(xcor1) '; ycor_a= ' num2str(ycor1) '; amp_a= ' num2str(a1) ';xcor_b=' num2str(xcor2) '; ycor_b= ' num2str(ycor2) '; amp_b= ' num2str(a2) '; error= ' num2str(noise) ])
    
%figure;
%set(gcf,'position',[10,10,950,370])
subplot(3,1,3)
    colorbar; 
    hold on
    scatter(xx,yy.*1000,30,residuals,'filled')
    xlabel('distance [km]')
    ylabel('depth [m]')
    set(gca,'fontsize',14)
    %caxis([-0.8 0.8])
    title([ 'xcor_a=' num2str(xcor1) '; ycor_a= ' num2str(ycor1) '; amp_a= ' num2str(a1) ';xcor_b=' num2str(xcor2) '; ycor_b= ' num2str(ycor2) '; amp_b= ' num2str(a2) '; error= ' num2str(noise) ])

% Define the file name for the PNG file
%timestamp = datenum(datetime('now'));
%file_name = strcat('objmap_',num2str(round(timestamp)),'.png');
file_name = strcat('objmap_',id,'.png');


% Save the plot as a PNG file in the specified folder
saveas(gcf, file_name);

%% Plot to look at residuals and skill of objective map in matching original data values
figure; 
set(gcf,'position',[10,10,650,220])

subplot(1,2,1)
scatter(data+bgrnd_ctd,objmap_ctdlocations,40,yy.*1000,'filled')
hold on
plot([-2:6],[-2:6],'--k')
grid on
colorbar
%caxis([-900 0])

subplot(1,2,2)
scatter(data+bgrnd_ctd,objmap_ctdlocations,40,xx,'filled')
hold on
plot([-2:6],[-2:6],'--k')
grid on
colorbar
%caxis([0 104])

%file_name = strcat('residuals_',num2str(round(timestamp)),'.png');
file_name = strcat('residuals_',id,'.png');

% Save the plot as a PNG file in the specified folder
saveas(gcf, file_name);


%% Write out stats of residuals to log file

%I want to add here a TS plot of original data and objective map data
%Stats about the original data - the mean, std, max, and min, and
%respective coordinates of max and min values. 

mask = ~isnan(residuals);
residuals_nonan = residuals(mask);

datamean = nanmean(data+bgrnd_ctd);

rstd = nanstd(residuals);
rmean = nanmean(residuals);
rmax = nanmax(residuals);
rmin = nanmin(residuals);
rmse = sqrt(nansum(residuals_nonan.^2)/length(residuals_nonan));
r2 = 1 - (nansum(residuals.^2)/nansum((data+bgrnd_ctd - datamean).^2));

disp(strcat('rstd: ',num2str(rstd)));
disp(strcat('rmean: ',num2str(rmean)));
disp(strcat('rmax: ',num2str(rmax)));
disp(strcat('rmin: ',num2str(rmin)));
disp(strcat('rmse: ',num2str(rmse)));
disp(strcat('r2: ',num2str(r2)));
max_err = max(max(errmap));
disp(strcat('max_err: ', num2str(max_err)));

% Values to write into file 

to_write = [datenum(datetime('now')), xcor1,ycor1,a1...
    xcor2,ycor2,a2,err,condnum, max_err,rstd, rmean, rmax, rmin, rmse, r2];

% Define the file name
file_name = strcat('logfile_',id,'.txt');
% Open the file in append mode

fileID = fopen(file_name, 'w');

% Write the new values (including timestamp) to the file as a comma-separated row
for i = 1:length(to_write)-1
    if isnumeric(to_write(i))
        fprintf(fileID, '%11.7f,', to_write(i));
    else
        fprintf(fileID, '%s,', char(to_write(i)));
    end
end
% Write the last value and start a new line
if isnumeric(to_write(end))
    fprintf(fileID, '%d\n', to_write(end));
else
    fprintf(fileID, '%s\n', char(to_write(end)));
end

% Close the file
fclose(fileID);

disp(['Stats written to log file:  "', file_name, '"']);

[XG,YG] = meshgrid(xg,yg);
objmap = objmap';
errmap = errmap';
residuals_mat = reshape(residuals,(size(xx_mat)));

end

