function [objmapT_corr,residualsT_corr] = correctFreezing(objmapT,objmapS,xg,yg,dataT,xd,yd,lat)

%correctFreezing takes in the objectively mapped grid for temperature and
%salinity and adjusts the temperature grid so that no temperature values
%are below the freezing point for a given density and pressure. Then the
%residuals for the temperature objective map have to be recalculated with
%this correction of the temperature. 

%This function uses the Gibbs Sea Water toolbox of TEOS-10 calculations.
%This must be downloaded and in the path prior to using this function. 

%INPUT
%objmapT = The objectively mapped gridded temperature field that needs
%correction
%objmapS = The objectively mapped gridded salinity field that matches the
%exact dimensions of the temperature field
%xg = x-dimensions horizontal of objmapT and objmapS. Units are in km.
%Dimensions match the dimensions are objmapT and objmapS. In meshgrid
%format. 
%yg = y-dimensions vertical coordinatores of objmapT and objmap S. Units
%are in depth (m).Dimensions match the dimensions are objmapT and objmapS. In meshgrid
%format. 
%dataT = original data used as input to create the objective map
%temperature field
%xd = horizontal coordinates of dataT. Units are in km. In meshgrid format,
%so dimensions match dataT. 
%yd = vertical coordinates of dataT. Units are in depth (m).In meshgrid format,
%so dimensions match dataT. 

[ng,mg] = size(objmapT);

%Change depth coordinates to pressure. Using GSW Toolbox function and a
%representative latitude for the fjord region. 
presg = gsw_p_from_z(yg,lat);

%Calculate freezing temperature at each point in objective map grid
TF_CT = gsw_CT_freezing(objmapS,presg); %Units conservative temperature


%loop through each cell of temperature values in objective map and use 
%corresponding salinity values at each point to then calculate the freezing
%temperature. If the value of the temperature grid cell is below the
%freezing value, then change it to the freezing value.

objmapTc = objmapT;

%variable to keep track number of corrected data points
numtracker = 0;

for i = 1:ng
    for j = 1:mg
        if objmapT(i,j)<TF_CT(i,j)
            objmapTc(i,j)=TF_CT(i,j);
            disp('Corrected temperature value to freezing point')
            numtracker = numtracker + 1;
        end
    end
end

objmapT_corr = objmapTc; %save output variable of corrected map. 

%Recalculate temperature residuals. Residuals = originalTdata - objmap

%First interpolate objective map values to original data locations
%Get objective map values at CTD profile locations
objmap_ctdlocs = interp2(xg,yg,objmapT,xd,yd);

residualsT_corr = dataT-objmap_ctdlocs;

disp(strcat('Number of data points corrected to freezing temp:  ',num2str(numtracker)))


end

