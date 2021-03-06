%% write bathy double gyre

   clear all,clf, close all
   iodir='/Net/yucatan/abozec/BB86_PACKAGE/MATLAB/';
   addpath(genpath([iodir,'/UTILITIES/']));   
   
   %% PATH
   io = [iodir,'/../topo/'];
   file_bat_new = 'depth_BB86_01' ; %% !! without .a or .b !!
   pl = 0;  %% 1 or 0 for plot or not

   %% size of the domain
   idm = 101 ;
   jdm = 101 ;

   %% definition of the bathymetry
   depth = 5000. ;%% constant bathy everywhere

   %% closed boundary (latbdy=0) or open boundaries (latbdy=1) or cyclic (latbdy =2)
   latbdy = 0;

   %%%%%%; END of the USER inputs %%%%%%%%%%%%%%;

   idm_hd = idm;
   jdm_hd = jdm;
   vmiss = 2.^100;  %% HYCOM missing values

   %% get the Depth 
   bathy_hd = zeros(jdm_hd, idm_hd);
   bathy_hd(:, :) = depth;
%    bathy_hd(:,1)=0.;
%    for i=2:50
%      bathy_hd(:,i)=bathy_hd(:,i-1) + 10.;
%    end 
   disp('Depth Ok')

   %% Plot
   if (pl == 1) 
     figure(1)
     pcolor(bathy_hd(:,:));colormap(jet(length(1:32)));
     colorbar;shading flat     
   end 

   %% Mask by missing values if any
   ind = find(bathy_hd == 0.);
   bathy_hd (ind) = vmiss;


   %% Get the edge of the domain right
   switch latbdy
      case 0
         bathy_hd(jdm_hd, :) = vmiss ;
         bathy_hd(:, idm_hd) = vmiss ;
         bathy_hd(1, :) = vmiss ;
         bathy_hd(:, 1) = vmiss ;

      case 1
         bathy_hd(jdm_hd, :) = vmiss ;
         bathy_hd(:, idm_hd) = vmiss ;

      case 2
         bathy_hd(jdm_hd, :) = vmiss ;
         bathy_hd(     1, :) = vmiss ;
   end


   %% write the file
   write_depth_hycom(idm_hd, jdm_hd, [io,file_bat_new], bathy_hd)
   disp('Writing depth file done ')

