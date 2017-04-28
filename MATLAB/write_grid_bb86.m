%% write grid double gyre

   clear all,clf, close all
   iodir='/Net/yucatan/abozec/BB86_PACKAGE/MATLAB/';
   addpath(genpath(['',iodir,'/UTILITIES/']));   
   
   %% PATH
   io = [iodir,'/../topo/'];
   file_grid_new = 'regional.grid.BB86'; %% !! without .a or .b !!

   %% size of the domain
   idm = 101 ;
   jdm = 101 ;

   %% longitude/latitude starting point + resolution  (in degrees) (NB: those
   %% variables are  not used in HYCOM)
   ini_lon = 0. ;
   ini_lat = 0. ;
   res = 0.20 ;

   %% scale dx and dy (in m) (used in HYCOM)
   dx = 20e3 ;
   dy = dx ;
   
   %% grid type
   mapflg=0 ; %% uniform or mercator


   %%%%%%; END of the USER inputs %%%%%%%%%%%%%%;

   %% missing value in HYCOM
   vmiss = 2.^100 ;

   %% Get the p-point the grid
   plon = zeros(jdm, idm) ;
   plat = zeros(jdm, idm) ;

   %% Longitude
   plon(:, 1) = ini_lon ;
   for i = 2:idm
     plon(:, i) = plon(:, i-1) + res;
   end
   %% Latitude
   plat(1, :) = ini_lat ;
   for j = 2:jdm
     plat(j, :) = plat(j-1, :) + res;
   end 

   disp('p-points grid OK')


   %% Declaration of the grid tabs
   idm_hd = idm ;
   jdm_hd = jdm ;
   
   plon_hd = plon ;
   plat_hd = plat ;

   qlon_hd = zeros(jdm_hd, idm_hd) ;
   qlat_hd = zeros(jdm_hd, idm_hd) ;
   ulon_hd = zeros(jdm_hd, idm_hd) ;
   ulat_hd = zeros(jdm_hd, idm_hd) ;
   vlon_hd = zeros(jdm_hd, idm_hd) ;
   vlat_hd = zeros(jdm_hd, idm_hd) ;

   pang_hd = zeros(jdm_hd, idm_hd) ;

   pscx_hd = zeros(jdm_hd, idm_hd) ;
   pscy_hd = zeros(jdm_hd, idm_hd) ;

   qscx_hd = zeros(jdm_hd, idm_hd) ;
   qscy_hd = zeros(jdm_hd, idm_hd) ;

   uscx_hd = zeros(jdm_hd, idm_hd) ;
   uscy_hd = zeros(jdm_hd, idm_hd) ;

   vscx_hd = zeros(jdm_hd, idm_hd) ;
   vscy_hd = zeros(jdm_hd, idm_hd) ;

   cori_hd = zeros(jdm_hd, idm_hd) ;
   pasp_hd = zeros(jdm_hd, idm_hd) ;

   
   %% Longitude/Latitude for each point (Q-points stay at the same location)
   %% longitude
   vlon_hd = plon_hd ;
   for i = 2:idm 
      qlon_hd(:, i) = 0.5*(plon_hd(:, i)+plon_hd(:, i-1));
   end
   diff = plon_hd(1, 3)-plon_hd(1, 2) ;
   qlon_hd(:, 1) = plon_hd(:, 1)-0.5*diff ;
   ulon_hd = qlon_hd ;
   
   %% latitude
   ulat_hd = plat_hd ;
   for j = 2:jdm
      qlat_hd(j, :) = 0.5*(plat_hd(j, :)+plat_hd(j-1, :));
   end 
   diff = plat_hd(3, 1)-plat_hd(2, 1) ;
   qlat_hd(1, :) = plat_hd(1, :)-0.5*diff ;
   vlat_hd = qlat_hd ;
 
   %% simplified grid with prescribed dx and dy
   pscx_hd(:, :) = dx;
   pscy_hd(:, :) = dy;
   uscx_hd(:, :) = dx; 
   uscy_hd(:, :) = dy;
   vscx_hd(:, :) = dx; 
   vscy_hd(:, :) = dy;
   qscx_hd(:, :) = dx; 
   qscy_hd(:, :) = dy;
   
   %% Coriolis
   beta = 2.e-11 ;
   for  j= 1:jdm_hd 
     for  i= 1:idm_hd 
      cori_hd(j, i)=.93e-4+double(j-(jdm_hd-1)/2)*dx*beta  ; 
     end 
   end  
    
   %% pang ( p-angle  for rotated grid)
   %% here uniform or mercator so pang = 0.


   %% pasp= paspect: pscx/pscy
   
   for j= 1:jdm_hd
     for i= 1:idm_hd
       pasp_hd(j,i) = pscx_hd(j,i)/pscy_hd(j,i)  ;
     end  
   end 
   
   
   %% Writing new grid file
   write_grid_hycom(idm_hd, jdm_hd, io, file_grid_new, plon_hd, plat_hd, ulon_hd,  ulat_hd,  vlon_hd, vlat_hd, ...
		    qlon_hd, qlat_hd, pang_hd,pscx_hd, pscy_hd, qscx_hd,  qscy_hd, uscx_hd, uscy_hd,  vscx_hd, vscy_hd, cori_hd, pasp_hd,mapflg)

   disp('Writing grid file done ')


   
   
