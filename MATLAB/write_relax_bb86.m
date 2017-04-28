%% write relax files (initial condition of T, S and interface
%depth)


   clear all,clf, close all
   iodir='/Net/yucatan/abozec/BB86_PACKAGE/MATLAB/';
   addpath(genpath(['',iodir,'/UTILITIES/']));   
   
   %% PATH
   io = [iodir,'../relax/010/'];
   file_topo = 'depth_BB86_01.a';
   file_int_new = 'relax_int_BB86' ; %% !! without .a or .b !!
   file_tem_new = 'relax_tem_BB86' ;
   file_sal_new = 'relax_sal_BB86' ;
   pl = 0; %% 1 or 0 for plot or not

   %% size of the domain
   idm = 101 ;
   jdm = 101 ;
   kdm = 3 ; %% number of layers (first layer very thin)
   tdm = 12 ; %% 12 month climatology

   %% interface depth (!!!first interface depth always 0. !!!)
   id = zeros(kdm) ;
   id = [0., 1., 500.] ;


   %% target density
   %% from bb86: rho=27.01037, rho=27.22136
   d = zeros(kdm) ;
   d = [27.0100, 27.01037,27.22136] ;
   
   %% salinity profile
   sa = zeros(kdm) ;
   sa = [37., 37.,  37.] ;

   %%%%%%; END of the USER inputs %%%%%%%%%%%%%%;

   %% constants
   rho = 1000.;
   g = 9.806;
   vmiss = 2.^100;

   %% temperature profile for sigma0
   sigma = 0 ; 
   te = 1:kdm ;
   for k = 1:kdm 
     te(k) = tofsig(sa(k), d(k), sigma) ;
   end 
   disp(te)

   
  %% read bathy for mask
   bathy=read_depth_hycom(idm, jdm, [io,'../../topo/',file_topo]);
   ind = isnan(bathy) ;
   bathy(ind) = 0. ;
   %% create mask
   mask2d = ones(jdm, idm) ;
   mask2d(ind) = 0. ;
   mask = zeros(jdm, idm, kdm, tdm) ;
   for t = 1:tdm
   for k = 1:kdm
     mask(:, :, k, t) = mask2d ;
   end 
   end 
   
   %% interface depth files (first layer always the surface i.e. 0.)
  int = zeros(jdm, idm, kdm, tdm) ;
  for k= 2:kdm
    int(:, :, k,:) = rho*g * id(k) ; %% 1st layer tiny to fake a two layer config
  end 
    
  %% make sure that the interface depths are not lower than the bathy
  for t = 1:tdm
  for k = 1:kdm
     for i = 1:idm 
        for j = 1:jdm
           if (int(j, i, k)/9806. > bathy(j, i)) 
	     int(j, i, k) = bathy(j,i)*9806. ;
	   end 
        end 
     end 
  end 
  end 

  %% mask 
  ind3d = find(mask == 0.) ;
  int(ind3d) = vmiss ;

  
  %% write the field in relax file
  write_relax_hycom(idm, jdm, kdm, io, file_int_new, 'intf', d, int) ;
  disp('Interface depth OK')

   %% Temperature
  tem = zeros(jdm, idm, kdm, tdm) ;
  for k = 1:kdm 
    tem(:, :, k,:) = te(k) ;
  end 
 
  %% mask 
  tem(ind3d) = vmiss;

  %% write the field in relax file
  write_relax_hycom(idm, jdm, kdm, io, file_tem_new, 'temp', d, tem) ;
  disp('Temperature OK')

  %% Salinity
  sal = zeros(jdm, idm, kdm, tdm) ;
  for k = 1:kdm  
    sal(:,:, k,:) = sa(k) ;
  end 

  %% mask 
  sal(ind3d) = vmiss ;

  %% write the field in relax file
  write_relax_hycom(idm, jdm, kdm, io, file_sal_new, 'saln', d, sal);
  disp('Salinity OK')



   %% Plot
   if (pl == 1) 
     figure(1)
     pcolor(tem(:,:,1));colormap(jet(length(1:32)));
     colorbar;shading flat     
   end 
