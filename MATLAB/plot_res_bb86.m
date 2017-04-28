%% plot results from bb86 and HYCOM-bb86

   clear all,clf, close all
   iodir='/Net/yucatan/abozec/BB86_PACKAGE/MATLAB/';
   addpath(genpath(['',iodir,'/UTILITIES/']));   

   %% PATH
   io = [iodir,'../'];

   %% size of the domain
   idm = 101 ;  %% size of the domain
   jdm = 101 ;  %% size of the domain
   kdm = 2 ;    %% number of layer in bb86
   dp0 = 500.;  %% thickness of the first layer (m)
   tdm = 1800 ; %% 12 month climatology
   tplot1 = 5 ; tplot2 = 5 ; %% time-stamp to plot (starts from 1) 
   
   %% Postscript
   PS=0;
   file_ps='../PS/test.eps'
   
   %% constants
   rho = 1000.;   %% reference density
   g = 9.806  ;   %% gravity
   

   %% Read grid
   file_grid = 'regional.grid.BB86.a'
   [plon,plat]=read_grid_hycom(idm, jdm, ([io,'topo/']), file_grid);
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% READ the HYCOM files  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   io_hycom = [io,'expt_01.0/data/output/'];
   
   for t = tplot1:tplot2 
      if (t <= 344) 
            year = '0001' ;
            day = sprintf('%3.3d',t+16)  ; %% 16 u and v == 0. so we start at day 17.
      elseif (t > 344) && (t <= 704) 
            year = '0002' ;
            day = sprintf('%3.3d',t-344) ;
      elseif (t > 704) && (t <=  1064) 
            year = '0003' ;
            day = sprintf('%3.3d',t-704) ;
      elseif (t > 1064) && (t <=  1424) 
            year = '0004' ;
            day = sprintf('%3.3d',t-1064) ;
      elseif (t > 1424) && (t <=  1784) 
            year = '0005' ;
            day = sprintf('%3.3d',t-1424) ;
      elseif (t > 1784) && (t <=  2144) 
            year = '0006' ;
            day = sprintf('%3.3d',t-1784) ;
      end  
      file = ['archv.',year,'_',day,'_00.a']
      

      %% archive files 2D variable units (see archv.*.b for the list of available var.)
      %%      montg1  : *1./g to get in m
      %%      srfhgt  : *1./g to get in m
      %%      surflx  : W/m2
      %%      salflx  : kg/m2/s
      %%      bl_dpth : *1./(rho*g) to get in m
      %%      mix_dpth: *1./(rho*g) to get in m
      %%      u_btrop : m/s
      %%      v_btrop : m/s

      %% extract barotropic velocities (ubarot)
      num2 = 8   ;  %% number of 2D variables in the archive file
      ivar_ub = 7;  %% index of ubaro (7th 2D variable)
      ivar_vb = 8;  %% index of vbaro (8th 2D variable)
      ubaro=sub_var2([io_hycom,file],[jdm, idm], num2, ivar_ub) ;
      vbaro=sub_var2([io_hycom,file],[jdm, idm], num2, ivar_vb) ;  

      
      
      %% archive files 3D variable units (see archv.*.b for the list of available var.)
      %%      u-vel.  : m/s
      %%      v-vel.  : m/s
      %%      thknss  : *1./(rho*g) to get in m
      %%      temp    : C
      %%      salin   : psu
      
      %% extract velocities (ubaroc)
      num3 = 5  ; %% number of 3D variables in the archive file
      ivar_u = 1; %% index of u baroclinic
      ivar_v = 2; %% index of v baroclinic
      ubac=sub_var3([io_hycom,file],[jdm, idm, kdm+1], num2, num3, ivar_u) ;
      vbac=sub_var3([io_hycom,file],[jdm, idm, kdm+1], num2, num3, ivar_v) ;   

      %% Get utot & vtot
      uhyc=zeros(jdm,idm,kdm+1,tdm) ;
      for k = 1:kdm+1
	uhyc(:, :, k, t) = ubaro(:, :)+ubac(:, :,k) ;
      end
      vhyc=zeros(jdm,idm,kdm+1,tdm) ;
      for k = 1:kdm+1
	vhyc(:, :, k, t) = vbaro(:, :)+vbac(:, :, k) ;
      end


      %% extract the layer thickness
      ivar_dp = 3 ;
      dp=sub_var3([io_hycom,file],[jdm, idm, kdm+1], num2, num3, ivar_dp) ;
      dphyc=zeros(jdm,idm,kdm+1,tdm) ;
      dphyc(:, :, :, t) = dp/(rho*g) ;

      %% put the u vel on the p-grid
      %%get the umask
      maskt = ones(jdm, idm, kdm+1, tdm) ;
      maskt(find(uhyc == 0.)) = 0. ;
      maxval = maskt+circshift(maskt, [0, -1, 0, 0]) ;
      maxval(find(maxval == 0.)) = NaN ;

      %% get the average value of u at the p-point
      uthyc = (uhyc+circshift(uhyc, [0, -1, 0, 0]))./maxval ;
      uthyc(:,idm, :, :) = uthyc(:,idm-1, :, :) ;
   
      %% put the v vel on the p-grid
      %%get the vmask
      maskt = ones(jdm, idm, kdm+1, tdm) ;
      maskt(find(vhyc == 0.)) = 0. ;
      maxval = maskt+circshift(maskt, [-1, 0, 0, 0]) ;
      maxval(find(maxval == 0.)) = NaN ;

      %% get the average value of v at the p-point
      vthyc = (vhyc+circshift(vhyc, [-1, 0, 0, 0]))./maxval ;
      vthyc(jdm, :, :, :) = vthyc(jdm-1,:, :, :) ;
   
    end    
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %% plot 
      for t = tplot1:tplot2 

         %% norme, min and max of the plots
         normeref_hyc  = 0.02 ; %% vector norm of HYCOM-BB86 (in m)
         onevectout = 2       ; %% plot one vector out of 'onevectorout'
         min_dp = -5.       ;
	 max_dp = 5.        ; %% layer thickness anomaly in m (min and max)

	 %% Define the figure dimensions
	 set(gcf,'Units','normalized');
	 set(gcf,'position',[0.1 0.08,0.5,0.7])
	 set(gcf,'PaperType', 'usletter');
	 set(gcf,'PaperUnits','normalized');
	 set(gcf,'PaperPosition',[0,0,1,1]);
	 set(gcf,'PaperOrientation', 'landscape');

	 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %% HYCOM-BB86 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
	 subplot(2,2,1)
         uu = uthyc(1:onevectout:jdm, 1:onevectout:idm, 2, t);
         vv = vthyc(1:onevectout:jdm, 1:onevectout:idm, 2, t);
         lon = plon(1:onevectout:jdm, 1:onevectout:idm);
         lat = plat(1:onevectout:jdm, 1:onevectout:idm);
	 [m,n] = size(uu);

         %% we normalize the vectors
         normeref = normeref_hyc ; %% in m
         norme =sqrt(uu.^2+vv.^2) ;
         normemax = nanmax(nanmax(norme))/normeref ;

	 quiver(lon, lat, uu,vv,3*normemax,'k')   %% ref vector over 3*dx
	 axis;
	 set(gca,'xlim',[0, 20],'ylim',[0,20]);
	 title(['HYCOM-BB86  Velocity field  t = ',sprintf('%4.4d',t)]);
	 xlabel('Longitude (E)')
	 ylabel('Latitude (N)')

	 %% draw a reference vector
	 framelim = get(gca,'Position');	 
	 x=framelim(1) ; y=framelim(2) ; %% (x,y) position of the arrow in the figure(normalized) units
	 xpl=1.2       ; ypl=-1.8      ; %% (x,y) position of the arrow in the plot units
	
	 x1=x+xpl*framelim(3)/(lon(1,n)-lon(1,1));
	 y1=y+ypl*framelim(4)/(lat(m,1)-lat(1,1));
	 dx=3*(lon(1,2)-lon(1,1))*framelim(3)/(lon(1,n)-lon(1,1));
	 dy=0.;
	 annotation(figure(1),'arrow','Position',[x1,y1,dx,dy],'headstyle','vback3','headwidth',5,'headlength',5)
	 text(xpl,ypl-1.,[sprintf('%3.2f',normeref),'m'])
	

         
	 subplot(2,2,2)

         levels = 40 ;
         Minss = min_dp ;
	 Maxss = max_dp ;
         step = (Maxss - Minss) / levels ;
	 toto = 1:levels;
         num_level = toto * step + Minss ;

	 caxis([min_dp max_dp]);hold on;
         diff_dp = dphyc(:, :, 2, t) - (dp0+1.) ;
	 contourf(plon,plat,diff_dp,num_level);hold on;
	 colorbar;
	 axis;
	 title(['HYCOM-BB86 thickness ano.   t = ',sprintf('%4.4d',t)]);
	 xlabel('Longitude (E)')
	 ylabel('Latitude (N)')
	
	 end
	 if (PS == 1)
	   print('-depsc2',file_ps);
	 end
   % end 
	 
	 
	 
