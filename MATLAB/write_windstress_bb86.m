%% write windstress double gyre

   clear all,clf, close all
   iodir='/Net/yucatan/abozec/BB86_PACKAGE/MATLAB/';
   addpath(genpath(['',iodir,'/UTILITIES/']));   
   
   %% PATH
   io = [iodir,'/../topo/'];
   file_grid = 'regional.grid.BB86.a'

   %% domain
   idm = 101 ;
   jdm = 101 ;
   tdm = 12 ;%% monthly files

   %% name the new files
   file_E = 'forcing.tauewd.BB86';  %% without .a or .b %%
   file_N = 'forcing.taunwd.BB86';


   %% Read grid
   [plon,plat]=read_grid_hycom(idm, jdm, ([io,'../topo/']), file_grid);

   %% Calculation of the analytical wind-stress (N/m2)
   ustress = zeros(jdm);
   stressa = -1.; 
   sconv = 1.e-1;   %% scale factor form dyn/cm2 to N/m2

   %% BB86 formulation
   for j = 1:jdm
     ustress(j) = stressa*cos(double(j-1)/double(jdm-1)*6.28318530718)*sconv;
   end

   %% plot the stress 
   figure(1)
   plot(ustress,plat);
   ylim([0 20]);
   xlim([-0.2 0.2]);
   grid on;
   

   %% Write the wind-stress files
   tte = zeros(jdm, idm, tdm);
   ttn = zeros(jdm, idm, tdm);
  
   IJDM=idm*jdm;
   npad=4096-mod(IJDM,4096);
   toto=zeros(npad,1);

 
   %% we apply  a taux , no tauy
   for j = 1:jdm 
     tte(j, :, :) = ustress(j);
   end  

   %% Taux
   taux_fid=fopen([io,'../force/',file_E,'.a'],'w');   

   for ll = 1:tdm 
     A=tte(:,:,ll)';
     A=reshape(A,IJDM,1);
     %% Writing the field
     fwrite(taux_fid,A,'float32','ieee-be');
     %% Writing the padding at the end of the record:
     fwrite(taux_fid,toto,'float32','ieee-be');
   end 
   fclose(taux_fid);
   
   %% Tauy
   tauy_fid=fopen([io,'../force/',file_N,'.a'],'w');   

   for ll = 1:tdm 
     A=ttn(:,:,ll)';
     A=reshape(A,IJDM,1);
     %% Writing the field
     fwrite(tauy_fid,A,'float32','ieee-be');
     %% Writing the padding at the end of the record:
     fwrite(tauy_fid,toto,'float32','ieee-be');
   end 
   fclose(tauy_fid);
   

   %% create .b file 
   time = 1:12;
   
   %% Write Eastward Wind stress file
   fldb=[io,'../force/',file_E,'.b'];
   fid1=fopen(fldb,'wt');

   fprintf(fid1, 'Analytical Eastward Wind-stress\n');
   fprintf(fid1, '\n');
   fprintf(fid1, '\n');
   fprintf(fid1, '\n');
   fprintf(fid1, 'i/jdm =  %i  %i \n',idm,jdm);
   for m = 1:tdm 
      fprintf(fid1,' tauewd: month,range =  %2.2i  %10.5E  %10.5E\n', ...
		   time(m), min(min(tte(:, :, m)')), max(max(tte(:, :, m)')));
   end 
        
   fclose(fid1);

   %% Write Northward Wind stress file
   fldb=[io,'../force/',file_N,'.b'];
   fid1=fopen(fldb,'wt');

   fprintf(fid1, 'Analytical Northward Wind-stress\n');
   fprintf(fid1, '\n');
   fprintf(fid1, '\n');
   fprintf(fid1, '\n');
   fprintf(fid1, 'i/jdm =  %i  %i \n',idm,jdm);
   for m = 1:tdm 
      fprintf(fid1,' taunwd: month,range =  %2.2i  %10.5E  %10.5E\n', ...
		   time(m), min(min(ttn(:, :, m)')), max(max(ttn(:, :, m)')));
   end 
        
   fclose(fid1);

  