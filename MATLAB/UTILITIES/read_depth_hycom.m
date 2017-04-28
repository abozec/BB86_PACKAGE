function bathy=read_depth_hycom(im,jm,file)

  %% Script to read the HYCOM bathymetry file
  %% A. Bozec & D. Dukhovskoy Aug, 2011

  %% Get the id of the bathymetry file
  depth_fid=fopen(file,'r');
  IDM=im;
  JDM=jm;
  IJDM=IDM*JDM;
  npad=4096-mod(IJDM,4096);
  
  %% Read bathymetry 
  [bathy,count]=fread(depth_fid,IJDM,'float32','ieee-be');

  %% Mask bathymetry
  y=find(bathy>1e20);
  bathy(y)=NaN;
  bathy=reshape(bathy,IDM,JDM)';

  fclose(depth_fid);
 
 
 
 




