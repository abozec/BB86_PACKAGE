function [plon, plat, ulon,  ulat,  vlon, vlat, qlon, qlat, pang, ...
	  pscx, pscy, qscx,  qscy, uscx, uscy,   vscx, vscy, cori, pasp]=read_grid_hycom(im, jm,  io, file)

  %% Script to read the HYCOM grid file
  %% A. Bozec & D. Dukhovskoy Aug, 2011
    
  %% Get the id of the grid file
  flgr=[io,file];
  grid_fid=fopen(flgr,'r');
  IDM=im;
  JDM=jm;
  IJDM=IDM*JDM;
  npad=4096-mod(IJDM,4096);
  
  %% Read the grid

  [plon,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  fseek(grid_fid,4*(npad+IJDM),-1);
  [plat,count]=fread(grid_fid,IJDM,'float32','ieee-be');

  fseek(grid_fid,2*4*(npad+IJDM),-1);
  [qlon,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  fseek(grid_fid,3*4*(npad+IJDM),-1);
  [qlat,count]=fread(grid_fid,IJDM,'float32','ieee-be');

  fseek(grid_fid,4*4*(npad+IJDM),-1);
  [ulon,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  fseek(grid_fid,5*4*(npad+IJDM),-1);
  [ulat,count]=fread(grid_fid,IJDM,'float32','ieee-be');

  fseek(grid_fid,6*4*(npad+IJDM),-1);
  [vlon,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  fseek(grid_fid,7*4*(npad+IJDM),-1);
  [vlat,count]=fread(grid_fid,IJDM,'float32','ieee-be');

  fseek(grid_fid,8*4*(npad+IJDM),-1);
  [pang,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  
  fseek(grid_fid,9*4*(npad+IJDM),-1);
  [pscx,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  fseek(grid_fid,10*4*(npad+IJDM),-1);
  [pscy,count]=fread(grid_fid,IJDM,'float32','ieee-be');

  fseek(grid_fid,11*4*(npad+IJDM),-1);
  [qscx,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  fseek(grid_fid,12*4*(npad+IJDM),-1);
  [qscy,count]=fread(grid_fid,IJDM,'float32','ieee-be');

  fseek(grid_fid,13*4*(npad+IJDM),-1);
  [uscx,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  fseek(grid_fid,14*4*(npad+IJDM),-1);
  [uscy,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  
  fseek(grid_fid,15*4*(npad+IJDM),-1);
  [vscx,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  fseek(grid_fid,16*4*(npad+IJDM),-1);
  [vscy,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  
  fseek(grid_fid,17*4*(npad+IJDM),-1);
  [cori,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  
  fseek(grid_fid,18*4*(npad+IJDM),-1);
  [pasp,count]=fread(grid_fid,IJDM,'float32','ieee-be');
  
  
  plon=reshape(plon,IDM,JDM)';
  plat=reshape(plat,IDM,JDM)';

  qlon=reshape(qlon,IDM,JDM)';
  qlat=reshape(qlat,IDM,JDM)';

  ulon=reshape(ulon,IDM,JDM)';
  ulat=reshape(ulat,IDM,JDM)';

  vlon=reshape(vlon,IDM,JDM)';
  vlat=reshape(vlat,IDM,JDM)';

  pang=reshape(pang,IDM,JDM)';

  pscx=reshape(pscx,IDM,JDM)';
  pscy=reshape(pscy,IDM,JDM)';
  qscx=reshape(qscx,IDM,JDM)';
  qscy=reshape(qscy,IDM,JDM)';

  uscx=reshape(uscx,IDM,JDM)';
  uscy=reshape(uscy,IDM,JDM)';
  vscx=reshape(vscx,IDM,JDM)';
  vscy=reshape(vscy,IDM,JDM)';

  cori=reshape(cori,IDM,JDM)';

  pasp=reshape(pasp,IDM,JDM)';
  
  %% Close file
  fclose(grid_fid);
  
