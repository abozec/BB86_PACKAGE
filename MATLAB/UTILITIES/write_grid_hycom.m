function write_grid_hycom(im, jm,  io, file,plon, plat, ulon,  ulat,  vlon, vlat, qlon, qlat, pang, ...
	  pscx, pscy, qscx,  qscy, uscx, uscy, vscx, vscy, cori, pasp, mapflg)

  %% Script to write the grid in the HYCOM format
  %% file is the name of the file WITHOUT the extension .a or .b 
  %% A. Bozec & D. Dukhovskoy Aug, 2011
    
  %% Get the id of the .a grid file
  flga=[io,file,'.a'];
  grid_fid=fopen(flga,'w');
  IDM=im;
  JDM=jm;
  IJDM=IDM*JDM;
  npad=4096-mod(IJDM,4096);
  if (npad ~= 4096) 
    toto=zeros(npad,1);
  end
  
  %% Writing .a file:
  A=plon';
  A=reshape(A,IJDM,1);
% Writing the field
  fwrite(grid_fid,A,'float32','ieee-be');
% Writing the padding at the end of the record:
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=plat';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=qlon';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=qlat';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=ulon';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=ulat';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=vlon';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=vlat';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  
  A=pang';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  
  A=pscx';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=pscy';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=qscx';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=qscy';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  
  A=uscx';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=uscy';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=vscx';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  A=vscy';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  
  A=cori';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

  
  A=pasp';
  A=reshape(A,IJDM,1);
  fwrite(grid_fid,A,'float32','ieee-be');
  if (npad ~= 4096) 
    fwrite(grid_fid,toto,'float32','ieee-be');
  end

    
  fclose(grid_fid);

  
  
  %% Writing .b file;
  fldb=[io,file,'.b'];
  fid1=fopen(fldb,'wt');
  fprintf(fid1,'  %i  ''idm   '' = longitudinal array size\n',IDM);
  fprintf(fid1,'  %i  ''jdm   '' = latitudinal array size\n', JDM);
  fprintf(fid1,'  %i  ''mapflg'' = map flag (0=mercator,10=panam,12=ulon-panam)\n',mapflg);
  fprintf(fid1,'plon:  min,max =     %12.5f  %12.5f\n', min(min(plon')), max(max(plon'))  );
  fprintf(fid1,'plat:  min,max =     %12.5f  %12.5f\n', min(min(plat')), max(max(plat'))  );
  fprintf(fid1,'qlon:  min,max =     %12.5f  %12.5f\n', min(min(qlon')), max(max(qlon'))  );
  fprintf(fid1,'qlat:  min,max =     %12.5f  %12.5f\n', min(min(qlat')), max(max(qlat'))  );
  fprintf(fid1,'ulon:  min,max =     %12.5f  %12.5f\n', min(min(ulon')), max(max(ulon'))  );
  fprintf(fid1,'ulat:  min,max =     %12.5f  %12.5f\n', min(min(ulat')), max(max(ulat'))  );
  fprintf(fid1,'vlon:  min,max =     %12.5f  %12.5f\n', min(min(vlon')), max(max(vlon'))  );
  fprintf(fid1,'vlat:  min,max =     %12.5f  %12.5f\n', min(min(vlat')), max(max(vlat'))  );
  fprintf(fid1,'pang:  min,max =     %12.5f  %12.5f\n', min(min(pang')), max(max(pang'))  );
  fprintf(fid1,'pscx:  min,max =     %12.5f  %12.5f\n', min(min(pscx')), max(max(pscx'))  );
  fprintf(fid1,'pscy:  min,max =     %12.5f  %12.5f\n', min(min(pscy')), max(max(pscy'))  );
  fprintf(fid1,'qscx:  min,max =     %12.5f  %12.5f\n', min(min(qscx')), max(max(qscx'))  );
  fprintf(fid1,'qscy:  min,max =     %12.5f  %12.5f\n', min(min(qscy')), max(max(qscy'))  );
  fprintf(fid1,'uscx:  min,max =     %12.5f  %12.5f\n', min(min(uscx')), max(max(uscx'))  );
  fprintf(fid1,'uscy:  min,max =     %12.5f  %12.5f\n', min(min(uscy')), max(max(uscy'))  );
  fprintf(fid1,'vscx:  min,max =     %12.5f  %12.5f\n', min(min(vscx')), max(max(vscx'))  );
  fprintf(fid1,'vscy:  min,max =     %12.5f  %12.5f\n', min(min(vscy')), max(max(vscy'))  );
  fprintf(fid1,'cori:  min,max =     %10.5E  %10.5E\n', min(min(cori')), max(max(cori'))  );
  fprintf(fid1,'pasp:  min,max =     %12.5f  %12.5f\n', min(min(pasp')), max(max(pasp'))  );
  
  fclose(fid1);
  
