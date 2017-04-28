function write_depth_hycom(im, jm, file, bathy)

  %% Script to write the bathymetry in the HYCOM format
  %% file is the name of the file WITHOUT the extension .a or .b 
  %% A. Bozec & D. Dukhovskoy Aug, 2011

  %% Get the id of the bathymetry file
  flda=[file,'.a'];
  depth_fid=fopen(flda,'w');
  IDM=im;
  JDM=jm;
  IJDM=IDM*JDM;
  npad=4096-mod(IJDM,4096);
  toto=zeros(npad,1);

  %% Writing .a file:
  A=bathy';
  A=reshape(A,IJDM,1);
% Writing the bathy:
  fwrite(depth_fid,A,'float32','ieee-be');
% Writing the padding at the end of the record:
  fwrite(depth_fid,toto,'float32','ieee-be');
  fclose(depth_fid);

  %%find mask value
  ind=find(A > 1e10);
  A(ind)=NaN;

  %% Writing .b file;
  fldb=[file,'.b'];
  fid1=fopen(fldb,'wt');
  fprintf(fid1,'Bathymetry\n');
  fprintf(fid1,'i/jdm =  %i  %i \n',IDM,JDM);
  fprintf(fid1,'\n');
  fprintf(fid1,'\n');
  fprintf(fid1,'\n');
%  fprintf(fid1,'min,max depth =      %10.5f  %10.5f \n',nanmin(nanmin(A)), nanmax(nanmax(A)) );
  fprintf(fid1,'min,max depth =      %10.5f  %10.5f \n',min(A(find(~isnan(A)))), max(A(find(~isnan(A)))) );
  fclose(fid1);
