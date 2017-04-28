function write_relax_hycom(im, jm, km, io, file, field, d, res)

   %% Script to write relax files in the HYCOM format (initial conditions)  
   %% file is the name of the file WITHOUT the extension .a or .b 
   %% field='intf'
   %% field='temp'
   %% field='saln'
   %% d(kdm); target density
 
   switch field
      case 'intf'
         long_name = 'Interface Depths' ;
         short_name = 'int' ;
       
      case 'temp'
         long_name = 'Potential Temperature' ;
         short_name = 'tem' ;
       
      case 'saln' 
         long_name = 'Salinity' ; 
         short_name = 'sal' ;
    otherwise 
         disp('Should be intf, temp or saln');
   end  

   %% Dimensions of the domain
   kdm = km ;
   idm = im ;
   jdm = jm ;
   tdm = 12 ;
   ijdm = idm*jdm ;

   %% NPAD size
   npad=4096-mod(ijdm,4096);
   if (npad ~= 4096) 
     toto=zeros(npad,1);
   end

   %% Grid Directory and file 
   file1 = [io,file];

   %% WRITING .a binary file
   flga=[file1,'.a'];
   relax_fid=fopen(flga,'w');

   for t=1:tdm  
     for k = 1:kdm 
       A=res(:,:,k,t)';
       A=reshape(A,ijdm,1);
       %% Writing the field
       fwrite(relax_fid,A,'float32','ieee-be');
       %% Writing the padding at the end of the record:
       if (npad ~= 4096) 
	 fwrite(relax_fid,toto,'float32','ieee-be');
       end
     end 
   end
   fclose(relax_fid);

   %% mask the field by nan
   index = find(res > 1e10) ;
   res(index) = NaN ;

   
   %% WRITING .b text file
   density = d ;

   fldb=[file1,'.b'];
   fid1=fopen(fldb,'wt');
   fprintf( fid1, ' %s \n',long_name );
   fprintf( fid1, '\n');
   fprintf( fid1, '\n');
   fprintf( fid1, '\n');
   fprintf( fid1, 'i/jdm = %i %i \n', idm, jdm);   
   for m = 1:tdm 
      for k = 1:kdm 
         fprintf( fid1, ' %s : month,layer,dens,range =  %2.2i %2.2i %7.3f %10.5E  %10.5E\n', short_name, m, k, density(k), nanmin(nanmin(res(:, :, k,m)')), nanmax(nanmax(res(:,:, k,m)')));
         
      end 
   end 

   fclose(fid1);
