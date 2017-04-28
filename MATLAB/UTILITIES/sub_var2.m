%
 function [aout]=sub_var2(fName,dims,num2,ivar,jj,ii)
%
%   Read/extract a 2-D variable from hycom archive file 
%----------------------------------------------------------------------
%   Input:
%          fName = [file Name]
%          dims  = dimensions jdm idm [kdm]
%          num2  = numer of 2d variable [check .b file]
%          ivar  = index of var (2d)
%
%   optional subset parameters 
%          jj =  j index 
%          ii =  i index 
%
%   Output:
%          aout
%----------------------------------------------------------------------
%   Xiaobiao Xu, 2007 
%
  aout = [];
  if  (nargin~=4 & nargin~=6)
      display('W r o n g  n a r g i n ');
  elseif (ivar>num2 | ivar<1);
      display('1 <= ivar <= num2 ');
  else
      fid = fopen(fName,'r','ieee-be');
      if (fid<0); display(['error in opening ' fName]);
      else

         jdm = dims(1);
         idm = dims(2);
         nn  = ceil(jdm*idm/4096)*4096;
%============================================================
         if nargin==6  % subregion 
%============================================================
            nj  = length(jj);
            ni  = length(ii);

            offs= (ivar-1);
            aout = nan*ones(nj,ni);
            for j=1:nj
                offset = (offs*nn + (jj(j)-1)*idm + ii(1)-1)*4;
                status = fseek(fid,offset,'bof');
                if status~=0; display('error in seek'); end
                aout(j,:) = fread(fid,[1,ni],'real*4');
            end
%============================================================
         else         % whole domain 
%============================================================
            offs = (ivar-1)*nn*4;
            status = fseek(fid,offs,'bof');
            if status~=0; display('error in seek'); end
            aone = fread(fid,[1,jdm*idm],'real*4');
            aout = reshape(aone,idm,jdm)';
            clear status aone
         end
         aout(aout>=2^100)=nan;
         fclose(fid);
      end
      clear fid

  end
%
