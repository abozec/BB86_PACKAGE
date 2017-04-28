%
 function [aout] = sub_var3(fName,dims,num2,num3,ivar,jj,ii,kk)
%
%   read/extract a 3-D variable from hycom archive file 
%----------------------------------------------------------------------
%   Input:
%        fName     : file name '.a'
%        dims      : dimension for full domain [jdm idm kdm] 
%        num2/3    : number of 2d/3d variables -- .b 
%        ivar      : index of the 3d variable  (1-num3)
%        [jj ii kk]: dimension for sub domain [optional] 
%----------------------------------------------------------------------
%   Xiaobiao Xu, 2007 
%

  aout = [];
  if  (nargin~=5 & nargin~=8)
      display('W r o n g  n a r g i n ');
  elseif (ivar>num3 | ivar<1);
      display('1 <= ivar <= num3 ');
  else

      fid = fopen(fName,'r','ieee-be');
      if  (fid<0); display(['error in opening ' fName]);
      else
          jdm = dims(1);
          idm = dims(2);
          kdm = dims(3);
          nn = ceil(jdm*idm/4096)*4096; 
% ========================================================================
          if  nargin==8 % subset of the variable --
% ========================================================================
              nj = length(jj);
              ni = length(ii);
              nk = length(kk);

              aout = nan(nj,ni,nk);
              for k=1:nk
                  offs = num2+num3*(kk(k)-1)+ivar-1;
                  for j=1:nj
                      offset = (offs*nn + (jj(j)-1)*idm + ii(1)-1)*4;
                      status = fseek(fid,offset,'bof');
                      if status~=0; display('error in seek'); end
                      aout(j,:,k) = fread(fid,[1,ni],'real*4');
                  end
              end
              clear aone k offset status
% ========================================================================
          else %full domain
% ========================================================================
              aout = nan(jdm,idm,kdm);
              for k=1:kdm
                  offset = (num2+num3*(k-1)+ivar-1)*nn*4;
                  status = fseek(fid,offset,'bof');
                  if status~=0; display('error in seek'); end
                  aone = fread(fid,[1,jdm*idm],'real*4');
                  aout(:,:,k) = reshape(aone,idm,jdm)';
              end
              clear aone k offset status
          end
          aout(aout>=2^100)=nan; 
          fclose(fid);
      end
      clear fid

  end
%
