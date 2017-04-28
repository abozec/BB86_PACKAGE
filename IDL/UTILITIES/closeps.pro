;+
;
; @file_comments
; Close the Postscript mode
;
; when archive_ps ne 0, we add the name and the date at the bottom left corner
; of the postscript page.
; If the postscript is called idl.ps we change its name to number.ps
; (number automatically found to be 1 larger that any of the existing ps file)
;
; @keyword INFOWIDGET {type=long integer}
; id of the information widget (created by <pro>openps</pro>)
; that we have to destroy at the end of closeps (when the postscript is done).
;
; @uses
; <pro>cm_4ps</pro>
;
; @history
; Sebastien Masson (smasson\@lodyc.jussieu.fr)
;                       21/12/98
; June 2005: Sebastien Masson, english version with new commons
;
; @version
; $Id: closeps.pro 371 2008-08-07 09:32:02Z pinsard $
;
;-
PRO closeps, INFOWIDGET=infowidget
;
compile_opt idl2, strictarrsubs
;
IF lmgr(/demo) EQ 1 THEN return
;
;@cm_4ps
  IF NOT keyword_set(key_forgetold) THEN BEGIN
;@updatenew
  ENDIF
;
  IF !d.name NE 'PS' THEN GOTO, last_part
;------------------------------------------------------------
; if archive_ps /= 0 we will add its name and the date at the bottom
; left corner of the page (in case if the postscript will be archived
; in printps
;------------------------------------------------------------
   IF keyword_set(archive_ps) THEN BEGIN
;------------------------------------------------------------
; we get the name of the latest created postscript.
;------------------------------------------------------------
     psdir = isadirectory(psdir, title = 'Select psdir')
     nameps = file_search(psdir+'*.ps' $
                          , /test_regular, /test_write, /nosort)
     dates = (file_info(nameps)).mtime
     lastdate = (reverse(sort(temporary(dates))))[0]
     nameps = nameps[lastdate]
     nameps = file_basename(nameps, '.ps')
; If this name is idl.ps then we change it to the number.ps
     IF nameps EQ 'idl' then BEGIN
; get the name of all the *.ps or *.ps.gz files available in psdir
       allps = file_search(psdir+'*[.ps|.ps.gz|.pdf]', /test_regular, /nosort)
       allps = file_basename(file_basename(allps,'.gz'),'.ps')
       allps = file_basename(allps,'.pdf')
; find which of these names corresponds to numbers...
; get ascii codes of the names
       testnumb = byte(allps)
; longest name
       maxstrlen = (size(testnumb, /dimensions))[0]
; ascii codes can be 0 or between byte('0') and byte('9')
       testnumb = testnumb EQ 0 OR $
                  (testnumb GE (byte('0'))[0] AND testnumb LE (byte('9'))[0])
       testnumb = where(total(testnumb, 1) EQ maxstrlen, count)
       IF count NE 0 THEN BEGIN
; get the largest number
         psnumber = fix(allps[testnumb])
         psnumber = (psnumber[reverse(sort(psnumber))])[0] + 1
       ENDIF ELSE psnumber = 0
       nameps = strtrim(psnumber, 2)
     ENDIF
;------------------------------------------------------------
; we annotate the postscript
;------------------------------------------------------------
     date = byte(systime(0))    ; we get the date
     xyouts, !d.x_px_cm, !d.y_px_cm $
             , nameps+') '+string(date[4:10])+string(date[20:23]) $
             , /device, charsize = .75
   ENDIF
;------------------------------------------------------------
; close the postscript mode
   device, /close
;
last_part:
;
   thisOS = strupcase(strmid(!version.os_family, 0, 3))
   CASE thisOS of
     'MAC': SET_PLOT, thisOS
     'WIN': SET_PLOT, thisOS
     ELSE: SET_PLOT, 'X'
   ENDCASE
;    def_myuniquetmpdir
;    colorfile = myuniquetmpdir + 'original_colors.dat'
;    IF file_test(colorfile, /regular) THEN BEGIN
;      restore, colorfile
;      file_delete, colorfile, /quiet
; ; reload the original colors
;      tvlct, red, green, blue
;    ENDIF
   !p.font = -1
; force background color to the last color (white)
   ; !p.BACKGROUND=(!d.n_colors-1) < 255
   ; !p.color=0
   ; if !d.n_colors gt 256 then !p.background='ffffff'x
;------------------------------------------------------------
   if keyword_set(infowidget) then $
    widget_control, long(infowidget), bad_id = toto, /destroy
;------------------------------------------------------------
   return
end
