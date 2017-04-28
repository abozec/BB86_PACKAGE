;+
;
; @file_comments
; switch to postscript mode and define it
;
; @param namepsin {in}{optional}
; name of the postscript file.
; Extension '.ps' is added if missing. It will be stored in the psdir directory.
;
; @keyword FILENAME
; to define the name of the postscript file through a keyword rather than with
; namepsin input argument (in this case the keyword can be pass through
; different routines via _EXTRA keyword).
;
; @keyword INFOWIDGET
; If INFOWIDGET is present, it specifies a named variable into which the id of
; the widget giving information about the postscript creation is stored as a
; long integer.
; This id is needed by <pro>closeps</pro> to kill the information widget.
;
; @keyword KEEP_PFONT
; activate to suppress the modification of !p.font
; (by default we force !p.font = 0  to make smaller postscripts)
;
; @keyword PORTRAIT
;
; @keyword LANDSCAPE
;
; @keyword KEEPPFONT
; same as keep_pfont
;
; @keyword LIGHTNESS
; a scalar used to change the Lightness of the color palette to be able to
; adjust according to the printer we use, the media (paper or slide)...
; > 1 to get darker colors
;
; @keyword _EXTRA
; Used to pass keywords to <proidl>DEVICE</proidl>.
;
; @uses
; <pro>cm_4ps</pro>
;
; @history
; Sebastien Masson (smasson\@lodyc.jussieu.fr)
; 21/12/98
; 1/2/98: ajout de nameps en input
; 1/9/1999: ajout du mot cle FILENAME et du widget
; June 2005: Sebastien Masson, cleaning, english version with new commons
;
; @version
; $Id: openps.pro 371 2008-08-07 09:32:02Z pinsard $
;
;-
PRO openps, namepsin, FILENAME=filename, PAGE_SIZE=page_size $
            , KEEPPFONT=keeppfont, KEEP_PFONT=keep_pfont $
            , PORTRAIT=key_portrait, LANDSCAPE=landscape $
            , LIGHTNESS=Lightness, _EXTRA=ex
;
;
  compile_opt idl2, strictarrsubs
;
   IF lmgr(/demo) EQ 1 THEN BEGIN
      dummy = report('impossible to create a PS in demo mode')
      return
   ENDIF
;
;@cm_4ps
IF NOT keyword_set(key_forgetold) THEN BEGIN
;@updatenew
ENDIF
;------------------------------------------------------------
; close the postscript device if we are already in postscript mode
   IF !d.name EQ 'PS' THEN device, /close
; switch to postscript mode
   set_plot,'ps'
;------------------------------------------------------------
; if we use  keyword Lightness
; save the actual color palette in a temporary file
; (to be restored when calling closeps
;------------------------------------------------------------
   IF n_elements(Lightness) NE 0 THEN BEGIN
     IF Lightness NE 1 THEN BEGIN
       tvlct, red, green, blue, /get
       def_myuniquetmpdir
       save, red, green, blue, filename = myuniquetmpdir + 'original_colors.dat'
       palit, Lightness, red, green, blue
     ENDIF
   ENDIF
;------------------------------------------------------------
; we define the name of the file
;------------------------------------------------------------
   CASE 1 OF
     n_params() EQ 1:nameps = namepsin
     keyword_set(filename): nameps = filename
     ELSE:nameps = xquestion('Name of the postscript file?', 'idl.ps', /chkwid)
   ENDCASE
; make sure that nameps ends with '.ps'
   nameps = file_dirname(nameps, /mark_directory) + $
            file_basename(nameps, '.ps') + '.ps'
; add path (psdir) and check that nameps is ok
;   nameps = isafile(nameps, iodir = psdir, /new)
;------------------------------------------------------------
; we define xsize, ysize, xoffset and yoffset
;------------------------------------------------------------
   IF n_elements(portrait) NE 0 OR n_elements(landscape) NE 0 THEN $
     key_portrait = keyword_set(portrait) * (1 - keyword_set(landscape))

   if key_portrait EQ 1 then begin
      xs = min(page_size)
      ys = max(page_size)
      xoff = 0.
      yoff = 0.
   ENDIF ELSE BEGIN
      xs = max(page_size)
      ys = min(page_size)
      xoff = 0.
      yoff = max(page_size)
   ENDELSE
;------------------------------------------------------------
; We define the device of the postscript mode
;------------------------------------------------------------
   device, /color, /Helvetica, filename = strcompress(nameps, /remove_all) $
           , LANDSCAPE = 1 - key_portrait, PORTRAIT = key_portrait $
           , xsize = xs, ysize = ys, xoffset = xoff, yoffset = yoff $
           , bits_per_pixel = 8, language_level = 2, _EXTRA = ex
; to make smaller postcripts
   IF NOT (keyword_set(keeppfont) OR keyword_set(keep_pfont)) $
   THEN !p.font = 0
   RETURN
END
