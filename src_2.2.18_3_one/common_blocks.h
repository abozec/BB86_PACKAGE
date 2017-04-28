c-----------------------------------------------------------------------------
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2) ::
     & u,v,           ! velocity components
     & dp,dpu,dpv,    ! layer thickness
     & temp,          ! temperature
     & saln,          ! salinity
     & th3d,          ! potential density
     & thstar,        ! virtual potential density
     & montg          ! montgomery potential

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm,2,mxtrcr) ::
     & tracer         ! inert tracers

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
     & p,pu,pv        ! interface pressure

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     & dpold,         ! layer thickness
     & dpoldm,        ! layer thickness
     & theta,         ! isopycnal layer target densties - thbase
     & diaflx         ! time integral of diapyc.flux

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & corio,         ! coriolis parameter
     & potvor,        ! potential vorticity
     & srfhgt,        ! sea surface height, g*ssh(m)
     & steric,        ! steric sea surface height, g*sssh(m)
     & sshgmn,        !   mean sea surface height, g*mssh(m)
     & thmean,        !   mean depth averaged density
     & montg1,        ! layer 1 montgomery potential
     & skap           ! thermobaric scale factor between reference states

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     & psikk,         ! montg.pot. in bottom layer
     & thkk           ! virtual potential density in bottom layer


      common/hycom1r/ u,v,dp,dpold,dpoldm,dpu,dpv,p,pu,pv,
     &                corio,psikk,thkk,potvor,
     &                srfhgt,steric,sshgmn,thmean,montg1,
     &                temp,saln,th3d,thstar,skap,theta,diaflx,tracer
      save  /hycom1r/
c                                                                   
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & kapi           ! thermobaric reference state index (1 or 3)
      
      common/hycom1i/ kapi
      save  /hycom1i/
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) ::
     & uflx,vflx,     ! mass fluxes
     & uflxav,vflxav, ! average fluxes
     & dpav           ! average fluxes

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,3) ::
     & ubavg,vbavg,   ! barotropic velocity
     & pbavg          ! barotropic pressure

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & defor1,defor2, ! deformation components
     & ubrhs, vbrhs,  ! rhs of barotropic u,v eqns.
     & utotm, vtotm,  ! total (barotrop.+baroclin.)..
     & utotn, vtotn,  ! ..velocities at 2 time levels
     & uflux, vflux,  ! horizontal mass fluxes
     & uflux1,vflux1, ! more mass fluxes
     & uflux2,vflux2, ! more mass fluxes
     & uflux3,vflux3  ! more mass fluxes

      common/hycom2r/ montg,defor1,defor2,ubavg,vbavg,pbavg,
     &                ubrhs,vbrhs,utotm,vtotm,utotn,vtotn,
     &                uflux, vflux, uflux1,vflux1,
     &                uflux2,vflux2,uflux3,vflux3,
     &                uflx,vflx,uflxav,vflxav,dpav
      save  /hycom2r/
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & util1,util2,   ! arrays for temporary storage
     & util3,util4,   ! arrays for temporary storage
     & util5,util6,   ! arrays for temporary storage
     & plon, plat,    ! lon,lat at p pts
     & ulon, ulat,    ! lon,lat at u pts
     & vlon, vlat,    ! lon,lat at v pts
     & scux, scuy,    ! mesh size at u pts in x,y dir.
     & scvx, scvy,    ! mesh size at v pts in x,y dir.
     & scpx, scpy,    ! mesh size at p pts in x,y dir.
     & scqx, scqy,    ! mesh size at q pts in x,y dir.
     & scu2, scv2,    ! grid box area at u,v pts
     & scp2, scq2,    ! grid box area at p,q pts
     & scp2i,scq2i,   ! inverses of scp2,scq2
     & scuxi,scvyi,   ! inverses of scux,scvy
     & aspux,aspuy,   ! u-grid aspect ratios for diffusion
     & aspvx,aspvy,   ! v-grid aspect ratios for diffusion
     & veldf2u,       ! u-grid laplacian  diffusion coefficient
     & veldf2v,       ! v-grid laplacian  diffusion coefficient
     & veldf4u,       ! u-grid biharmonic diffusion coefficient
     & veldf4v,       ! v-grid biharmonic diffusion coefficient
     & thkdf4u,       ! u-grid biharmonic diffusion coefficient
     & thkdf4v,       ! v-grid biharmonic diffusion coefficient
     & pgfx, pgfy,    ! horiz. presssure gradient
     & gradx,grady,   ! horiz. presssure gradient
     & depthu,depthv, ! bottom pres. at u,v points
     & pvtrop,        ! pot.vort. of barotropic flow
     & depths,        ! water depth
     & drag,          ! bottom drag
     & dragrh,        ! tidal bottom drag roughness (r*H)
     & topiso,        ! shallowest depth for isopycnal layers (pressure units)
     & diwlat         ! spacially varying background/internal wave diffusivity

      common/hycom3r/ util1,util2,util3,util4,util5,util6,
     &                plon,plat,ulon,ulat,vlon,vlat,
     &                scux,scuy,scvx,scvy,scuxi,scvyi,
     &                scpx,scpy,scqx,scqy,
     &                scu2,scv2,scp2,scq2,scp2i,scq2i,
     &                aspux,aspuy,aspvx,aspvy,
     &                veldf2u,veldf2v,veldf4u,veldf4v,thkdf4u,thkdf4v,
     &                pgfx,pgfy,gradx,grady,
     &                depthu,depthv,pvtrop,depths,drag,dragrh,topiso,
     &                diwlat
      save  /hycom3r/
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & uja,   ujb,    ! velocities at lateral ..
     & via,   vib,    !       .. neighbor points
     & pbot,          ! bottom pressure at t=0
     & sgain,         ! salin.changes from diapyc.mix.
     & surtx,         ! surface net x-stress on p-grid
     & surty,         ! surface net y-stress on p-grid
     & surflx,        ! surface net thermal energy flux
     & sswflx,        ! surface swv thermal energy flux
     & mixflx,        ! mixed layer thermal energy flux
     & sstflx,        ! surface thermal  flux from sst relax
     & sssflx,        ! surface salinity flux from sss relax
     & salflx,        ! surface salinity flux
     & buoflx,        ! mixed layer buoyancy flux
     & bhtflx,        ! mixed layer buoyancy flux from heat
     & ustar,         ! friction velocity
     & ustarb,        ! bottom friction velocity
     & turgen,        ! turb.kin.energ. generation
     & thkice,        ! grid-cell avg. ice thknss (m)
     & covice,        ! ice coverage (rel.units)
     & temice,        ! ice surface temperature
     & flxice,        ! heat flux under ice
     & fswice,        ! swv  flux under ice
     & sflice,        ! salt flux under ice
     &   si_c,        ! ice concentration   on p-grid from coupler
     &   si_h,        ! ice thickness       on p-grid from coupler
     &   si_t,        ! ice temperature     on p-grid from coupler
     &   si_u,        ! ice u-velocity      on p-grid from coupler
     &   si_v,        ! ice v-velocity      on p-grid from coupler
     &  si_tx,        ! x-stress  under ice on p-grid from coupler
     &  si_ty         ! y-stesss  under ice on p-grid from coupler

      common/hycom4r/ uja,ujb,via,vib,pbot,
     &                sgain,surtx,surty,surflx,sswflx,mixflx,
     &                sstflx,sssflx,salflx,buoflx,bhtflx,
     &                ustar,ustarb,turgen,
     &                thkice,covice,temice,
     &                flxice,fswice,sflice,
     &                si_c,si_h,si_t,si_u,si_v,si_tx,si_ty
      save  /hycom4r/
c
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & klist,         ! k-index
     & jerlov         ! jerlov water type 1-5

      common/hycom4i/ klist,jerlov
      save  /hycom4i/
c
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     & dpmixl,        ! mixed layer depth
     & t1sav,         ! upper sublayer temperature
     & s1sav,         ! upper sublayer salinity
     & tmlb,          ! temp in lyr. containing mlb.
     & smlb           ! saln in lyr. containing mlb

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & hekman,        ! ekman layer thickness
     & hmonob,        ! monin-obukhov length
     & dpbl,          ! turbulent boundary layer depth
     & dpbbl,         ! bottom turbulent boundary layer depth
     & dpmold,        ! mixed layer depth at old time step
     & tmix,          ! mixed layer temperature
     & smix,          ! mixed layer salinity
     & thmix,         ! mixed layer potential density
     & umix,  vmix    ! mixed layer velocity

      real, dimension (kdm) ::
     & dp0k,          ! minimum deep    z-layer separation
     & ds0k,          ! minimum shallow z-layer separation
     & dssk           ! sigma depth scale factor

      common/hycom5r/ hekman,hmonob,dpbl,dpbbl,
     &                dpmold,dpmixl,tmix,smix,thmix,
     &                t1sav,s1sav,tmlb,smlb,umix,vmix,
     &                dp0k,ds0k,dssk
      save  /hycom5r/
c
      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     & nmlb           ! layer containing mlb.

      common/hycom5i/ nmlb
      save  /hycom5i/
c
c ---  s w i t c h e s    (if set to .true., then...)
c --- btrlfr      leapfrog barotropic time step
c --- btrmas      barotropic is mass conserving
c --- diagno      output model fields and diagnostic messages
c --- thermo      use thermodynamic forcing (flxflg>0)
c --- windf       use wind stress   forcing (wndflg>0)
c --- pcipf       use evap-precip surface salinity flux
c --- epmass      treat evap-precip as a mass exchange
c --- priver      use river precip bogas
c --- rivera      annual-only river precip bogas
c --- kparan      annual-only kpar
c --- relax       activate lateral boundary T/S/p  climatological nudging
c --- srelax      activate surface salinity        climatological nudging
c ---              (sssflg==1)
c --- trelax      activate surface temperature     climatological nudging
c ---              (sstflg==1)
c --- trcrlx      activate lateral boundary tracer climatological nudging
c --- relaxf      input T/S/p   relaxation fields
c --- relaxs      input surface relaxation fields only
c --- relaxt      input tracer  relaxation fields
c --- locsig      use locally-referenced potential density for stability
c --- vsigma      use spacially varying target densities
c --- hybrid      use hybrid vertical coordinates
c --- isopyc      use isopycnic vertical coordinates (MICOM mode)
c --- icegln      use energy loan ice model (iceflg==1)
c --- isopcm      HYBGEN: use PCM to remap isopycnal layers
c --- mxlkta      KT:    activate    original mixed layer model (mlflag==2)
c --- mxlktb      KT:    activate alternative mixed layer model (mlflag==3)
c --- mxlkrt      KT:    activate MICOM or HYCOM Kraus-Turner (mlflag==2,3)
c --- pensol      KT:    activate penetrating solar radiation
c --- mxlkpp      KPP:   activate mixed layer model (mlflag==1)
c --- bblkpp      KPP:   activate bottom boundary layer
c --- shinst      KPP:   activate shear instability mixing
c --- dbdiff      KPP:   activate double diffusion  mixing
c --- nonloc      KPP:   activate nonlocal b. layer mixing
c --- latdiw      KPROF: activate lat.depen.int.wav mixing
c --- botdiw      GISS:  activate bot.enhan.int.wav mixing
c --- difout      KPROF: output visc/diff coeffs in archive
c --- mxlmy       MY2.5: activate mixed layer model (mlflag==5)
c --- mxlpwp      PWP:   activate mixed layer model (mlflag==4)
c --- mxlgiss     GISS:  activate mixed layer model (mlflag==6)
c --- flxoff      add a net heat flux offset
c --- flxsmo      activate smoothing of surface fluxes
c --- trcrin      initialize tracer from restart file
c --- trcout      advect tracer and save results in history/restart file
c --- dsur1p      single point only surface diagnostics
c --- arcend      always write a 3-d archive at the end of the run
c
      logical       btrlfr,btrmas,diagno,thermo,windf,
     &              pcipf,epmass,priver,rivera,kparan,
     &              relax,srelax,trelax,trcrlx,relaxf,relaxs,relaxt,
     &              locsig,vsigma,hybrid,isopyc,icegln,isopcm,
     &              mxlkta,mxlktb,mxlkrt,pensol,
     &              mxlkpp,bblkpp,shinst,dbdiff,nonloc,
     &              latdiw,botdiw,difout,
     &              mxlmy,mxlpwp,mxlgiss,flxoff,flxsmo,trcrin,trcout,
     &              dsur1p,arcend
      common/swtchs/btrlfr,btrmas,diagno,thermo,windf,
     &              pcipf,epmass,priver,rivera,kparan,
     &              relax,srelax,trelax,trcrlx,relaxf,relaxs,relaxt,
     &              locsig,vsigma,hybrid,isopyc,icegln,isopcm,
     &              mxlkta,mxlktb,mxlkrt,pensol,
     &              mxlkpp,bblkpp,shinst,dbdiff,nonloc,
     &              latdiw,botdiw,difout,
     &              mxlmy,mxlpwp,mxlgiss,flxoff,flxsmo,trcrin,trcout,
     &              dsur1p,arcend
      save  /swtchs/
c
c ---  t e x t
c ---  ctitle     four lines describing the simulation
c
      character*80    ctitle
      common/hycom1c/ ctitle(4)
      save  /hycom1c/
c
c --- atmospheric forcing fields
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) ::
     & taux,          ! wind stress in x direction
     & tauy,          ! wind stress in y direction
     & wndspd,        ! wind speed
     & ustara,        ! ustar (tke source)
     & airtmp,        ! air temperature
     & vapmix,        ! atmosph. vapor mixing ratio
     & precip,        ! precipitation
     & radflx,        ! net solar radiation
     & swflx,         ! net shortwave radiation
     & surtmp,        ! surface temp. used to calculate input lw radiation
     & seatmp,        ! best available SST from observations
     & akpar,         ! photosynthetically available radiation coefficent
     & rivers         ! river inflow bogused to surface precipitation

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & offlux         ! net heat flux offset

      real, dimension (5) ::
     & betard,        ! red  extinction coefficient
     & betabl,        ! blue extinction coefficient
     & redfac         ! fract. of penetr. red light

      real, dimension (0:24,-91:91) ::
     & diurnl         ! hourly vs latitude shortwave scale factor table

      common/frcing/ taux,tauy,wndspd,ustara,
     &               airtmp,vapmix,precip,radflx,swflx,surtmp,seatmp,
     &               akpar,rivers,
     &               offlux,
     &               betard,betabl,redfac,
     &               diurnl
      save  /frcing/
c
c --- surface and sidewall and nestwall boundary fields
c ---  (kkwall and kknest are either kdm or, if inactive, 1).
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkwall,4) ::
     & pwall,         ! pressure b.c. at sidewalls
     & swall,         ! salinity b.c. at sidewalls
     & twall          ! temp.    b.c. at sidewalls

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kkwall,4,
     &                                                   mxtrcr) ::
     & trwall         ! tracer   b.c. at sidewalls

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kknest,2) ::
     & pnest,         ! pressure b.c. at nestwalls
     & snest,         ! salinity b.c. at nestwalls
     & tnest,         ! temp.    b.c. at nestwalls
     & unest,         ! u-vel.   b.c. at nestwalls
     & vnest          ! v-vel.   b.c. at nestwalls
 
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,2) ::
     & ubnest,        ! barotropic u-velocity at nestwalls
     & vbnest,        ! barotropic v-velocity at nestwalls
     & ubpnst,        ! barotropic u-velocity at nestwalls on p-grid
     & vbpnst,        ! barotropic v-velocity at nestwalls on p-grid
     & pbnest         ! barotropic pressure   at nestwalls

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & rmu,           ! weights for   s.w.b.c. relax
     & rmunp,         ! weights for p.n.w.b.c. relax
     & rmunv,         ! weights for v.n.w.b.c. relax
     & rmutra         ! weights for tracr.b.c. relax (maximum of all tracers)

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,mxtrcr) ::
     & rmutr          ! weights for tracr.b.c. relax

      common/wall1r/ rmu,  pwall,swall,twall,
     &               rmutra,rmutr,trwall,
     &               rmunp,rmunv, pnest,snest,tnest,unest,vnest,
     &               ubnest,vbnest,ubpnst,vbpnst,pbnest
      save  /wall1r/

      integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     & maskbc         ! mask for nested barotropic boundary condition

      common/wall1i/ maskbc
      save  /wall1i/
c
c --- pwp variables
      real ::
     & rigc           ! PWP: critical gradient richardson number
     &,ribc           ! PWP: critical bulk richardson number

      common/pwpr/ rigc,ribc
      save  /pwpr/
c
c --- m-y 2.5 variables
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1,2) ::
     & q2             !  tke
     &,q2l            !  tke * turbulent length scale

      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,0:kkmy25+1) ::
     & difqmy         !  tke diffusivity
     &,vctymy         !  viscosity on mellor-yamada vertical grid
     &,diftmy         !  temperature diffusivity on mellor-yamada vertical grid

      real ::
     & ghc            !  constant for calculating tke production
     &,sef            !  constant for calculating tke production
     &,smll           !  constant for calculating tke
     &,const1         !  coefficient for estimating surface and bottom bc's
     &,coef4          !  coefficient for calculating  viscosity/diffusivity
     &,coef5          !  coefficient for calculating  viscosity/diffusivity
     &,a1my           !  coefficient for calculating  viscosity/diffusivity
     &,b1my           !  coefficient for calculating  viscosity/diffusivity
     &,a2my           !  coefficient for calculating  viscosity/diffusivity
     &,b2my           !  coefficient for calculating  viscosity/diffusivity
     &,c1my           !  coefficient for calculating  viscosity/diffusivity
     &,e1my           !  coefficient for calculating  viscosity/diffusivity
     &,e2my           !  coefficient for calculating  viscosity/diffusivity
     &,e3my           !  coefficient for calculating  viscosity/diffusivity
c
      common/myr/ q2,q2l,difqmy,vctymy,diftmy,
     &            ghc,sef,smll,const1,
     &            coef4,coef5,a1my,b1my,a2my,b2my,c1my,e1my,e2my,e3my
      save  /myr/
c
c --- kpp variables
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm+1) ::
     & zgrid          !  grid levels in meters
     &,vcty           !  vert. viscosity coefficient
     &,difs           !  vert. scalar diffusivity
     &,dift           !  vert. temperature diffusivity
     &,ghats          !  nonlocal transport

      real ::
     & vonk           !  von karman constant
     &,zmin,zmax      !  zehat limits for table
     &,umin,umax      !  ustar limits for table
     &,epsilon        ! vertical coordinate scale factor
     &,cmonob         ! KPP: scale factor for Monin-Obukov length
     &,cekman         ! KPP: scale factor for Ekman depth
     &,qrinfy         ! KPP: 1/max grad.rich.no. for shear instability
     &,difm0          ! KPP: max viscosity    due to shear instability
     &,difs0          ! KPP: max diffusivity  due to shear instability
     &,difmiw         ! KPP: background/internal wave viscosity   (m^2/s)
     &,difsiw         ! KPP: background/internal wave diffusivity (m^2/s)
     &,qdif0          ! KPP: difm0 /difs0
     &,qdifiw         ! KPP: difmiw/difsiw
     &,dsfmax         ! KPP: salt fingering diffusivity factor    (m^2/s)
     &,rrho0          ! KPP: salt fingering rp=(alpha*delT)/(beta*delS)
     &,ricr           ! KPP: critical bulk richardson number
     &,cs             ! KPP: value for nonlocal flux term
     &,cstar          ! KPP: value for nonlocal flux term
     &,cv             ! KPP: buoyancy frequency ratio (0.0 to use a fn. of N)
     &,c11            ! KPP: value for turb velocity scale
     &,deltaz         ! delta zehat in table
     &,deltau         ! delta ustar in table
     &,vtc            ! constant for estimating background shear in rib calc.
     &,cg             ! constant for estimating nonlocal flux term of diff. eq.
     &,dp0enh         ! dist. for tapering diff. enhancement at interface nbl-1

      common/kppr/ zgrid,vcty,dift,difs,ghats,
     &             vonk,zmin,zmax,umin,umax,
     &             qrinfy,difm0,difs0,difmiw,difsiw,qdif0,qdifiw,
     &             rrho0,dsfmax,ricr,epsilon,
     &             cmonob,cekman,cs,cstar,cv,c11,
     &             deltaz,deltau,vtc,cg,dp0enh
      save  /kppr/
c
      integer ::
     & niter          ! KPP: iterations for semi-implicit soln. (2 recomended)
     &,hblflg         ! KPP: b.layer interpolation flag (0=con.1=lin.,2=quad.)
      common/kppi/ niter,hblflg
      save  /kppi/
c
c --- nasa giss variables
c
      integer nextrtbl0,ifexpabstable,nextrtbl1,
     &        nextrtbl,nposapprox,mt0,mt,ntbl,
     &        mt_ra_r,n_theta_r_oct,nbig

      common/gissi1/
     &        nextrtbl0,ifexpabstable,nextrtbl1,
     &        nextrtbl,nposapprox,mt0,mt,ntbl,
     &        mt_ra_r,n_theta_r_oct,nbig
      save  /gissi1/
c
      real    deltheta_r,pidbl,rri

      common/gissr1/
     &        deltheta_r,pidbl,rri
      save  /gissr1/
c
      integer
     & irimax(-762:762)
     &,nb

      common/gissi2/
     &        irimax,nb
      save  /gissi2/
c
      real
     & ribtbl(-762:762)
     &,ridb(  -762:762)
     &,slq2b( -762:762,-nlgiss:nlgiss)
     &,dri
     &,smb(   -762:762,-nlgiss:nlgiss)
     &,shb(   -762:762,-nlgiss:nlgiss)
     &,ssb(   -762:762,-nlgiss:nlgiss)
     &,back_ra_r(-39:117)
     &,sisamax(  -39:117)
     &,ra_rmax(  -39:117)
     &,c_y_r0(   -39:117)
     &,sm_r1(    -39:117)
     &,sh_r1(    -39:117)
     &,ss_r1(    -39:117)
     &,slq2_r1(  -39:117)
     &,b1,visc_cbu_limit,diff_cbt_limit
     &,theta_rcrp,theta_rcrn

      common/gissr2/
     &        ribtbl,ridb,slq2b,dri,smb,shb,ssb,back_ra_r,
     &        sisamax,ra_rmax,c_y_r0,sm_r1,sh_r1,ss_r1,slq2_r1,
     &        b1,visc_cbu_limit,diff_cbt_limit,
     &        theta_rcrp,theta_rcrn
      save  /gissr2/
c
      integer ifback,ifsali,ifepson2,ifrafgmax,
     &        ifsalback,ifchengcon,ifunreal,idefmld,
     &        ifpolartablewrite,ifbg_theta_interp

      common/gissi3/
     &        ifback,ifsali,ifepson2,ifrafgmax,
     &        ifsalback,ifchengcon,ifunreal,idefmld,
     &        ifpolartablewrite,ifbg_theta_interp
      save  /gissi3/
c
      real   back_ph_0,adjust_gargett,back_k_0,back_del_0,back_s2,
     &       ri0,ebase,epson2_ref,
     &       eps_bot0,scale_bot,      !for bottom-enhanced
     &       eplatidepmin,wave_30,    !and latitude dependent mixing
     &       deltemld,delrhmld,
     &       back_sm2,v_back0,t_back0,
     &       s_back0,ri_internal,backfrac,backfact,ako,tpvot0,sgmt,
     &       tptot0,tpcot0,ttot0,tcot0,tctot0,tpvot,tptot,tpcot,
     &       ttot,tcot,tctot,back_l_0

      common/gissr3/
     &       back_ph_0,adjust_gargett,back_k_0,back_del_0,back_s2,
     &       ri0,ebase,epson2_ref,
     &       eps_bot0,scale_bot,      !for bottom-enhanced
     &       eplatidepmin,wave_30,    !and latitude dependent mixing
     &       deltemld,delrhmld,
     &       back_sm2,v_back0,t_back0,
     &       s_back0,ri_internal,backfrac,backfact,ako,tpvot0,sgmt,
     &       tptot0,tpcot0,ttot0,tcot0,tctot0,tpvot,tptot,tpcot,
     &       ttot,tcot,tctot,back_l_0
      save  /gissr3/
c
      real*8           area,avgbot,watcum,empcum
      common/varblsd/  area,avgbot,watcum,empcum
      save  /varblsd/
c
      real            time,delt1,dlt,
     &                w0, w1, w2, w3,  ! wind  interp. scale factors
     &                wk0,wk1,wk2,wk3, ! kpar  interp. scale factors
     &                wr0,wr1,wr2,wr3, ! river interp. scale factors
     &                wc0,wc1,wc2,wc3, ! clim. interp. scale factors
     &                wn0,wn1,         ! nest  interp. scale factors
     &                wb0,wb1          ! baro. interp. scale factors
      common/varblsr/ time,delt1,dlt,
     &                w0, w1, w2, w3, wk0,wk1,wk2,wk3,
     &                wr0,wr1,wr2,wr3,wc0,wc1,wc2,wc3,
     &                wn0,wn1,wb0,wb1
      save  /varblsr/
c
      integer         nstep,nstep1,nstep2,lstep,
     &                l0, l1, l2, l3,  ! wind  indexes
     &                lk0,lk1,lk2,lk3, ! kpar  indexes
     &                lr0,lr1,lr2,lr3, ! river indexes
     &                lc0,lc1,lc2,lc3, ! clim. indexes
     &                ln0,ln1,         ! nest  indexes
     &                lb0,lb1          ! baro. indexes
      common/varblsi/ nstep,nstep1,nstep2,lstep,
     &                l0, l1, l2, l3, lk0,lk1,lk2,lk3,
     &                lr0,lr1,lr2,lr3,lc0,lc1,lc2,lc3,
     &                ln0,ln1,lb0,lb1
      save  /varblsi/
c
c --- 'sigma ' = isopyncnal layer target densities (sigma units)
c --- 'thbase' = reference density (sigma units)
c --- 'saln0'  = initial salinity value
c --- 'baclin' = baroclinic time step
c --- 'batrop' = barotropic time step
c ---'qhybrlx' = HYBGEN: relaxation coefficient (inverse baroclinic time steps)
c --- 'hybiso' = HYBGEN: Use PCM if layer is within hybiso of target density
c --- 'visco2' = deformation-dependent Laplacian  viscosity factor
c --- 'visco4' = deformation-dependent biharmonic viscosity factor
c --- 'facdf4' =       speed-dependent biharmonic viscosity factor
c --- 'veldf2' = diffusion velocity (m/s) for Laplacian  momentum dissipation
c --- 'veldf4' = diffusion velocity (m/s) for biharmonic momentum dissipation
c --- 'temdf2' = diffusion velocity (m/s) for Laplacian  temp/saln diffusion
c --- 'temdfc' = temp diffusion conservation (0.0 all density, 1.0 all temp)
c --- 'thkdf2' = diffusion velocity (m/s) for Laplacian  thickness diffusion
c --- 'thkdf4' = diffusion velocity (m/s) for biharmonic thickness diffusion
c --- 'vertmx' = diffusion velocity (m/s) for mom.mixing across mix.layr.base
c --- 'tofset' = temperature anti-drift offset (degC/century)
c --- 'sofset' = salnity     anti-drift offset  (psu/century)
c --- 'diapyc' = KT: diapycnal diffusivity x buoyancy freq. (m**2/s**2)
c --- 'dtrate' = KT: maximum permitted m.l. detrainment rate (m/day)
c --- 'slip'   = +1 for free-slip, -1  for non-slip boundary conditions
c --- 'cb'     = coefficient of quadratic bottom friction
c --- 'cbar'   = rms flow speed (m/s) for linear bottom friction law
!!Alex add linear bottom drag cbar2
c --- 'cbar2'  = linear bottom drag
c --- 'drglim' = limiter for explicit friction (1.0 no limiter, 0.0 implicit)
c --- 'drgscl' = scale factor for tidal drag   (0.0 for no tidal drag)
c --- 'thkdrg' = thickness of bottom boundary layer for tidal drag (m)
c --- 'dsurfq' = number of days between model diagnostics at the surface
c --- 'diagfq' = number of days between model diagnostics
c --- 'tilefq' = number of days between model diagnostics on some tiles
c --- 'meanfq' = number of days between model diagnostics (time averaged)
c --- 'rstrfq' = number of days between model restart output
c --- 'bnstfq' = number of days between baro. nesting archive input
c --- 'nestfq' = number of days between 3-d   nesting archive input
c --- 'wuv1/2' = weights for time smoothing of u,v field
c --- 'wts1/2' = weights for time smoothing of t,s field
c --- 'wbaro'  = weight for time smoothing of barotropic u,v,p field
c --- 'thkmls' = reference mixed-layer thickness for SSS relaxation (m)
c --- 'thkmlt' = reference mixed-layer thickness for SST relaxation (m)
c --- 'thkriv' = nominal thickness of river inflow (m)
c --- 'thkfrz' = maximum thickness of near surface freezing zone (m)
c --- 'tfrz_0' = ENLN: ice melting point (degC) at S=0psu
c --- 'tfrz_s' = ENLN: gradient of ice melting point (degC/psu)
c --- 'ticegr' = ENLN: vertical temperature gradient inside ice (deg/m)
c ---                   (0.0 to get ice surface temp. from atmos. surtmp)
c --- 'hicemn' = ENLN: minimum ice thickness (m)
c --- 'hicemx' = ENLN: maximum ice thickness (m)
c --- 'thkmin' = KT/PWP: minimum mixed-layer thickness (m)
c --- 'bldmin' = KPP:    minimum surface boundary layer thickness (m)
c --- 'bldmax' = KPP:    maximum surface boundary layer thickness (m)
c --- 'thkbot' = thickness of bottom boundary layer (m)
c --- 'sigjmp' = minimum density jump across interfaces   (theta units)
c --- 'tmljmp' = equivalent temperature jump across the mixed layer (degC)
c --- 'salmin' = minimum salinity allowed in an isopycnic layer
c --- 'dp00'   = deep    z-level spacing minimum thickness (m)
c --- 'dp00x'  = deep    z-level spacing maximum thickness (m)
c --- 'dp00f'  = deep    z-level spacing stretching factor (1.0=const.z)
c --- 'ds00'   = shallow z-level spacing minimum thickness (m)
c --- 'ds00x'  = shallow z-level spacing maximum thickness (m)
c --- 'ds00f'  = shallow z-level spacing stretching factor (1.0=const.z)
c --- 'dp00i'  = deep iso-pycnal spacing minimum thickness (m)
c --- 'isotop' = shallowest depth for isopycnal layers     (m), <0 from file
c --- 'nhybrd' = number of hybrid levels (0=all isopycnal)
c --- 'nsigma' = number of sigma  levels (nhybrd-nsigma z-levels)
c --- 'hybmap' = HYBGEN:  remapper  flag (0=PCM,1=PLM,2=PPM,-ve:isoPCM)
c --- 'hybflg' = HYBGEN:  generator flag (0=T&S,1=th&S,2=th&T)
c --- 'advflg' = thermal  advection flag (0=T&S,1=th&S,2=th&T)
c --- 'advtyp' = scalar   advection type (0=PCM,1=MPDATA,2=FCT2,4=FCT4)
c --- 'momtyp' = momentum advection type (2=2nd order, 4=4th order)
c --- 'kapref' = thermobaric reference state (-1=input,0=none,1,2,3=constant)
c --- 'kapnum' = number of thermobaric reference states (1 or 2)
c --- 'tsofrq' = number of time steps between anti-drift offset calcs
c --- 'mixfrq' = KT: number of time steps between diapycnal mixing calcs
c --- 'icefrq' = number of time steps between sea-ice updates
c --- 'ntracr' = number of tracers (<=mxtrcr)
c --- 'trcflg' = tracer type flag (one per tracer)
c --- 'clmflg' = climatology frequency flag (6=bimonthly,12=monthly)
c --- 'dypflg' = KT: diapycnal mixing flag (0=none,1=KPP,2=explicit)
c --- 'iniflg' = initial state flag (0=level,1=zonal,2=climatology)
c --- 'lbflag' = lateral barotropic bndy flag (0=none,1=port,2=nest,3=flather)
c --- 'mapflg' = map flag (0=mercator,2=uniform,3=beta-plane,4=input)
c --- 'yrflag' = days in year flag   (0=360,1=366,2=366Jan1,3=actual)
c --- 'sshflg' = diagnostic SSH flag (0=SSH,1=SSH&stericSSH)
c --- 'iversn' = hycom version number x10
c --- 'iexpt'  = experiment number x10
c --- 'jerlv0' = initial jerlov water type (1 to 5; 0 to use kpar)
c --- 'iceflg' = sea ice model flag   (0=none,1=energy loan,2=coupled/esmf)
c --- 'wndflg' = wind str. input flag (0=none,1=on u/v grid,2,3=on p grid)
c --- 'ustflg' = ustar forcing flag          (3=input,1=wndspd,2=stress)
c --- 'flxflg' = thermal forcing flag (0=none,3=net-flux,1,2,4=sst-based)
c --- 'empflg' = E-P     forcing flag (0=none,3=net_E-P, 1,2,4=sst-based_E)
c --- 'dswflg' = diurnal shortwv flag (0=none,1=daily to diurnal correction)
c --- 'sssflg' = SSS relaxation  flag (0=none,1=clim)
c --- 'lwflag' = longwave corr.  flag (0=none,1=clim,2=atmos), sst-based
c --- 'sstflg' = SST relaxation  flag (0=none,1=clim,2=atmos,3=obs)
c --- 'icmflg' = ice mask        flag (0=none,1=clim,2=atmos,3=obs)
c --- 'difsmo' = KPROF: number of layers with horiz smooth diff coeffs
c
!!Alex cbar2
      real           sigma,thbase,saln0,baclin,batrop,
     &               qhybrlx,hybiso,
     &               visco2,visco4,veldf2,veldf4,facdf4,
     &               temdf2,temdfc,thkdf2,thkdf4,vertmx,diapyc,
     &               tofset,sofset,dtrate,slip,cb,cbar,
     &               drglim,drgscl,thkdrg,
     &               dsurfq,diagfq,tilefq,meanfq,rstrfq,bnstfq,nestfq,
     &               wuv1,wuv2,wts1,wts2,wbaro,
     &               thkmls,thkmlt,thkriv,thkmin,bldmin,bldmax,thkbot,
     &               thkfrz,tfrz_0,tfrz_s,ticegr,hicemn,hicemx,
     &               dp00,dp00f,dp00x,ds00,ds00f,ds00x,dp00i,isotop,
     &               sigjmp,tmljmp,salmin,
     &               cbar2
!!Alex cbar2
      common/parms1r/sigma(kdm),thbase,saln0,baclin,batrop,
     &               qhybrlx,hybiso,
     &               visco2,visco4,veldf2,veldf4,facdf4,
     &               temdf2,temdfc,thkdf2,thkdf4,vertmx,diapyc,
     &               tofset,sofset,dtrate,slip,cb,cbar,
     &               drglim,drgscl,thkdrg,
     &               dsurfq,diagfq,tilefq,meanfq,rstrfq,bnstfq,nestfq,
     &               wuv1,wuv2,wts1,wts2,wbaro,
     &               thkmls,thkmlt,thkriv,thkmin,bldmin,bldmax,thkbot,
     &               thkfrz,tfrz_0,tfrz_s,ticegr,hicemn,hicemx,
     &               dp00,dp00f,dp00x,ds00,ds00f,ds00x,dp00i,isotop,
     &               sigjmp,tmljmp,salmin(kdm),
     &               cbar2

      save  /parms1r/
c
      integer        tsofrq,mixfrq,icefrq,nhybrd,nsigma,
     &               hybmap,hybflg,advflg,advtyp,momtyp,
     &               kapref,kapnum,
     &               ntracr,trcflg(mxtrcr),
     &               clmflg,dypflg,iniflg,lbflag,mapflg,yrflag,sshflg,
     &               iversn,iexpt,jerlv0,
     &               iceflg,icmflg,wndflg,ustflg,
     &               flxflg,empflg,dswflg,lwflag,sstflg,sssflg,
     &               difsmo
      common/parms1i/tsofrq,mixfrq,icefrq,nhybrd,nsigma,
     &               hybmap,hybflg,advflg,advtyp,momtyp,
     &               kapref,kapnum,
     &               ntracr,trcflg,
     &               clmflg,dypflg,iniflg,lbflag,mapflg,yrflag,sshflg,
     &               iversn,iexpt,jerlv0,
     &               iceflg,icmflg,wndflg,ustflg,
     &               flxflg,empflg,dswflg,lwflag,sstflg,sssflg,
     &               difsmo
      save  /parms1i/
c
c --- 'tenm,onem,...' = pressure thickness values corresponding to 10m,1m,...
c --- 'qonem'  = 1/onem
c --- 'g'      = gravity acceleration
c --- 'thref'  = reference value of specific volume (m**3/kg)
c --- 'qthref' = 1/thref
c --- 'spcifh' = specific heat of sea water (j/kg/deg)
c --- 'epsil'  = small nonzero number used to prevent division by zero
c --- 'huge'   = large number used to indicate land points
c
      real          tenm,onem,tencm,onecm,onemm,qonem,
     &              g,thref,qthref,spcifh,epsil,huge,radian,pi
      common/consts/tenm,onem,tencm,onecm,onemm,qonem,
     &              g,thref,qthref,spcifh,epsil,huge,radian,pi
      save  /consts/
c
c --- grid point where detailed diagnostics are desired:
      integer       itest,jtest,ittest,jttest
      common/testpt/itest,jtest,ittest,jttest
      save  /testpt/
c
c --- filenames.
      character*80  flnmdep,flnmgrd,flnmrsi,flnmrso, flnmflx,
     &              flnmarc,flnmovr,flnmfor,flnmforw,flnminp,
     &              flnmarcm,
     &              flnmarct
      common/iovars/flnmdep,flnmgrd,flnmrsi,flnmrso, flnmflx,
     &              flnmarc,flnmovr,flnmfor,flnmforw,flnminp,
     &              flnmarcm,
     &              flnmarct
      save  /iovars/
c
c --- CCSM3 variables (should be in CCSM3 modules)
      logical       dosstbud,doovtn,chk_ovtn,dodump,dohist,dorestart
      common/ccsmhl/dosstbud,doovtn,chk_ovtn,dodump,dohist,dorestart
      save  /ccsmhl/
c
      integer       istrt_mo,icurrent_mo,istrt_yr,icurrent_yr
      common/ccsmhi/istrt_mo,icurrent_mo,istrt_yr,icurrent_yr
      save  /ccsmhi/
c --- diurnal cycle factor for short wave heat flux
      integer       nsteps_per_day,nsteps_today
      common/dicyci/nsteps_per_day,nsteps_today
      save  /dicyci/
c
      real, dimension (nsteps_baclin) ::
     &              diurnal_cycle_factor
      common/dicycr/diurnal_cycle_factor
      save  /dicycr/

!!Alex add pstrsi for BB86 config
c --- 'pstrsi' = depth over which the wind is apply (m) (> 0 for bb86 only )
      real          pstrsi  
      common/bb86/pstrsi
      save  /bb86/

c
c> Revision history:
c>
c> Feb. 2001 - added halo and converted to f90 declarations
c> Aug. 2001 - added bnstfq,nestfq
c> Jan. 2002 - added curvilinear grid arrays and deleted /pivot/
c> May  2002 - added d[ps]00*
c> May  2002 - added PWP and MY2.5 mixed layer models
c> Aug  2002 - added ntracr and trcflg
c> Nov  2002 - added thkmls and thkmlt
c> Apr  2003 - added dp00i, vsigma, and priver
c> May  2003 - added bldmin, bldmax, flxsmo, and icmflg
c> Jun  2003 - added locsig
c> Nov  2003 - added advtyp
c> Jan  2004 - added latdiw
c> Jan  2004 - added bblkpp
c> Jan  2004 - added hblflg
c> Feb  2004 - added botdiw
c> Feb  2004 - added temdfc
c> Mar  2004 - added thkriv and epmass
c> Mar  2004 - added isotop
c> Mar  2005 - added tfrz_0, tfrz_s, ticegr, hicemn, and hicemx
c> Mar  2005 - added tsofrq, tofset, and sofset
c> Mar  2005 - added empflg
c> May  2005 - added kapref and kapnum
c> Jun  2006 - added surtx,surty
c> Jun  2006 - added icefrq,txice,tyice,uice,vice,flxice,fswice,sflice
c> Jun  2006 - added thkfrz,
c> Jan  2007 - added si_t; renamed si_[chuv] and si_t[xy]
c> Feb  2007 - added CCSM3-only variables
c> Apr  2007 - added dragrh,drglim,drgscl
c> Apr  2007 - added srfhgt,montg1
c> Apr  2007 - added btrlfr and btrmas
c> Jun  2007 - added momtyp and facdf4.
c> Sep  2007 - added hybmap, hybiso and isopcm.
c> Feb  2008 - added thkdrg.
c> Feb  2008 - added sshflg and steric,sshgmn,thmean.
c> Jun  2008 - added tilefq.
c> Mar  2008 - added dswflg and diurnl.
c> Dec  2008 - difsmo is now an integer number of layers.
c> Jan  2009 - added arcend.
c-----------------------------------------------------------------------------
