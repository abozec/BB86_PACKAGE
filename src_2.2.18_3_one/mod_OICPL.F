#define OI_REGRID_CHOICE 0  /* ocean and seaice have the same grid and layout */

      module mod_OICPL
c
c --- ESMF Framework module
      use ESMF_Mod
c
      implicit none
      private
c
      public OICPL_SetServices
c
c --- phase info
      integer, parameter, public :: ice2ocn_phase = 1
      integer, parameter, public :: ocn2ice_phase = 2
c
c --- VM and PET info
      type(ESMF_VM), save :: vm
      integer,       save :: petCount,localPet
c --- Misc field stuff
      type(ESMF_Bundle),         save :: ocnBundle,    iceBundle
      integer(ESMF_KIND_I4),     save :: numOcnFields, numIceFields
      type(ESMF_Field),          save :: ocnField,     iceField
      real(ESMF_KIND_R4),pointer,save :: ocnData(:,:), iceData(:,:)
#if OI_REGRID_CHOICE != 0
c
c --- Route handles for regridding and redistribution
      type(ESMF_RouteHandle), allocatable, save ::
     &    i2oRouteHandle(:), o2iRouteHandle(:)
#endif /* OI_REGRID_CHOICE != 0 */

      contains

      subroutine OICPL_SetServices(cplComp, rc)
c
      type(ESMF_CplComp) :: cplComp
      integer :: rc
c
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_SETINIT,
     &     OICPL_Init,
     &     ESMF_SINGLEPHASE,
     &     rc)
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_SETRUN,
     &     OICPL_Run_I2O,
     &     ice2ocn_phase,
     &     rc)
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_SETRUN,
     &     OICPL_Run_O2I,
     &     ocn2ice_phase,
     &     rc)
      call ESMF_CplCompSetEntryPoint(
     &     cplComp,
     &     ESMF_SETFINAL,
     &     OICPL_Final,
     &     ESMF_SINGLEPHASE,
     &     rc)
c
      end subroutine OICPL_SetServices


      subroutine OICPL_Init(cplComp, impState, expState, extClock, rc)
c
c --- Calling parameters
      type(ESMF_CplComp) :: cplComp
      type(ESMF_State)   :: impState
      type(ESMF_State)   :: expState
      type(ESMF_Clock)   :: extClock
      integer :: rc
c
c --- Locals
      integer :: i
      character(ESMF_MAXSTR) :: msg
      type(ESMF_State) :: oiState, oeState, iiState, ieState
c
c --- Report
      call ESMF_LogWrite("OICPL initialize routine called",
     &     ESMF_LOG_INFO)
!-----call ESMF_LogFlush
c
c --- Get VM
      call ESMF_CplCompGet(cplComp, vm=vm, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "Get VM failed",
     &   rc)) call ESMF_Finalize(rc=rc)
c
c --- Get PET info
      call ESMF_VMGet(vm, petCount=petCount, localPET=localPet, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "Get VM info failed",
     &   rc)) call ESMF_Finalize(rc=rc)
c
c --- Get OCEAN and SEAICE import states
      CALL ESMF_StateGetState(impState, "OCEAN Import", oiState, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "Get OCEAN impState failed",
     &   rc)) call ESMF_Finalize(rc=rc)
      CALL ESMF_StateGetState(impState, "SEAICE Import", iiState, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "Get SEAICE impState failed",
     &   rc)) call ESMF_Finalize(rc=rc)
c
c --- Get OCEAN and SEAICE export states
      CALL ESMF_StateGetState(expState, "OCEAN Export", oeState, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "Get OCEAN expState failed",
     &   rc)) call ESMF_Finalize(rc=rc)
      CALL ESMF_StateGetState(expState, "SEAICE Export", ieState, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "Get SEAICE expState failed",
     &   rc)) call ESMF_Finalize(rc=rc)
c
c --- Initialize I2O
c
c --- Get bundle for ice
      call ESMF_StateGetBundle(ieState, "CICE",  iceBundle, rc=rc)
c
c --- Get number of ice fields from bundle
      call ESMF_BundleGet(iceBundle, fieldCount=numIceFields, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "BundleGet numIceFields failed",
     &   rc)) call ESMF_Finalize(rc=rc)
#if OI_REGRID_CHOICE != 0
c
c --- Allocate route handles for I2O
      allocate(i2oRouteHandle(numIceFields))
c
c --- Setup I2O route handles
      call ESMF_StateGetBundle(oiState, "HYCOM", ocnBundle, rc=rc)
      do i = 1,numIceFields
c       this works because ocn and ice fields are in the same order.
        call ESMF_BundleGetField(iceBundle, i, iceField, rc=rc)
        call ESMF_BundleGetField(ocnBundle, i, ocnField, rc=rc)
#if OI_REGRID_CHOICE == 1
        i2oRouteHandle(i) = ESMF_RouteHandleCreate(rc)
        call ESMF_FieldRedistStore(iceField, ocnField, vm,
     &       routeHandle=i2oRouteHandle(i), rc=rc)
        call ESMF_FieldRedist(iceField, ocnField,
     &       i2oRouteHandle(i), rc=rc)
#elif OI_REGRID_CHOICE == 2
        i2oRouteHandle(i) = ESMF_RouteHandleCreate(rc)
        call ESMF_FieldRegridStore(iceField, ocnField, vm,
     &       routeHandle=i2oRouteHandle(i),
     &       regridMethod=ESMF_REGRID_METHOD_BILINEAR, rc=rc)
        call ESMF_FieldRegrid(iceField, ocnField,
     &       i2oRouteHandle(i), rc=rc)
#endif
      enddo !i
#endif /* OI_REGRID_CHOICE != 0 */
c
c  Initialize O2I
c
c
c --- Get bundle for ocn
      call ESMF_StateGetBundle(oeState, "HYCOM", ocnBundle, rc=rc)
c
c --- Get number of ocn fields from bundle
      call ESMF_BundleGet(ocnBundle, fieldCount=numOcnFields, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "BundleGet numOcnFields failed",
     &   rc)) call ESMF_Finalize(rc=rc)
#if OI_REGRID_CHOICE != 0
c
c --- Allocate route handles for O2I
      allocate(o2iRouteHandle(numOcnFields))
c
c --- Setup O2I route handles
      call ESMF_StateGetBundle(iiState, "CICE",  iceBundle, rc=rc)
      do i = 1,numOcnFields
c       this works because ocn and ice fields are in the same order.
        call ESMF_BundleGetField(ocnBundle, i, ocnField, rc=rc)
        call ESMF_BundleGetField(iceBundle, i, iceField, rc=rc)
#if OI_REGRID_CHOICE == 1
        o2iRouteHandle(i) = ESMF_RouteHandleCreate(rc)
        call ESMF_FieldRedistStore(ocnField, iceField, vm,
     &       routeHandle=o2iRouteHandle(i),
     &       rc=rc)
        call ESMF_FieldRedist(ocnField, iceField,
     &       o2iRouteHandle(i), rc=rc)
#elif OI_REGRID_CHOICE == 2
        o2iRouteHandle(i) = ESMF_RouteHandleCreate(rc)
        call ESMF_FieldRegridStore(ocnField, iceField, vm,
     &       routeHandle=o2iRouteHandle(i),
     &       regridMethod=ESMF_REGRID_METHOD_BILINEAR,
     &       rc=rc)
        call ESMF_FieldRegrid(ocnField, iceField,
     &       o2iRouteHandle(i), rc=rc)
#endif
      enddo !i
#endif /* OI_REGRID_CHOICE != 0 */
c
      return
      end subroutine OICPL_Init

      subroutine OICPL_Run_I2O(cplComp, impState, expState, extClock,
     &                         rc)
c
c --- Calling parameters
      type(ESMF_CplComp) :: cplComp
      type(ESMF_State)   :: impState
      type(ESMF_State)   :: expState
      type(ESMF_Clock)   :: extClock
      integer :: rc
c
c --- Locals
      integer :: i
      character(ESMF_MAXSTR) :: msg
      type(ESMF_State) :: oiState, oeState, iiState, ieState
c
c --- Report
      write(msg,'(a,i3)')
     &   "OICPL I2O run routine called, numIceFields =",numIceFields
      call ESMF_LogWrite(trim(msg), ESMF_LOG_INFO) !OICPL I2O run
!-----call ESMF_LogFlush
c
c --- Get OCEAN import state
      CALL ESMF_StateGetState(impState, "OCEAN Import", oiState, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "Get OCEAN impState failed",
     &   rc)) call ESMF_Finalize(rc=rc)
c
c --- Get SEAICE export state
      CALL ESMF_StateGetState(expState, "SEAICE Export", ieState, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "Get SEAICE expState failed",
     &   rc)) call ESMF_Finalize(rc=rc)
c
c --- Get bundle for ice
      call ESMF_StateGetBundle(ieState, "CICE",  iceBundle, rc=rc)
c
c --- Transfer fields from ice state to ocn state
      call ESMF_StateGetBundle(oiState, "HYCOM", ocnBundle, rc=rc)
      do i = 1,numIceFields
c       this works because ocn and ice fields are in the same order.
        call ESMF_BundleGetField(iceBundle, i, iceField, rc=rc)
        IF (ESMF_LogMsgFoundError(rc, "BundleGetField iceEx failed",
     &   rc)) call ESMF_Finalize(rc=rc)
        call ESMF_BundleGetField(ocnBundle, i, ocnField, rc=rc)
        IF (ESMF_LogMsgFoundError(rc, "BundleGetField ocnIm failed",
     &   rc)) call ESMF_Finalize(rc=rc)
#if OI_REGRID_CHOICE == 0
        call ESMF_FieldGetDataPointer(iceField, iceData,
     &       copyFlag=ESMF_DATA_REF, rc=rc)
        IF (ESMF_LogMsgFoundError(rc, "FieldGetDataP iceEx failed",
     &   rc)) call ESMF_Finalize(rc=rc)
        call ESMF_FieldGetDataPointer(ocnField, ocnData,
     &       copyFlag=ESMF_DATA_REF, rc=rc)
        IF (ESMF_LogMsgFoundError(rc, "FieldGetDataP ocnIm failed",
     &   rc)) call ESMF_Finalize(rc=rc)
        ocnData(:,:) = iceData(:,:)
#elif OI_REGRID_CHOICE == 1
        call ESMF_FieldRedist(iceField, ocnField,
     &       i2oRouteHandle(i), rc=rc)
#elif OI_REGRID_CHOICE == 2
        call ESMF_FieldRegrid(iceField, ocnField,
     &       i2oRouteHandle(i), rc=rc)
#endif
      enddo !i
c
      return
      end subroutine OICPL_Run_I2O

      subroutine OICPL_Run_O2I(cplComp, impState, expState, extClock,
     &                         rc)
c
c --- Calling parameters
      type(ESMF_CplComp) :: cplComp
      type(ESMF_State) :: impState
      type(ESMF_State) :: expState
      type(ESMF_Clock) :: extClock
      integer :: rc
c
c --- Locals
      integer :: i
      character(ESMF_MAXSTR) :: msg
      type(ESMF_State) :: oiState, oeState, iiState, ieState
c
c --- Report
      write(msg,'(a,i3)')
     &   "OICPL O2I run routine called, numOcnFields =",numOcnFields
      call ESMF_LogWrite(trim(msg), ESMF_LOG_INFO) !OICPL O2I run
!-----call ESMF_LogFlush
c
c --- Get SEAICE import state
      CALL ESMF_StateGetState(impState, "SEAICE Import", iiState, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "Get SEAICE impState failed",
     &   rc)) call ESMF_Finalize(rc=rc)
c
c --- Get OCEAN export state
      CALL ESMF_StateGetState(expState, "OCEAN Export", oeState, rc=rc)
      IF (ESMF_LogMsgFoundError(rc, "Get OCEAN expState failed",
     &   rc)) call ESMF_Finalize(rc=rc)
c
c --- Get bundle for first nest of ocn
      call ESMF_StateGetBundle(oeState, "HYCOM", ocnBundle, rc=rc)
c
c --- Transfer fields from ocn state to ice state
      call ESMF_StateGetBundle(iiState, "CICE",  iceBundle, rc=rc)
      do i = 1,numOcnFields
c       this works because ocn and ice fields are in the same order.
        call ESMF_BundleGetField(ocnBundle, i, ocnField, rc=rc)
        IF (ESMF_LogMsgFoundError(rc, "BundleGetField ocnEx failed",
     &   rc)) call ESMF_Finalize(rc=rc)
        call ESMF_BundleGetField(iceBundle, i, iceField, rc=rc)
        IF (ESMF_LogMsgFoundError(rc, "BundleGetField iceIm failed",
     &   rc)) call ESMF_Finalize(rc=rc)
#if OI_REGRID_CHOICE == 0
        call ESMF_FieldGetDataPointer(ocnField, ocnData,
     &       copyFlag=ESMF_DATA_REF, rc=rc)
        IF (ESMF_LogMsgFoundError(rc, "FieldGetDataP ocnEx failed",
     &   rc)) call ESMF_Finalize(rc=rc)
        call ESMF_FieldGetDataPointer(iceField, iceData,
     &       copyFlag=ESMF_DATA_REF, rc=rc)
        IF (ESMF_LogMsgFoundError(rc, "FieldGetDataP iceIm failed",
     &   rc)) call ESMF_Finalize(rc=rc)
        iceData(:,:) = ocnData(:,:)
#elif OI_REGRID_CHOICE == 1
        call ESMF_FieldRedist(ocnField, iceField,
     &       o2iRouteHandle(i), rc=rc)
#elif OI_REGRID_CHOICE == 2
        call ESMF_FieldRegrid(ocnField, iceField,
     &       o2iRouteHandle(i), rc=rc)
#endif
      enddo !i
c
      return
      end subroutine OICPL_Run_O2I

      subroutine OICPL_Final(cplComp, impState, expState, extClock, rc)
c
c --- Calling parameters
      type(ESMF_CplComp) :: cplComp
      type(ESMF_State) :: impState
      type(ESMF_State) :: expState
      type(ESMF_Clock) :: extClock
      integer :: rc
c
c --- Locals
      integer :: i
      character(ESMF_MAXSTR) :: msg
      type(ESMF_State) :: oiState, oeState, iiState, ieState
c
c --- Report
      call ESMF_LogWrite("OICPL finalize routine called", ESMF_LOG_INFO)
!-----call ESMF_LogFlush
#if OI_REGRID_CHOICE != 0
c
c --- Release i2o regrid/redist route handles
      do i = 1,numIceFields
#if OI_REGRID_CHOICE == 1
        call ESMF_FieldRedistRelease(i2oRouteHandle(i), rc=rc)
#elif OI_REGRID_CHOICE == 2
        call ESMF_FieldRegridRelease(i2oRouteHandle(i), rc=rc)
#endif
      enddo
c
c --- Release o2i regrid/redist route handles
      do i = 1,numOcnFields
#if OI_REGRID_CHOICE == 1
        call ESMF_FieldRedistRelease(o2iRouteHandle(i), rc=rc)
#elif OI_REGRID_CHOICE == 2
        call ESMF_FieldRegridRelease(o2iRouteHandle(i), rc=rc)
#endif
      enddo
c
c --- Deallocate regrid/redist route handles
      deallocate(i2oRouteHandle,o2iRouteHandle)
#endif /* OI_REGRID_CHOICE != 0 */
c
      return
      end subroutine OICPL_Final

      end module mod_OICPL
