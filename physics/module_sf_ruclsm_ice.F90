#define LSMRUC_DBG_LVL 3000
!WRF:MODEL_LAYER:PHYSICS
!
MODULE module_sf_ruclsm_ice

! Notes for perturbations of soil properties (Judith Berner)
! Perturbations are applied in subroutine soilprob to array hydro;
! soilprop is called from subroutine SICETMP which is called from subroutine LSMRUC;
! subroutine LSMRUC had two new 3D fields: pattern_spp_lsm (in) and field_sf(inout);
!    their vertical dimension is number of atmospheric levels (kms:kme) - (suboptimal, but easiest hack)
!    field_sf is used to pass perturbed fields of hydrop up to model (and output) driver;
! in argument list to SICETMP the arrays are passed as pattern_spp_lsm(i,1:nzs,j), and exist henceforth as
! column arrays;
! in the subroutines below SICETMP (SNOW and SNOWSOIL) the fields are called rstochcol,fieldcol_sf
! to reflect their dimension rstochcol (1:nzs)


  USE module_model_constants
  USE module_wrf_error

! VEGETATION PARAMETERS
        INTEGER :: LUCATS , BARE, NATURAL, CROP, URBAN
        integer, PARAMETER :: NLUS=50
        CHARACTER*8 LUTYPE
        INTEGER, DIMENSION(1:NLUS) :: IFORTBL
        real, dimension(1:NLUS) ::  SNUPTBL, RSTBL, RGLTBL, HSTBL, LAITBL,         &
                                    ALBTBL, Z0TBL, LEMITBL, PCTBL, SHDTBL, MAXALB
        REAL ::   TOPT_DATA,CMCMAX_DATA,CFACTR_DATA,RSMAX_DATA
! SOIL PARAMETERS
        INTEGER :: SLCATS
        INTEGER, PARAMETER :: NSLTYPE=30
        CHARACTER*8 SLTYPE
        REAL, DIMENSION (1:NSLTYPE) :: BB,DRYSMC,HC,                           &
        MAXSMC, REFSMC,SATPSI,SATDK,SATDW, WLTSMC,QTZ

! LSM GENERAL PARAMETERS
        INTEGER :: SLPCATS
        INTEGER, PARAMETER :: NSLOPE=30
        REAL, DIMENSION (1:NSLOPE) :: SLOPE_DATA
        REAL ::  SBETA_DATA,FXEXP_DATA,CSOIL_DATA,SALP_DATA,REFDK_DATA,           &
                 REFKDT_DATA,FRZK_DATA,ZBOT_DATA,  SMLOW_DATA,SMHIGH_DATA,        &
                        CZIL_DATA

        CHARACTER*256  :: err_message


CONTAINS

!-----------------------------------------------------------------
    SUBROUTINE LSMRUC_ice(spp_lsm,                               &
#if (EM_CORE==1)
                   pattern_spp_lsm,field_sf,                     &
#endif
                   DT,KTAU,NSL,                                  &
#if (EM_CORE==1)
                   lakemodel,lakemask,                           &
                   graupelncv,snowncv,rainncv,                   &
#endif
                   ZS,RAINBL,SNOW,SNOWH,SNOWC,FRZFRAC,frpcpn,    &
                   rhosnf,precipfr,                              & ! pass it out to module_diagnostics
                   Z3D,P8W,T3D,QV3D,QC3D,RHO3D,                  & !p8W in [PA]
                   GLW,GSW,EMISS,CHKLOWQ, CHS,                   & 
                   FLQC,FLHC,ALB,ZNT,                            &
                   Z0,SNOALB,ALBBCK,mminlu,                      &  !new
                   QSFC,QSG,QVG,QCG,DEW,SOILT1,TSNAV,            &
                   TBOT,IVGTYP,ISLTYP,XLAND,                     &
                   ISWATER,ISICE,XICE,XICE_THRESHOLD,            &
                   CP,ROVCP,G0,LV,STBOLT,                        &
                   SOILMOIS,SH2O,TSO,SOILT,HFX,QFX,LH,           &
                   SFCRUNOFF,ACRUNOFF,                           &
                   SFCEVP,GRDFLX,SNOWFALLAC,ACSNOW,SNOM,         &
                   SMFR3D,KEEPFR3DFLAG,myjpbl,                   &
                   ids,ide, jds,jde, kds,kde,                    &
                   ims,ime, jms,jme, kms,kme,                    &
                   its,ite, jts,jte, kts,kte                     &
                         )                      

!-----------------------------------------------------------------
   IMPLICIT NONE
!-----------------------------------------------------------------
!
! The RUC LSM model is described in:
!  Smirnova, T.G., J.M. Brown, and S.G. Benjamin, 1997: 
!     Performance of different soil model configurations in simulating 
!     ground surface temperature and surface fluxes. 
!     Mon. Wea. Rev. 125, 1870-1884.
!  Smirnova, T.G., J.M. Brown, and D. Kim, 2000: Parameterization of 
!     cold-season processes in the MAPS land-surface scheme. 
!     J. Geophys. Res. 105, 4077-4086.
!-----------------------------------------------------------------
!-- DT            time step (second)
!        ktau - number of time step
!        NSL  - number of seaice layers
!        NZS  - number of levels in seaice
!        ZS   - depth of seaice levels (m)
!-- RAINBL    - accumulated rain in [mm] between the PBL calls
!-- RAINNCV         one time step grid scale precipitation (mm/step)
!        SNOW - snow water equivalent [mm]
!        FRAZFRAC - fraction of frozen precipitation
!-- PRECIPFR (mm) - time step frozen precipitation
!-- SNOWC       flag indicating snow coverage (1 for snow cover)
!-- Z3D         heights (m)
!-- P8W         3D pressure (Pa)
!-- T3D         temperature (K)
!-- QV3D        3D water vapor mixing ratio (Kg/Kg)
!        QC3D - 3D cloud water mixing ratio (Kg/Kg)
!       RHO3D - 3D air density (kg/m^3)
!-- GLW         downward long wave flux at ground surface (W/m^2)
!-- GSW         absorbed short wave flux at ground surface (W/m^2)
!-- EMISS       surface emissivity (between 0 and 1)
!        FLQC - surface exchange coefficient for moisture (kg/m^2/s)
!        FLHC - surface exchange coefficient for heat [W/m^2/s/degreeK]     
!         ALB - surface albedo (between 0 and 1)
!      SNOALB - maximum snow albedo (between 0 and 1)
!      ALBBCK - snow-free albedo (between 0 and 1)
!         ZNT - roughness length [m]
!-- TBOT        seaice temperature at lower boundary (K)
!-- XLAND       land mask (1 for land, 2 for water)
!-- CP          heat capacity at constant pressure for dry air (J/kg/K)
!-- G0          acceleration due to gravity (m/s^2)
!-- LV          latent heat of melting (J/kg)
!-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
!    SOILMOIS - seaice moisture content (volumetric fraction=1.)
!         TSO - seaice temp (K)
!-- SOILT       surface temperature (K)
!-- HFX         upward heat flux at the surface (W/m^2)
!-- QFX         upward moisture flux at the surface (kg/m^2/s)
!-- LH          upward latent heat flux (W/m^2)
!   SFCRUNOFF - ground surface runoff [mm]
!   ACRUNOFF - run-total surface runoff [mm]
!   SFCEVP - total evaporation in [kg/m^2]
!   GRDFLX - seaice heat flux (W/m^2: negative, if downward from surface)
!   SNOWFALLAC - run-total snowfall accumulation [m]   
!   ACSNOW - run-toral SWE of snowfall [mm]   
!-- CHKLOWQ - is either 0 or 1 (so far set equal to 1).
!--           used only in MYJPBL. 
!-- tice - sea ice temperture (C)
!-- rhosice - sea ice density (kg m^-3)
!-- capice - sea ice volumetric heat capacity (J/m^3/K)
!-- thdifice - sea ice thermal diffusivity (m^2/s)
!--
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-------------------------------------------------------------------------

   INTEGER,     PARAMETER            ::     nvegclas=24+3

   REAL,       INTENT(IN   )    ::     DT
   LOGICAL,    INTENT(IN   )    ::     myjpbl,frpcpn
   INTEGER,    INTENT(IN   )    ::     spp_lsm
   INTEGER,    INTENT(IN   )    ::     ktau, nsl, isice, iswater, &
                                       ims,ime, jms,jme, kms,kme, &
                                       ids,ide, jds,jde, kds,kde, &
                                       its,ite, jts,jte, kts,kte

#if (EM_CORE==1)
   REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ),OPTIONAL::    pattern_spp_lsm
   REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ),OPTIONAL::    field_sf
#endif
   REAL,    DIMENSION( ims:ime, 1  :nsl, jms:jme )         ::    field_sf_loc

   REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )            , &
            INTENT(IN   )    ::                           QV3D, &
                                                          QC3D, &
                                                           p8w, &
                                                         rho3D, &
                                                           T3D, &
                                                           z3D

   REAL,       DIMENSION( ims:ime , jms:jme ),                   &
               INTENT(IN   )    ::                       RAINBL, &
                                                            GLW, &
                                                            GSW, &
                                                         ALBBCK, &
                                                           CHS , &
                                                           XICE, &
                                                          XLAND, &
                                                           TBOT

   REAL,       DIMENSION( ims:ime , jms:jme ),                   &
               INTENT(INOUT   )    ::                 FLQC,FLHC


#if (EM_CORE==1)
   REAL,       OPTIONAL, DIMENSION( ims:ime , jms:jme ),         &
               INTENT(IN   )    ::                   GRAUPELNCV, &
                                                        SNOWNCV, &
                                                        RAINNCV
   REAL,       DIMENSION( ims:ime , jms:jme ),                   &
               INTENT(IN   )    ::                     lakemask
   INTEGER,    INTENT(IN   )    ::                    LakeModel
#endif

   REAL,       DIMENSION( 1:nsl), INTENT(IN   )      ::      ZS

   REAL,       DIMENSION( ims:ime , jms:jme ),                   &
               INTENT(INOUT)    ::                               &
                                                           SNOW, &
                                                          SNOWH, &
                                                          SNOWC, &
                                                         SNOALB, &
                                                            ALB, &
                                                          EMISS, &
                                                            Z0 , &
                                                            ZNT

   REAL,       DIMENSION( ims:ime , jms:jme ),                   &
               INTENT(IN   )    ::                               &
                                                        FRZFRAC

   INTEGER,    DIMENSION( ims:ime , jms:jme ),                   &
               INTENT(IN   )    ::                       IVGTYP, &
                                                         ISLTYP
   CHARACTER(LEN=*), INTENT(IN   )    ::                 MMINLU

   REAL, INTENT(IN   )          ::         CP,ROVCP,G0,LV,STBOLT,XICE_threshold
 
   REAL,       DIMENSION( ims:ime , 1:nsl, jms:jme )           , &
               INTENT(INOUT)    ::                 SOILMOIS,SH2O,TSO

   REAL,       DIMENSION( ims:ime, jms:jme )                   , &
               INTENT(INOUT)    ::                        SOILT, &
                                                            HFX, &
                                                            QFX, &
                                                             LH, &
                                                         SFCEVP, &
                                                      SFCRUNOFF, &
                                                       ACRUNOFF, &
                                                         GRDFLX, &
                                                         ACSNOW, &
                                                           SNOM, &
                                                            QVG, &
                                                            QCG, &
                                                            DEW, &
                                                           QSFC, &
                                                            QSG, &
                                                        CHKLOWQ, &
                                                         SOILT1, &
                                                          TSNAV

   REAL,       DIMENSION( its:ite, jts:jte )    ::               &
                                                             PC, &
                                                        RUNOFF1, &
                                                        RUNOFF2, &
                                                         EMISSL, &
                                                           ZNTL, &
                                                          SMELT, &
                                                           SNOH, &
                                                          SNFLX, &
                                                         SUBLIM, &
                                                           sflx, &
                                                          EVAPL, &
                                                          PRCPL, &
                                                         SEAICE
! Energy and water budget variables:
   REAL,       DIMENSION( its:ite, jts:jte )    ::               &
                                                         budget, &
                                                       acbudget, &
                                                    waterbudget, &
                                                        snowold, & 
                                                  acwaterbudget


   REAL,       DIMENSION( ims:ime, 1:nsl, jms:jme)               &
                                             ::    KEEPFR3DFLAG, &
                                                         SMFR3D


   REAL,       DIMENSION( ims:ime, jms:jme ), INTENT(OUT)     :: &
                                                         RHOSNF, & !RHO of snowfall
                                                       PRECIPFR, & ! time-step frozen precip
                                                     SNOWFALLAC
!--- seaice/snow properties
   REAL                                                          &
                             ::                                  &
                                                       RHONEWSN, &
                                                          RHOSN, &
                                                      RHOSNFALL, &
                                                       SNOWFRAC, &
                                                          SNHEI, &
                                                           SNWE

   REAL                                      ::              CW, &
                                                           C1SN, &
                                                           C2SN


   REAL,     DIMENSION(1:NSL)                ::          ZSMAIN, &
                                                         ZSHALF, &
                                                         DZhalf, &
                                                         DZmain, &
                                                         DTDZS2

   REAL,     DIMENSION(1:2*(nsl-2))          ::           DTDZS

   REAL,     DIMENSION(1:5001)               ::             TBQ


   REAL,     DIMENSION( 1:nsl )              ::                  & 
                                                          TSO1D

   REAL                           ::                        RSM, &
                                                      SNWEPRINT, &
                                                     SNHEIPRINT

   REAL                           ::                     PRCPMS, &
                                                        NEWSNMS, &
                                                      prcpncliq, &
                                                       prcpncfr, &
                                                      prcpculiq, &
                                                       prcpcufr, &
                                                           PATM, &
                                                          PATMB, &
                                                           TABS, &
                                                          QVATM, &
                                                          QCATM, &
                                                          Q2SAT, &
                                                         CONFLX, &
                                                            RHO, &
                                                           QKMS, &
                                                           TKMS, &
                                                        snowrat, &
                                                       grauprat, &
                                                       graupamt, &
                                                         icerat, &
                                                          curat
   REAL      ::  cq,r61,r273,arp,brp,x,evs,eis
   real      ::  zsmain_dop, zshalf_dop, tso_dop

   REAL      ::  ac,as, wb
 
   INTEGER   ::  I,J,K,NZS,NZS1,NDDZS
   INTEGER   ::  k1,l,k2,kp,km
   CHARACTER (LEN=132) :: message

   REAL,DIMENSION(ims:ime,1:nsl,jms:jme) :: rstoch 

!-----------------------------------------------------------------
         ! thicknesses centered at soil levels
         DZhalf = (/ 0.005 , 0.02, 0.045 , 0.13 , 0.25 , 0.35 , 0.50 , 1.00 , 3.15 /) ! ruclsm
         !DZhalf = (/ 0.005 , 0.01, 0.05 , 0.07 , 0.33 , 0.27 , 0.53 , 0.67 , 2.23 /) ! prlsm

         NZS=NSL
         NDDZS=2*(nzs-2)
         !-----
         NZS1=NZS-1
         !--- Levels in soil
         zsmain(1)=0.
         zshalf(1)=0.
         zshalf(2)=dzhalf(1)
         do k=2,nzs
           zsmain(k)=zs(k)
           dzmain(k-1)=zs(k)-zs(k-1)
         enddo
         dzmain(nzs)=4.29

         do k=3,nzs
           zshalf(k)= zshalf(k-1) + dzhalf(k-1) !0.5*(zsmain(k-1) + zsmain(k))
         enddo
         zshalf_dop = zshalf(nzs) + dzhalf(nzs)
         zsmain_dop = zsmain(nzs) + dzmain(nzs)
     !-------------------------------------------------------------
     !-- DDZS and DSDZ1 are for implicit solution of termal diffusion
     !-------------------------------------------------------------
     IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,' DT,NZS1, ZSMAIN, ZSHALF --->', dt,nzs1,zsmain,zshalf
     ENDIF
         DTDZS2 = 0.
         DTDZS  = 0.
       DO  K=2,NZS1
         K1=2*K-3
         K2=K1+1
         X=DT/2./(ZSHALF(K+1)-ZSHALF(K))
         DTDZS(K1)=X/(ZSMAIN(K)-ZSMAIN(K-1))
         DTDZS2(K-1)=X
         DTDZS(K2)=X/(ZSMAIN(K+1)-ZSMAIN(K))
       END DO
         X=DT/2./(ZSHALF_dop-ZSHALF(nzs))
         DTDZS(2*nzs-3)=X/(ZSMAIN(nzs)-ZSMAIN(nzs1))
         DTDZS2(nzs1)=X
         DTDZS(2*nzs-2)=X/(ZSMAIN_dop-ZSMAIN(nzs))

         CW =4.183E6

!***********************************************************************
!--- Constants for snow density calculations C1SN and C2SN

         c1sn=0.026
         !c1sn=0.01
         c2sn=21.
         rhonewsn = 100.

!***********************************************************************


         rstoch(ims:ime,1:nsl,jms:jme)=0.0
         field_sf_loc(ims:ime,1:nsl,jms:jme)=0.0
!beka added
#if (EM_CORE==1)
       if (spp_lsm==1) then
         do J=jts,jte
           do i=its,ite
             do k=1,nsl
               rstoch(i,k,j) = pattern_spp_lsm(i,k,j)
               field_sf_loc(i,k,j)=field_sf(i,k,j)
             enddo
           enddo
         enddo 
  !if(ktau.eq.10 .and. ktau.eq.20 .and. ktau.eq.30) then
  ! print *,'ktau, min/max rstoch',ktau,minval(rstoch(:,1,:)),maxval(rstoch(:,1,:))
  !endif
       endif  
#endif
!---- table TBQ is for resolution of balance equation in VILKA
        CQ=173.15-.05
        R273=1./273.15
        R61=6.1153*0.62198
        ARP=77455.*41.9/461.525
        BRP=64.*41.9/461.525

        DO K=1,5001
          CQ=CQ+.05
!          TBQ(K)=R61*EXP(ARP*(R273-1./CQ)-BRP*LOG(CQ*R273))
        EVS=EXP(17.67*(CQ-273.15)/(CQ-29.65))
        EIS=EXP(22.514-6.15E3/CQ)
        if(CQ.ge.273.15) then
          ! tbq is in mb
          tbq(k) = R61*evs
        else
          tbq(k) = R61*eis
        endif

        END DO

   SEAICE = 0.

#if ( NMM_CORE == 1 )
     if(ktau+1.eq.1) then
#else
     if(ktau.eq.1) then
#endif
     DO J=jts,jte
         DO i=its,ite
#if (EM_CORE==1)
     if(lakemodel==1. .and. lakemask(i,j)==1.) goto 299
!Lakes
#endif
      IF((XLAND(I,J)-1.5).GT.0.)THEN
      !-- Water 
        SEAICE(i,j)=0.0
      ELSE
      ! LAND OR SEA ICE
        if ( ( XICE(I,J) .GE. XICE_THRESHOLD ) .and. ( XICE(i,j) .LE. 1.0 ) ) THEN
        !-- sea ice
           SEAICE(i,j)=1.
        else
           SEAICE(i,j)=0.0
        endif
      ENDIF

           IF(SEAICE(I,J).GT.0.5)THEN
           !-- Sea-ice case
           !--- Initialize seaice parameters at initial time step
           !--- initializing snow fraction, thereshold = 32 mm of snow water 
           !    or ~100 mm of snow height
           !
             if(snow(i,j) > 0. .and. snowc(i,j) <= 0.) then
               snowc(i,j) = min(1.,snow(i,j)/32.)
             endif
             !--- initializing inside snow temp if it is not defined
             IF((soilt1(i,j) .LT. 1200.) .or. (soilt1(i,j) .GT.350.)) THEN
               IF(snow(i,j).gt.32.) THEN
                 soilt1(i,j)=0.5*(soilt(i,j)+tso(i,1,j))
                IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
                  WRITE ( message , FMT='(A,F8.3,2I6)' ) &
                    'Temperature inside snow is initialized in RUCLSM ', soilt1(i,j),i,j
                  CALL wrf_debug ( 0 , message )
                ENDIF
               ELSE
                 soilt1(i,j) = tso(i,1,j)
               ENDIF
             ENDIF
             tsnav(i,j) =0.5*(soilt(i,j)+tso(i,1,j))-273.15
             patmb=P8w(i,kms,j)*1.e-2
             QSG  (i,j) = QSN(SOILT(i,j),TBQ)/PATMB
             IF((qvg(i,j) .LE. 0.) .or. (qvg(i,j) .GT.0.1)) THEN
               qvg  (i,j) = qv3d(i,1,j)
               qcg  (i,j) = 0.
               IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
                WRITE ( message , FMT='(A,3F8.3,2I6)' ) &
                  'QVG is initialized in RUCLSM ', qvg(i,j),qsg(i,j),i,j
                CALL wrf_debug ( 0 , message )
               ENDIF
             ENDIF
             qsfc(i,j) = qvg(i,j)/(1.+qvg(i,j))
             SMELT(i,j) = 0.
             SNOM (i,j) = 0.
             !ACSNOW(i,j) = 0.
             SNOWFALLAC(i,j) = 0.
             PRECIPFR(i,j) = 0.
             RHOSNF(i,j) = -1.e3 ! non-zero flag
             SNFLX(i,j) = 0.
             DEW  (i,j) = 0.
             zntl (i,j) = 0.
             RUNOFF1(i,j) = 0.
             RUNOFF2(i,j) = 0.
             SFCRUNOFF(i,j) = 0.
             ACRUNOFF(i,j) = 0.
             emissl (i,j) = 0.
             budget(i,j) = 0.
             acbudget(i,j) = 0.
             waterbudget(i,j) = 0.
             acwaterbudget(i,j) = 0.
            
             ! For RUC LSM CHKLOWQ needed for MYJPBL should 
             ! 1 because is actual specific humidity at the surface, and
             ! not the saturation value
             chklowq(i,j) = 1.
             snoh  (i,j) = 0.
             sublim(i,j) = 0.
             sflx  (i,j) = 0.
             evapl (i,j) = 0.
             prcpl (i,j) = 0.

           ENDIF ! sea ice

  299      continue ! lakes

         ENDDO  ! end i loop
       ENDDO  ! end j loop
     endif   ! end initizlize
        
!***********************************************************************
   ! Main integration loop for sea ice
   DO J=jts,jte
     DO i=its,ite

#if (EM_CORE==1)
     if(lakemodel==1. .and. lakemask(i,j)==1.) goto 2999 !Lakes
#endif
      IF((XLAND(I,J)-1.5).GT.0.)THEN
      !-- Water 
        SEAICE(i,j)=0.0
        SNOW(I,J)=0.0
        SNOWH(I,J)=0.0
        SNOWC(I,J)=0.0
        ! accumulated water equivalent of frozen precipitation over water [mm]
        ! acsnow(i,j)=acsnow(i,j)+precipfr(i,j)

        patmb=P8w(i,1,j)*1.e-2
        qvg  (i,j) = QSN(SOILT(i,j),TBQ)/PATMB
        qsfc(i,j) = qvg(i,j)/(1.+qvg(i,j))
        CHKLOWQ(I,J)=1.
        Q2SAT=QSN(TABS,TBQ)/PATMB

        DO K=1,NZS
          SOILMOIS(I,K,J)=1.0
          SH2O    (I,K,J)=1.0 
          TSO(I,K,J)= SOILT(I,J)
        ENDDO

      IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
          PRINT*,'  water point, I=',I,'J=',J,'SOILT=', SOILT(i,j)
      ENDIF

      ELSE 
      ! LAND OR SEA ICE
        if ( ( XICE(I,J) .GE. XICE_THRESHOLD ) .and. ( XICE(i,j) .LE. 1.0 ) ) then
        !-- sea ice
           SEAICE(i,j)=1.

      IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
          PRINT*,'  ice point, I=',I,'J=',J,xice(i,j),'SOILT=', SOILT(i,j)
      ENDIF
 
           ZNT(I,J)    = 0.011
           snoalb(i,j) = 0.75
           emissl(i,j) = 0.98

           patmb=P8w(i,1,j)*1.e-2
           qvg  (i,j) = QSN(SOILT(i,j),TBQ)/PATMB
           qsg  (i,j) = qvg(i,j)
           qsfc (i,j) = qvg(i,j)/(1.+qvg(i,j))

           DO K=1,NZS
             soilmois(i,k,j) = 1.
             smfr3d(i,k,j)   = 1.
             sh2o(i,k,j)     = 0.
             keepfr3dflag(i,k,j) = 0.
             tso(i,k,j) = min(271.4,tso(i,k,j))
           ENDDO

           DO k=1,nzs
              tso1d   (k) = tso(i,k,j)
           ENDDO

        else
        !-- land
           SEAICE(i,j)=0.
        endif ! land or sea ice
      ENDIF ! xland

      IF(SEAICE(I,J).GT.0.5)THEN
      !-- Sea-ice case
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
      print *,' IN LSMRUC_ICE ','its,ite,jts,jte,nzs,i,j', &
                its,ite,jts,jte,nzs, i,j
      print *,' SOILT,QVG,P8w',i,j,soilt(i,j),qvg(i,j),p8w(i,1,j)
      print *, 'LSMRUC_ICE, I,J,xland,seaice, QFX,HFX from SFCLAY',i,j,xland(i,j), &
                  seaice(i,j),qfx(i,j),hfx(i,j)
      print *, 'GSW, GLW =',i,j,gsw(i,j),glw(i,j)
      print *, 'SOILT, TSO start of time step =',soilt(i,j),(tso(i,k,j),k=1,nsl)
      print *, ' I,J=, after SFCLAY CHS,FLHC ',i,j,chs(i,j),flhc(i,j)
      print *, 'LSMRUC_ICE, ALB = ', alb(i,j),i,j
      print *, 'LSMRUC_ICE I,J,DT,RAINBL =',I,J,dt,RAINBL(i,j)
      print *, 'XLAND ---->, i,j',xland(i,j),i,j
    ENDIF

         TABS      = T3D(i,kms,j)
         QVATM     = QV3D(i,kms,j)
         QCATM     = QC3D(i,kms,j)
         PATM      = P8w(i,kms,j)*1.e-5
         !-- Z3D(1) is thickness between first full sigma level and the surface, 
         !-- but first mass level is at the half of the first sigma level 
         !-- (u and v are also at the half of first sigma level)
         CONFLX    = Z3D(i,kms,j)*0.5
         RHO       = RHO3D(I,kms,J)
         tso_dop   = tbot(i,j)

         !-- initialize snow, graupel and ice fractions in frozen precip
         PRCPMS = 0.
         newsnms = 0.
         prcpncliq = 0.
         prcpculiq = 0.
         prcpncfr = 0.
         prcpcufr = 0.

         snowrat = 0.
         grauprat = 0.
         icerat = 0.
         curat = 0.
       IF(FRPCPN) THEN
#if (EM_CORE==1)
         prcpncliq = rainncv(i,j)*(1.-frzfrac(i,j))
         prcpncfr = rainncv(i,j)*frzfrac(i,j)
         !- apply the same frozen precipitation fraction to convective precip
         !tgs - 31 mar17 - add safety temperature check in case Thompson MP produces
         !                 frozen precip at T > 273.
         if(frzfrac(i,j) > 0..and. tabs < 273.) then
           prcpculiq = max(0.,(rainbl(i,j)-rainncv(i,j))*(1.-frzfrac(i,j)))
           prcpcufr = max(0.,(rainbl(i,j)-rainncv(i,j))*frzfrac(i,j))
         else
           if(tabs < 273.) then
             prcpcufr = max(0.,(rainbl(i,j)-rainncv(i,j)))
             prcpculiq = 0.
           else
             prcpcufr = 0.
             prcpculiq = max(0.,(rainbl(i,j)-rainncv(i,j)))
           endif  ! tabs < 273.
         endif  ! frzfrac > 0.
         !--- 1*e-3 is to convert from mm/s to m/s
         PRCPMS   = (prcpncliq + prcpculiq)/DT*1.e-3
         NEWSNMS  = (prcpncfr + prcpcufr)/DT*1.e-3

         IF ( PRESENT( graupelncv ) ) THEN
             graupamt = graupelncv(i,j)
         ELSE
             graupamt = 0.
         ENDIF

         if((prcpncfr + prcpcufr) > 0.) then
         ! -- calculate snow, graupel and ice fractions in falling frozen precip
           snowrat=min(1.,max(0.,snowncv(i,j)/(prcpncfr + prcpcufr)))
           grauprat=min(1.,max(0.,graupamt/(prcpncfr + prcpcufr)))
           icerat=min(1.,max(0.,(prcpncfr-snowncv(i,j)-graupamt) &
                 /(prcpncfr + prcpcufr)))
           curat=min(1.,max(0.,(prcpcufr/(prcpncfr + prcpcufr))))
         endif
#else
         PRCPMS    = (RAINBL(i,j)/DT*1.e-3)*(1-FRZFRAC(I,J))
         NEWSNMS  = (RAINBL(i,j)/DT*1.e-3)*FRZFRAC(I,J)
         if(newsnms == 0.) then
           snowrat = 0.
         else
           snowrat = min(1.,newsnms/(newsnms+prcpms))
         endif
#endif

       ELSE  ! .not. FRPCPN
         if (tabs.le.273.15) then
           PRCPMS    = 0.
           NEWSNMS   = RAINBL(i,j)/DT*1.e-3
           !-- here no info about constituents of frozen precipitation,
           !-- suppose it is all snow
           snowrat = 1.
         else
           PRCPMS    = RAINBL(i,j)/DT*1.e-3
           NEWSNMS   = 0.
         endif
       ENDIF

         ! -- save time-step water equivalent of frozen precipitation 
         ! -- in PRECIPFR array to be used in module_diagnostics
         precipfr(i,j) = NEWSNMS * DT *1.e3

         !--- convert exchange coeff QKMS to [m/s]
         QKMS=FLQC(I,J)/RHO
         !TKMS=FLHC(I,J)/RHO/CP
         TKMS=FLHC(I,J)/RHO/(CP*(1.+0.84*QVATM))  ! mynnsfc uses CPM
         !--- convert incoming snow from mm to m
         SNWE=SNOW(I,J)*1.E-3
         SNHEI=SNOWH(I,J)

         SNOWFRAC=SNOWC(I,J)
         RHOSNFALL=RHOSNF(I,J)

         snowold(i,j)=snwe

         ! diagnose snow density from SWE and snow height
         if(SNOW(i,j).gt.0. .and. SNOWH(i,j).gt.0.) then
           RHOSN = SNOW(i,j)/SNOWH(i,j)
         else
           RHOSN = 300.
         endif

!#if ( NMM_CORE == 1 )
!     if(ktau+1.gt.1) then
!#else
!     if(ktau.gt.1) then
!#endif
      ! extract dew from the cloud water at the surface
      !30july13 QCG(I,J)=QCG(I,J)-DEW(I,J)/QKMS
      !endif

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
      print *,'SEA ICE, i,j,tso1d,PATM,TABS,QVATM,QCATM,RHO',  &
                        i,j,tso1d,PATM,TABS,QVATM,QCATM,RHO
      print *,'CONFLX =',CONFLX 
    ENDIF

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
      print *,'before SICETMP, spp_lsm, rstoch, field_sf_loc',      &
        i,j,spp_lsm,(rstoch(i,k,j),k=1,nzs),(field_sf_loc(i,k,j),k=1,nzs)
    ENDIF

!-----------------------------------------------------------------
        CALL SICETMP (spp_lsm,rstoch(i,:,j),field_sf_loc(i,:,j), & 
                dt,ktau,conflx,i,j,                              &
!--- input variables
                nzs,nddzs,                                       &   
                PRCPMS, NEWSNMS,SNWE,SNHEI,SNOWFRAC,             &
                RHOSN,RHONEWSN,RHOSNFALL,                        &
                snowrat,grauprat,icerat,curat,                   &
                PATM,TABS,QVATM,QCATM,RHO,                       &
                GLW(I,J),GSW(I,J),EMISSL(I,J),                   &
                QKMS,TKMS,                                       &
                alb(I,J),znt(I,J),snoalb(i,j),albbck(i,j),       &
                myjpbl,seaice(i,j),isice,                        &
!--- soil fixed fields
                zsmain,zshalf,DTDZS,DTDZS2,tbq,                  &
!--- constants
                cp,rovcp,g0,lv,stbolt,cw,c1sn,c2sn,              &
!--- output variables
                snweprint,snheiprint,rsm,tso1d,                  &
                soilt(I,J),soilt1(i,j),tsnav(i,j),dew(I,J),      &
                qvg(I,J),qsg(I,J),qcg(I,J),SMELT(I,J),           &
                SNOH(I,J),SNFLX(I,J),SNOM(I,J),SNOWFALLAC(I,J),  &
                ACSNOW(I,J),qfx(I,J),                            &
                lh(I,J),hfx(I,J),sflx(I,J),sublim(I,J),          &
                evapl(I,J),prcpl(I,J),budget(i,j),runoff1(i,j))

!-----------------------------------------------------------------

! Fill in field_sf to pass perturbed field of hydraulic cond. up to model driver and output
#if (EM_CORE==1)
       if (spp_lsm==1) then
         do k=1,nsl
           field_sf(i,k,j)=field_sf_loc(i,k,j)
         enddo
       endif
#endif

!***  DIAGNOSTICS
          !--- Convert the water unit into mm
          SFCRUNOFF(I,J) = SFCRUNOFF(I,J)+RUNOFF1(I,J)*DT*1000.0
          ACRUNOFF(I,J)  = ACRUNOFF(I,J)+RUNOFF1(I,J)*DT*1000.0

          do k=1,nzs
            tso(i,k,j) = min(271.4,tso1d(k))
          enddo

          !tgs add together dew and cloud at the ground surface
          !30july13        qcg(i,j)=qcg(i,j)+dew(i,j)/qkms

          Z0       (I,J) = ZNT (I,J)
          patmb=P8w(i,1,j)*1.e-2
          Q2SAT=QSN(TABS,TBQ)/PATMB
          QSFC(I,J) = QVG(I,J)/(1.+QVG(I,J))
          ! for MYJ PBL scheme
          IF((myjpbl).AND.(QVATM.GE.Q2SAT*0.95).AND.QVATM.LT.qvg(I,J))THEN
            CHKLOWQ(I,J)=0.
          ELSE
            CHKLOWQ(I,J)=1.
          ENDIF

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
      if(CHKLOWQ(I,J).eq.0.) then
         print *,'i,j,CHKLOWQ',  &
                  i,j,CHKLOWQ(I,J)
      endif
    ENDIF

          if(snow(i,j)==0.) EMISSL(i,j) = 0.98

! evan
          if (spp_lsm == 1) then
            EMISSL(i,j) = min(EMISSL(i,j) * (1. + 0.1*rstoch(i,1,j)), 1.)
          endif

          EMISS (I,J) = EMISSL(I,J)
          ! SNOW is in [mm], SNWE is in [m]; 
          SNOW   (i,j) = SNWE*1000.
          SNOWH  (I,J) = SNHEI 

      IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        print *,' LAND, I=,J=, QFX, HFX after SICETMP', i,j,lh(i,j),hfx(i,j)
      ENDIF
          SFCEVP (I,J) = SFCEVP (I,J) + QFX (I,J) * DT
          GRDFLX (I,J) = -1. * sflx(I,J)

          !--- SNOWC snow cover flag
          if(snowfrac > 0. .and. xice(i,j).ge.xice_threshold ) then
          ! sea ice - does not call FPE
            SNOWFRAC = SNOWFRAC*XICE(I,J)
          endif

          SNOWC(I,J)=SNOWFRAC

          !--- RHOSNF - density of snowfall
          RHOSNF(I,J)=RHOSNFALL

          ! Accumulated moisture flux [kg/m^2]
          SFCEVP (I,J) = SFCEVP (I,J) + QFX (I,J) * DT
 
          ac=0.
          as=0.
  
          as=max(0.,snwe-snowold(i,j))
          wb =rainbl(i,j)+smelt(i,j)*dt*1.e3   & ! source
                         -qfx(i,j)*dt          &
                         -runoff1(i,j)*dt*1.e3 &
                         -ac-as 

          waterbudget(i,j)=rainbl(i,j)+smelt(i,j)*dt*1.e3 & ! source
                         -qfx(i,j)*dt                     &
                         -runoff1(i,j)*dt*1.e3            &
                         -ac-as


          acwaterbudget(i,j)=acwaterbudget(i,j)+waterbudget(i,j)

        IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
          print *,'Budget',budget(i,j),i,j
          print *,'Water budget ', i,j,waterbudget(i,j)
          print *,'rainbl,qfx*dt,runoff1,smelt*dt*1.e3',          &
                i,j,rainbl(i,j),qfx(i,j)*dt,runoff1(i,j)*dt*1.e3, &
                smelt(i,j)*dt*1.e3
          print *,'SNOW,SNOWold',i,j,snwe,snowold(i,j)
          print *,'SNOW-SNOWold',i,j,max(0.,snwe-snowold(i,j))
          print *,'SEA ICE, i,j,tso1d,soilt - end of time step',  &
                            i,j,tso1d,soilt(i,j)
          print *,'SEA ICE, QFX, HFX after SICETMP', i,j,lh(i,j),hfx(i,j)
        ENDIF

      ENDIF ! sea ice

2999  continue ! lakes

      ENDDO

   ENDDO

!-----------------------------------------------------------------
   END SUBROUTINE LSMRUC_ice
!-----------------------------------------------------------------



   SUBROUTINE SICETMP (spp_lsm,rstochcol,fieldcol_sf,            &
                delt,ktau,conflx,i,j,                            &
!--- input variables
                nzs,nddzs,                                       &
                PRCPMS,NEWSNMS,SNWE,SNHEI,SNOWFRAC,              &
                RHOSN,RHONEWSN,RHOSNFALL,                        &
                snowrat,grauprat,icerat,curat,                   &
                PATM,TABS,QVATM,QCATM,rho,                       &
                GLW,GSW,EMISS,QKMS,TKMS,                         &
                ALB,ZNT,ALB_SNOW,ALB_SNOW_FREE,                  &
                MYJ,SEAICE,ISICE,                                &
!--- soil fixed fields
                zsmain,zshalf,DTDZS,DTDZS2,tbq,                  &
!--- constants
                cp,rovcp,g0,lv,stbolt,cw,c1sn,c2sn,              &
!--- output variables
                snweprint,snheiprint,rsm,                        &
                ts1d,soilt,soilt1,tsnav,dew,qvg,qsg,qcg,         &
                SMELT,SNOH,SNFLX,SNOM,SNOWFALLAC,ACSNOW,         &
                eeta,qfx,hfx,s,sublim,                           &
                evapl,prcpl,fltot,runoff1                        &
                 )
!-----------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------

!--- input variables

   INTEGER,  INTENT(IN   )   ::  isice,i,j,ktau,nzs ,      &
                                 nddzs                             !nddzs=2*(nzs-2)

   REAL,     INTENT(IN   )   ::  DELT,CONFLX
   REAL,     INTENT(IN   )   ::  C1SN,C2SN,CW

   LOGICAL,    INTENT(IN   ) ::     myj
!--- 3-D Atmospheric variables
   REAL                                                        , &
            INTENT(IN   )    ::                            PATM, &
                                                           TABS, &
                                                          QVATM, &
                                                          QCATM
!--- 2-D variables
   REAL                                                        , &
            INTENT(IN   )    ::                             GLW, &
                                                            GSW, &
                                                  ALB_SNOW_FREE, &
                                                         SEAICE, &
                                                            RHO, &
                                                           QKMS, &
                                                           TKMS
                                                             
   REAL                                                        , &
            INTENT(INOUT)    ::                           EMISS, &
                                                       SNOWFRAC, &
                                                       ALB_SNOW, &
                                                            ALB

!--- constants
   REAL,     INTENT(IN   )   ::                                  &
                                                             CP, &
                                                          ROVCP, &
                                                             G0, &
                                                             LV, &
                                                         STBOLT

   REAL,     DIMENSION(1:NZS), INTENT(IN)  ::            ZSMAIN, &
                                                         ZSHALF, &
                                                         DTDZS2 

   REAL,     DIMENSION(1:NZS), INTENT(IN)  ::          rstochcol
   REAL,     DIMENSION(1:NZS), INTENT(INOUT) ::        fieldcol_sf


   REAL,     DIMENSION(1:NDDZS), INTENT(IN)  ::           DTDZS

   REAL,     DIMENSION(1:5001), INTENT(IN)  ::              TBQ

!--- input/output variables
   REAL,     DIMENSION( 1:nzs )                                , &
             INTENT(INOUT)   ::                            TS1D

!-------- 2-d variables
   REAL                                                        , &
             INTENT(INOUT)   ::                             DEW, &
                                                           EETA, &
                                                          EVAPL, &
                                                          RHOSN, & 
                                                       RHONEWSN, &
                                                      rhosnfall, &
                                                        snowrat, &
                                                       grauprat, &
                                                         icerat, &
                                                          curat, &
                                                         SUBLIM, &
                                                          PRCPL, &
                                                            QVG, &
                                                            QSG, &
                                                            QCG, &
                                                            QFX, &
                                                            HFX, &
                                                          fltot, &
                                                              S, &  
                                                        RUNOFF1, &
                                                         ACSNOW, &
                                                     SNOWFALLAC, &
                                                           SNWE, &
                                                          SNHEI, &
                                                          SMELT, &
                                                           SNOM, &
                                                           SNOH, &
                                                          SNFLX, &
                                                          SOILT, &
                                                         SOILT1, &
                                                          TSNAV, &
                                                            ZNT

   REAL,     DIMENSION(1:NZS)              ::                    &
                                                           tice, &
                                                        rhosice, &
                                                         capice, &
                                                       thdifice, &
                                                          TS1DS
!-------- 1-d variables
   REAL :: &
                                                            DEWS, &
                                                           EETAs, &
                                                          EVAPLs, &
                                                          PRCPLS, &
                                                            QVGS, &
                                                            QSGS, &
                                                            QCGS, &
                                                            QFXS, &
                                                            HFXS, &
                                                          fltots, &
                                                        RUNOFF1S, &
                                                              SS, &
                                                          SOILTs
                     
   REAL,  INTENT(INOUT)                     ::              RSM, &  
                                                      SNWEPRINT, &
                                                     SNHEIPRINT
   INTEGER,   INTENT(IN)                    ::          spp_lsm     
!--- Local variables

   INTEGER ::  K,ILNB

   REAL    ::  BSN, XSN                                        , &
               RAINF, SNTH, NEWSN, PRCPMS, NEWSNMS             , &
               T3, UPFLUX, XINET
   REAL    ::  snhei_crit, snhei_crit_newsn, keep_snow_albedo, SNOWFRACnewsn
   REAL    ::  newsnowratio, dd1

   REAL    ::  rhonewgr,rhonewice

   REAL    ::  RNET,GSWNEW,GSWIN,EMISSN,ZNTSN,EMISS_snowfree
   REAL    ::  snow_mosaic, snfr
   real    ::  cice, albice, albsn
!-----------------------------------------------------------------
        integer,   parameter      ::      ilsnow=99 
        
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        print *,' in SICETMP',i,j,nzs,nddzs,      &
                 SNWE,RHOSN,SNOM,SMELT,TS1D
    ENDIF

        snow_mosaic=0.
        snfr = 1.
        NEWSN=0.
        newsnowratio = 0.
        snowfracnewsn=0.
        if(snhei == 0.) snowfrac=0.
        smelt = 0.
        RAINF = 0.
        RSM=0.

!---initialize local arrays for sea ice
        do k=1,nzs
          tice(k) = 0.
          rhosice(k) = 0. 
          cice = 0.
          capice(k) = 0.
          thdifice(k) = 0.
        enddo

        GSWnew=GSW
        GSWin=GSW/(1.-alb)
        ALBice=ALB_SNOW_FREE
        ALBsn=alb_snow
        EMISSN = 0.98
        EMISS_snowfree = 0.98

! evan
        if (spp_lsm == 1) then
          EMISS_snowfree = min(EMISS_snowfree * (1. + 0.1*rstochcol(1)), 1.)
        endif

        !--- sea ice properties
        !--- N.N Zubov "Arctic Ice"
        !--- no salinity dependence because we consider the ice pack
        !--- to be old and to have low salinity (0.0002)
        do k=1,nzs
          tice(k) = ts1d(k) - 273.15
          rhosice(k) = 917.6/(1-0.000165*tice(k))
          cice = 2115.85 +7.7948*tice(k)
          capice(k) = cice*rhosice(k)
          thdifice(k) = 2.260872/capice(k)
         enddo
        !-- SEA ICE ALB dependence on ice temperature. When ice temperature is
        !-- below critical value of -10C - no change to albedo.
        !-- If temperature is higher that -10C then albedo is decreasing.
        !-- The minimum albedo at t=0C for ice is 0.1 less.
        ALBice = MIN(ALB_SNOW_FREE,MAX(ALB_SNOW_FREE - 0.05,   &
                 ALB_SNOW_FREE - 0.1*(tice(1)+10.)/10. ))

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        print *,'alb_snow_free',ALB_SNOW_FREE
        print *,'GSW,GSWnew,GLW,SOILT,EMISS,ALB,ALBice,SNWE',&
                 GSW,GSWnew,GLW,SOILT,EMISS,ALB,ALBice,SNWE
    ENDIF

        if(snhei.gt.0.0081*1.e3/rhosn) then
        !*** Update snow density for current temperature (Koren et al. 1999)
          BSN=delt/3600.*c1sn*exp(0.08*min(0.,tsnav)-c2sn*rhosn*1.e-3)
            if(bsn*snwe*100..lt.1.e-4) goto 777
          XSN=rhosn*(exp(bsn*snwe*100.)-1.)/(bsn*snwe*100.)
          rhosn=MIN(MAX(58.8,XSN),500.)
 777    continue
        endif

        IF(PRCPMS > 0.) THEN
        ! PRCPMS is liquid precipitation rate
        ! RAINF is a flag used for calculation of rain water
        ! heat content contribution into heat budget equation. Rain's temperature
        ! is set equal to air temperature at the first atmospheric
        ! level.  
           RAINF=1.
       ENDIF

       !---- ACSNOW - run-total snowfall water [mm]
       ! acsnow=acsnow+newsn*1.e3

       newsn=newsnms*delt ! [m]
       IF(NEWSN.GT.0.) THEN

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
      print *, 'THERE IS NEW SNOW, newsn', newsn
    ENDIF

       newsnowratio = min(1.,newsn/(snwe+newsn))

         !*** Calculate fresh snow density (t > -15C, else MIN value)
         !--- old formulation from Koren (1999)
         !*** Eq. 10 from Koren et al. (1999)
         !if(tabs.lt.258.15) then
         !  rhonewsn=50.
         !  rhonewsn=100.
         !  rhonewsn=62.5
         !else
         !  rhonewsn=MIN(rhonewsn,400.)
         !endif
         !--- end of old formulation

         !--- new formulation
         !--- 27 Feb 2014 - empirical formulations from John M. Brown
         !rhonewsn=min(250.,rhowater/max(4.179,(13.*tanh((274.15-Tabs)*0.3333))))
         !--- 13 Mar 2018 - formulation from Trevor Alcott
         rhonewsn=min(125.,1000.0/max(8.,(17.*tanh((276.65-Tabs)*0.15))))
         rhonewgr=min(500.,rhowater/max(2.,(3.5*tanh((274.15-Tabs)*0.3333))))
         rhonewice=rhonewsn

         !--- compute density of "snowfall" from weighted contribution
         !                 of snow, graupel and ice fractions

         !13mar18-- rhosnfall = min(500.,max(76.9,(rhonewsn*snowrat +  &
         rhosnfall = min(500.,max(58.8,(rhonewsn*snowrat +  &
                     rhonewgr*grauprat + rhonewice*icerat + rhonewgr*curat)))

         ! from now on rhonewsn is the density of falling frozen precipitation
         rhonewsn=rhosnfall

         !*** Define average snow density of the snow pack considering
         !*** the amount of fresh snow (eq. 9 in Koren et al.(1999) 
         !*** without snow melt )
         xsn=(rhosn*snwe+rhonewsn*newsn)/                         &
             (snwe+newsn)
         rhosn=MIN(MAX(58.8,XSN),500.)

         !Update snow on the ground
         snwe=max(0.,snwe+newsn)
         snhei=snwe*rhowater/rhosn
         NEWSN=NEWSN*rhowater/rhonewsn ! [m] of snow depth

         ! fraction of new snow is needed only for albedo
         SNOWFRACnewsn=MIN(1.,SNHEI/SNHEI_CRIT_newsn)
       ENDIF ! end NEWSN > 0.


         ! SNHEI_CRIT is a threshold for fractional snow
         SNHEI_CRIT=0.01601*1.e3/rhosn
         SNHEI_CRIT_newsn=0.0005*1.e3/rhosn
         ! snowfrac from the previous time step
         SNOWFRAC=MIN(1.,SNHEI/(2.*SNHEI_CRIT))

     IF(SNHEI.GT.0.0) THEN
     !-- SNOW on ice including the fresh snow

         if(snowfrac < 0.75) snow_mosaic = 1.
       
         KEEP_SNOW_ALBEDO = 0.
       IF (NEWSN > 0. .and. snowfracnewsn > 0.99) THEN
         ! new snow
         KEEP_SNOW_ALBEDO = 1.
         snow_mosaic=0.  ! ???
       ENDIF

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
      print *,'SNHEI_CRIT,SNOWFRAC,SNHEI_CRIT_newsn,SNOWFRACnewsn', &
               SNHEI_CRIT,SNOWFRAC,SNHEI_CRIT_newsn,SNOWFRACnewsn
    ENDIF

     if( snow_mosaic == 1.) then
     !---- snow mosaic
       ALBsn=alb_snow
       Emiss= emissn
     else
       ALBsn   = MAX(keep_snow_albedo*alb_snow,               &
                 MIN((albice + (alb_snow - albice) * snowfrac), alb_snow))
       Emiss   = MAX(keep_snow_albedo*emissn,                 &
                 MIN((emiss_snowfree +                        &
                (emissn - emiss_snowfree) * snowfrac), emissn))
     endif

     IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'Snow on ice snow_mosaic,ALBsn,emiss',i,j,ALBsn,emiss,snow_mosaic
     ENDIF

     !-- ALB dependence on snow temperature. When snow temperature is
     !-- below critical value of -10C - no change to albedo.
     !-- If temperature is higher that -10C then albedo is decreasing.
     if(albsn.lt.alb_snow .or. keep_snow_albedo .eq.1.)then
        ALB=ALBsn
     else
     !-- change albedo when no fresh snow
        ALB = MIN(ALBSN,MAX(ALBSN - 0.15*ALBSN*(soilt - 263.15)/  &
                 (273.15-263.15), ALBSN - 0.1))
     endif

     if (snow_mosaic==1.) then 
     ! may 2014 - treat separately snow-free ice

     ! compute absorbed GSW for snow-free portion
         gswnew=GSWin*(1.-albice)
     !--------------
         T3      = STBOLT*SOILT*SOILT*SOILT
         UPFLUX  = T3 *SOILT
         XINET   = EMISS_snowfree*(GLW-UPFLUX)
         RNET    = GSWnew + XINET
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
     print *,'Fractional snow - snowfrac=',snowfrac
     print *,'Snowfrac<1 GSWin,GSWnew -',GSWin,GSWnew,'SOILT, RNET',soilt,rnet
    ENDIF
         do k=1,nzs
           ts1ds(k) = ts1d(k)
         enddo
         soilts = soilt
         qvgs = qvg
         qsgs = qsg
         qcgs = qcg
         smelt= 0.
         runoff1s=0.
 
         CALL SICE(                                             &
!--- input variables
            i,j,delt,ktau,conflx,nzs,nddzs,                     &
            PRCPMS,RAINF,PATM,QVATM,QCATM,GLW,GSWnew,           &
            0.98,RNET,QKMS,TKMS,rho,myj,                        &
!--- sea ice parameters
            tice,rhosice,capice,thdifice,                       &
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &
!--- constants
            lv,CP,rovcp,cw,stbolt,tabs,                         &
!--- output variable
            ts1ds,dews,soilts,qvgs,qsgs,qcgs,                   &
            eetas,qfxs,hfxs,ss,evapls,prcpls,fltots             &
                                                                )
         runoff1 = prcpms

     endif ! snow_mosaic=1.
                       
       !--- recompute absorbed solar radiation and net radiation
       !--- for updated value of snow albedo - ALB
       gswnew=GSWin*(1.-alb)
       !print *,'SNOW fraction GSWnew',gswnew,'alb=',alb
       !--------------
       T3      = STBOLT*SOILT*SOILT*SOILT
       UPFLUX  = T3 *SOILT
       XINET   = EMISS*(GLW-UPFLUX)
       RNET    = GSWnew + XINET
      IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        print *,'RNET=',rnet
        print *,'SNOW - I,J,newsn,snwe,snhei,GSW,GSWnew,GLW,UPFLUX,ALB',&
                 i,j,newsn,snwe,snhei,GSW,GSWnew,GLW,UPFLUX,ALB
      ENDIF

       !--- treat snow-covered ice
       if(snow_mosaic==1.)then
         snfr=1.
       else
         snfr=snowfrac
       endif

         CALL SNOWSEAICE (                                      &
            i,j,delt,ktau,conflx,nzs,nddzs,                     &    
            rhonewsn,SNHEI_CRIT,                                &
            PRCPMS,RAINF,NEWSN,snhei,SNWE,snfr,                 &    
            RHOSN,PATM,QVATM,QCATM,                             &    
            GLW,GSWnew,EMISS,RNET,                              &    
            QKMS,TKMS,RHO,myj,                                  &    
!--- sea ice parameters
            ALB,ZNT,                                            &
            tice,rhosice,capice,thdifice,                       &    
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &    
!--- constants
            lv,CP,rovcp,cw,stbolt,tabs,                         &    
!--- output variables
            ilnb,snweprint,snheiprint,rsm,ts1d,                 &    
            dew,soilt,soilt1,tsnav,qvg,qsg,qcg,                 &    
            SMELT,SNOH,SNFLX,SNOM,eeta,                         &    
            qfx,hfx,s,sublim,prcpl,fltot                        &    
                                                                )    
           runoff1 = smelt

       if(snhei.eq.0.) then
       !--- all snow is melted
       alb=alb_snow_free
       endif

       if (snow_mosaic==1.) then
       !--- Now combine fluxes for snow-free sea ice and snow-covered area
     IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'SOILT snow on ice', soilt
     ENDIF
        do k=1,nzs
          ts1d(k) = ts1ds(k)*(1.-snowfrac) + ts1d(k)*snowfrac
        enddo
          dew = dews*(1.-snowfrac) + dew*snowfrac
          soilt = soilts*(1.-snowfrac) + soilt*snowfrac
          qvg = qvgs*(1.-snowfrac) + qvg*snowfrac
          qsg = qsgs*(1.-snowfrac) + qsg*snowfrac
          qcg = qcgs*(1.-snowfrac) + qcg*snowfrac
          eeta = eetas*(1.-snowfrac) + eeta*snowfrac
          qfx = qfxs*(1.-snowfrac) + qfx*snowfrac
          hfx = hfxs*(1.-snowfrac) + hfx*snowfrac
          s = ss*(1.-snowfrac) + s*snowfrac
          sublim = eeta
          prcpl = prcpls*(1.-snowfrac) + prcpl*snowfrac
          fltot = fltots*(1.-snowfrac) + fltot*snowfrac
          !alb
          ALB   = MAX(keep_snow_albedo*alb,              &
                  MIN((albice + (alb - alb_snow_free) * snowfrac), alb))

          Emiss = MAX(keep_snow_albedo*emissn,           &
                  MIN((emiss_snowfree +                  &
              (emissn - emiss_snowfree) * snowfrac), emissn))

          runoff1 = runoff1s*(1.-snowfrac) + runoff1*snowfrac
          smelt = smelt * snowfrac
          snoh = snoh * snowfrac
          snflx = snflx * snowfrac
          snom = snom * snowfrac
      IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        print *,'SOILT combined on ice', soilt
      ENDIF
    endif ! snow_mosaic = 1.
 
      !  run-total accumulated snow based on snowfall and snowmelt in [m]
      snowfallac = snowfallac + max(0.,(newsn - rhowater/rhonewsn*smelt*delt*newsnowratio))

   ELSE
         !--- no snow
         snheiprint=0.
         snweprint=0.
         smelt=0.

         T3      = STBOLT*SOILT*SOILT*SOILT
         UPFLUX  = T3 *SOILT
         XINET   = EMISS*(GLW-UPFLUX)
         RNET    = GSWnew + XINET
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
     print *,'NO snow on the ground GSWnew -',GSWnew,'RNET=',rnet
    ENDIF

         ! If current ice albedo is not the same as from the previous time step, then
         ! update GSW, ALB and RNET for surface energy budget
         if(ALB.ne.ALBice) GSWnew=GSW/(1.-ALB)*(1.-ALBice)
         alb=albice
         RNET    = GSWnew + XINET

          CALL SICE(                                            &
!--- input variables
            i,j,delt,ktau,conflx,nzs,nddzs,   &
            PRCPMS,RAINF,PATM,QVATM,QCATM,GLW,GSWnew,           &
            EMISS,RNET,QKMS,TKMS,rho,myj,                       &
!--- sea ice parameters
            tice,rhosice,capice,thdifice,                       &
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &
!--- constants
            lv,CP,rovcp,cw,stbolt,tabs,                         &
!--- output variables
            ts1d,dew,soilt,qvg,qsg,qcg,                         &
            eeta,qfx,hfx,s,evapl,prcpl,fltot                          &
                                                                )
           runoff1 = prcpms

        ENDIF

!---------------------------------------------------------------
   END SUBROUTINE SICETMP
!---------------------------------------------------------------

        SUBROUTINE SICE (                                       &
!--- input variables
            i,j,delt,ktau,conflx,nzs,nddzs,                     &
            PRCPMS,RAINF,PATM,QVATM,QCATM,GLW,GSW,              &
            EMISS,RNET,QKMS,TKMS,rho,myj,                       &
!--- sea ice parameters
            tice,rhosice,capice,thdifice,                       &
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &
!--- constants
            xlv,CP,rovcp,cw,stbolt,tabs,                        &
!--- output variables
            tso,dew,soilt,qvg,qsg,qcg,                          &
            eeta,qfx,hfx,s,evapl,prcpl,fltot                    &
                                                                )

!*****************************************************************
!   Energy budget and  heat diffusion eqns. for
!   sea ice
!*************************************************************

        IMPLICIT NONE
!-----------------------------------------------------------------

!--- input variables

   INTEGER,  INTENT(IN   )   ::  ktau,nzs                , &
                                 nddzs                    !nddzs=2*(nzs-2)
   INTEGER,  INTENT(IN   )   ::  i,j
   REAL,     INTENT(IN   )   ::  DELT,CONFLX
   LOGICAL,  INTENT(IN   )   ::  myj
!--- 3-D Atmospheric variables
   REAL,                                                         &
            INTENT(IN   )    ::                            PATM, &
                                                          QVATM, &
                                                          QCATM
!--- 2-D variables
   REAL,                                                         &
            INTENT(IN   )    ::                             GLW, &
                                                            GSW, &
                                                          EMISS, &
                                                            RHO, &
                                                           QKMS, &
                                                           TKMS
!--- sea ice properties
   REAL,    DIMENSION(1:NZS)                                   , &
            INTENT(IN   )    ::                                  &
                                                           tice, &
                                                        rhosice, &
                                                         capice, &
                                                       thdifice


   REAL,     INTENT(IN   )   ::                                  &
                                                             CW, &
                                                            XLV


   REAL,     DIMENSION(1:NZS), INTENT(IN)  ::            ZSMAIN, &
                                                         ZSHALF, &
                                                         DTDZS2

   REAL,     DIMENSION(1:NDDZS), INTENT(IN)  ::           DTDZS

   REAL,     DIMENSION(1:5001), INTENT(IN)  ::              TBQ


!--- input/output variables
!----soil temperature
   REAL,     DIMENSION( 1:nzs ),  INTENT(INOUT)   ::        TSO
!-------- 2-d variables
   REAL,                                                         &
             INTENT(INOUT)   ::                             DEW, &
                                                           EETA, &
                                                          EVAPL, &
                                                          PRCPL, &
                                                            QVG, &
                                                            QSG, &
                                                            QCG, &
                                                           RNET, &
                                                            QFX, &
                                                            HFX, &
                                                              S, &
                                                          SOILT

!--- Local variables
   REAL    ::  x,x1,x2,x4,tn,denom
   REAL    ::  RAINF,  PRCPMS                                  , &
               TABS, T3, UPFLUX, XINET

   REAL    ::  CP,rovcp,G0,LV,STBOLT,xlmelt,dzstop             , &
               epot,fltot,ft,fq,hft,ras,cvw                    

   REAL    ::  FKT,D1,D2,D9,D10,DID,R211,R21,R22,R6,R7,D11     , &
               PI,H,FKQ,R210,AA,BB,PP,Q1,QS1,TS1,TQ2,TX2       , &
               TDENOM,QGOLD,SNOH

   REAL    ::  AA1,RHCS, icemelt


   REAL,     DIMENSION(1:NZS)  ::   cotso,rhtso

   INTEGER ::  nzs1,nzs2,k,k1,kn,kk

!-----------------------------------------------------------------

        !-- define constants
        !STBOLT=5.670151E-8
        XLMELT=3.35E+5
        cvw=cw

        prcpl=prcpms

        NZS1=NZS-1
        NZS2=NZS-2
        dzstop=1./(zsmain(2)-zsmain(1))
        RAS=RHO*1.E-3

        do k=1,nzs
          cotso(k)=0.
          rhtso(k)=0.
        enddo

        cotso(1)=0.
        rhtso(1)=TSO(NZS)

        DO K=1,NZS2
          KN=NZS-K
          K1=2*KN-3
          X1=DTDZS(K1)*THDIFICE(KN-1)
          X2=DTDZS(K1+1)*THDIFICE(KN)
          FT=TSO(KN)+X1*(TSO(KN-1)-TSO(KN))                             &
             -X2*(TSO(KN)-TSO(KN+1))
          DENOM=1.+X1+X2-X2*cotso(K)
          cotso(K+1)=X1/DENOM
          rhtso(K+1)=(FT+X2*rhtso(K))/DENOM
        ENDDO

!************************************************************************
!--- THE HEAT BALANCE EQUATION (Smirnova et al., 1996, EQ. 21,26)
        RHCS=CAPICE(1)
        H=1.
        FKT=TKMS
        D1=cotso(NZS1)
        D2=rhtso(NZS1)
        TN=SOILT
        D9=THDIFICE(1)*RHCS*dzstop
        D10=TKMS*CP*RHO
        R211=.5*CONFLX/DELT
        R21=R211*CP*RHO
        R22=.5/(THDIFICE(1)*DELT*dzstop**2)
        R6=EMISS *STBOLT*.5*TN**4
        R7=R6/TN
        D11=RNET+R6
        TDENOM=D9*(1.-D1+R22)+D10+R21+R7                              &
              +RAINF*CVW*PRCPMS
        FKQ=QKMS*RHO
        R210=R211*RHO
        AA=XLS*(FKQ+R210)/TDENOM
        BB=(D10*TABS+R21*TN+XLS*(QVATM*FKQ                            &
          +R210*QVG)+D11+D9*(D2+R22*TN)                               &
          +RAINF*CVW*PRCPMS*max(273.15,TABS)                          &
           )/TDENOM
        AA1=AA
        PP=PATM*1.E3
        AA1=AA1/PP
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        PRINT *,' VILKA-SEAICE'
        print *,'D10,TABS,R21,TN,QVATM,FKQ',                          &
                 D10,TABS,R21,TN,QVATM,FKQ
        print *,'RNET, EMISS, STBOLT, SOILT',RNET, EMISS, STBOLT, SOILT
        print *,'R210,QVG,D11,D9,D2,R22,RAINF,CVW,PRCPMS,TDENOM',     &
                 R210,QVG,D11,D9,D2,R22,RAINF,CVW,PRCPMS,TDENOM
        print *,'tn,aa1,bb,pp,fkq,r210',                              &
                 tn,aa1,bb,pp,fkq,r210
    ENDIF

        QGOLD=QSG
        CALL VILKA(TN,AA1,BB,PP,QS1,TS1,TBQ,KTAU,i,j)
        !--- it is saturation over sea ice
        QVG=QS1
        QSG=QS1
        TSO(1)=min(271.4,TS1)
        QCG=0.
        !--- sea ice melting is not included in this simple approach
        !--- SOILT - skin temperature
        SOILT=TSO(1)
        !---- Final solution for soil temperature - TSO
        DO K=2,NZS
          KK=NZS-K+1
          TSO(K)=min(271.4,rhtso(KK)+cotso(KK)*TSO(K-1))
        END DO
        !--- CALCULATION OF DEW USING NEW VALUE OF QSG OR TRANSP IF NO DEW
        DEW=0.

        !--- THE DIAGNOSTICS OF SURFACE FLUXES 
        T3      = STBOLT*TN*TN*TN
        UPFLUX  = T3 *0.5*(TN+SOILT)
        XINET   = EMISS*(GLW-UPFLUX)
        !HFT=-TKMS*CP*RHO*(TABS-SOILT)
        HFX=-TKMS*CP*RHO*(TABS-SOILT)                        &
               *(P1000mb*0.00001/Patm)**ROVCP
        Q1=-QKMS*RAS*(QVATM - QSG)
        IF (Q1.LE.0.) THEN
        ! ---  condensation
          if(myj) then
          !-- moisture flux for coupling with MYJ PBL
            EETA=-QKMS*RAS*(QVATM/(1.+QVATM) - QSG/(1.+QSG))*1.E3
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'MYJ EETA',eeta
    ENDIF
          else ! myj
          !-- actual moisture flux from RUC LSM
            DEW=QKMS*(QVATM-QSG)
            EETA= - RHO*DEW
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'RUC LSM EETA',eeta
    ENDIF
          endif ! myj
          QFX= XLS*EETA
          EETA= - RHO*DEW
        ELSE
        ! ---  evaporation
          if(myj) then
          !-- moisture flux for coupling with MYJ PBL
            EETA=-QKMS*RAS*(QVATM/(1.+QVATM) - QVG/(1.+QVG))*1.E3
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'MYJ EETA',eeta
    ENDIF
          else ! myj
          ! to convert from m s-1 to kg m-2 s-1: *rho water=1.e3************
          !-- actual moisture flux from RUC LSM
            EETA = Q1*1.E3
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'RUC LSM EETA',eeta
    ENDIF
          endif ! myj

          QFX= XLS * EETA
          EETA = Q1*1.E3
        ENDIF

        EVAPL=EETA

        S=THDIFICE(1)*CAPICE(1)*DZSTOP*(TSO(1)-TSO(2))
        ! heat storage in surface layer
        SNOH=0.
        ! There is ice melt
        X= (cp*rho*r211+rhcs*zsmain(2)*0.5/delt)*(SOILT-TN) +   &
            XLS*rho*r211*(QSG-QGOLD)
        X=X &
        ! "heat" from rain
          -RAINF*CVW*PRCPMS*(max(273.15,TABS)-SOILT)

        !-- excess energy spent on sea ice melt
        icemelt=RNET-XLS*EETA -HFT -S -X
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        print *,'icemelt=',icemelt
    ENDIF

        FLTOT=RNET-XLS*EETA-HFT-S-X-icemelt
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'SICE - FLTOT,RNET,HFT,QFX,S,SNOH,X=', &
                       FLTOT,RNET,HFT,XLS*EETA,s,icemelt,X
    ENDIF

!-------------------------------------------------------------------
   END SUBROUTINE SICE
!-------------------------------------------------------------------

           SUBROUTINE SNOWSEAICE(                               &
            i,j,delt,ktau,conflx,nzs,nddzs,                     &
            rhonewsn,SNHEI_CRIT,                                &
            PRCPMS,RAINF,NEWSNOW,snhei,SNWE,snowfrac,           &
            RHOSN,PATM,QVATM,QCATM,                             &
            GLW,GSW,EMISS,RNET,                                 &
            QKMS,TKMS,RHO,myj,                                  &
!--- sea ice parameters
            ALB,ZNT,                                            &
            tice,rhosice,capice,thdifice,                       &
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &
!--- constants
            xlv,CP,rovcp,cw,stbolt,tabs,                        &
!--- output variables
            ilnb,snweprint,snheiprint,rsm,tso,                  &
            dew,soilt,soilt1,tsnav,qvg,qsg,qcg,                 &
            SMELT,SNOH,SNFLX,SNOM,eeta,                         &
            qfx,hfx,s,sublim,prcpl,fltot                        &
                                                                )
!***************************************************************
!   Solving energy budget for snow on sea ice and heat diffusion 
!   eqns. in snow and sea ice
!***************************************************************


        IMPLICIT NONE
!-------------------------------------------------------------------
!--- input variables

   INTEGER,  INTENT(IN   )   ::  ktau,nzs     ,                  &
                                 nddzs                         !nddzs=2*(nzs-2)
   INTEGER,  INTENT(IN   )   ::  i,j

   REAL,     INTENT(IN   )   ::  DELT,CONFLX,PRCPMS            , &
                                 RAINF,NEWSNOW,RHONEWSN,         &
                                 snhei_crit
   real                      ::  rhonewcsn

   LOGICAL,  INTENT(IN   )   ::  myj
!--- 3-D Atmospheric variables
   REAL,                                                         &
            INTENT(IN   )    ::                            PATM, &
                                                          QVATM, &
                                                          QCATM
!--- 2-D variables
   REAL                                                        , &
            INTENT(IN   )    ::                             GLW, &
                                                            GSW, &
                                                            RHO, &
                                                           QKMS, &
                                                           TKMS

!--- sea ice properties
   REAL,     DIMENSION(1:NZS)                                  , &
            INTENT(IN   )    ::                                  &
                                                           tice, &
                                                        rhosice, &
                                                         capice, &
                                                       thdifice

   REAL,     INTENT(IN   )   ::                                  &
                                                             CW, &
                                                            XLV

   REAL,     DIMENSION(1:NZS), INTENT(IN)  ::            ZSMAIN, &
                                                         ZSHALF, &
                                                         DTDZS2

   REAL,     DIMENSION(1:NDDZS), INTENT(IN)  ::           DTDZS

   REAL,     DIMENSION(1:5001), INTENT(IN)  ::              TBQ

!--- input/output variables
!-------- 3-d soil moisture and temperature
   REAL,     DIMENSION(  1:nzs )                               , &
             INTENT(INOUT)   ::                             TSO

!-------- 2-d variables
   REAL                                                        , &
             INTENT(INOUT)   ::                             DEW, &
                                                           EETA, &
                                                          RHOSN, &
                                                         SUBLIM, &
                                                          PRCPL, &
                                                            ALB, &
                                                          EMISS, &
                                                            ZNT, &
                                                            QVG, &
                                                            QSG, &
                                                            QCG, &
                                                            QFX, &
                                                            HFX, &
                                                              S, &
                                                           SNWE, &
                                                          SNHEI, &
                                                          SMELT, &
                                                           SNOM, &
                                                           SNOH, &
                                                          SNFLX, &
                                                          SOILT, &
                                                         SOILT1, &
                                                       SNOWFRAC, &
                                                          TSNAV

   INTEGER, INTENT(INOUT)    ::                            ILNB

   REAL,     INTENT(OUT)                    ::              RSM, &
                                                      SNWEPRINT, &
                                                     SNHEIPRINT
!--- Local variables


   INTEGER ::  nzs1,nzs2,k,k1,kn,kk
   REAL    ::  x,x1,x2,dzstop,ft,tn,denom

   REAL    ::  SNTH, NEWSN                                     , &
               TABS, T3, UPFLUX, XINET                         , &
               BETA, SNWEPR,EPDT,PP
   REAL    ::  CP,rovcp,G0,LV,xlvm,STBOLT,xlmelt               , &
               epot,fltot,fq,hft,q1,ras,rhoice,ci,cvw          , &
               RIW,DELTSN,H

   REAL    ::  rhocsn,thdifsn,                                   &
               xsn,ddzsn,x1sn,d1sn,d2sn,d9sn,r22sn

   REAL    ::  cotsn,rhtsn,xsn1,ddzsn1,x1sn1,ftsnow,denomsn
   REAL    ::  fso,fsn,                                          &
               FKT,D1,D2,D9,D10,DID,R211,R21,R22,R6,R7,D11,      &
               FKQ,R210,AA,BB,QS1,TS1,TQ2,TX2,                   &
               TDENOM,AA1,RHCS,H1,TSOB, SNPRIM,                  &
               SNODIF,SOH,TNOLD,QGOLD,SNOHGNEW
   REAL,     DIMENSION(1:NZS)  ::  cotso,rhtso

   REAL                   :: RNET,rsmfrac,soiltfrac,hsn,icemelt,rr
   integer                ::      nmelt


!-----------------------------------------------------------------
        XLMELT=3.35E+5
        !-- heat of sublimation of water vapor
        XLVm=XLV+XLMELT
        !STBOLT=5.670151E-8

!--- SNOW flag -- ISICE
!         ILAND=isice

!--- DELTSN - is the threshold for splitting the snow layer into 2 layers.
!--- With snow density 400 kg/m^3, this threshold is equal to 7.5 cm,
!--- equivalent to 0.03 m SNWE. For other snow densities the threshold is
!--- computed using SNWE=0.03 m and current snow density.
!--- SNTH - the threshold below which the snow layer is combined with
!--- the top sea ice layer. SNTH is computed using snwe=0.016 m, and
!--- equals 4 cm for snow density 400 kg/m^3.

! increase thickness of top snow layer from 3 cm SWE to 5 cm SWE
!           DELTSN=5.*SNHEI_CRIT
!           snth=0.4*SNHEI_CRIT

        DELTSN=0.05*1.e3/rhosn
        snth=0.01*1.e3/rhosn
        !snth=0.01601*1.e3/rhosn 

        ! For 2-layer snow model when the snow depth is marginlly higher than DELTSN,
        ! reset DELTSN to half of snow depth.
        IF(SNHEI.GE.DELTSN+SNTH) THEN
          if(snhei-deltsn-snth.lt.snth) deltsn=0.5*(snhei-snth)
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        print *,'DELTSN ICE is changed,deltsn,snhei,snth', &
                                  i,j, deltsn,snhei,snth
    ENDIF
        ENDIF

        RHOICE=900.
        CI=RHOICE*2100.
        RAS=RHO*1.E-3
        RIW=rhoice*1.e-3
        RSM=0.

        XLMELT=3.35E+5
        RHOCSN=2090.* RHOSN
        !18apr08 - add rhonewcsn
        RHOnewCSN=2090.* RHOnewSN
        THDIFSN = 0.265/RHOCSN
        RAS=RHO*1.E-3

        SOILTFRAC=SOILT

        SMELT=0.
        SOH=0.
        SNODIF=0.
        SNOH=0.
        SNOHGNEW=0.
        RSM = 0.
        RSMFRAC = 0.
        fsn=1.
        fso=0.
        cvw=cw

        NZS1=NZS-1
        NZS2=NZS-2

        QGOLD=QSG
        TNOLD=SOILT
        DZSTOP=1./(ZSMAIN(2)-ZSMAIN(1))

        snweprint=0.
        snheiprint=0.
        prcpl=prcpms

        !*** DELTSN is the depth of the top layer of snow where
        !*** there is a temperature gradient, the rest of the snow layer
        !*** is considered to have constant temperature

        H=1.
        SMELT=0.

        FQ=QKMS
        SNHEI=SNWE*1.e3/RHOSN
        SNWEPR=SNWE

        ! check if all snow can evaporate during DT
        BETA=1.
        EPOT = -FQ*(QVATM-QSG)
        EPDT = EPOT * RAS *DELT
        IF(EPDT.GT.0. .and. SNWEPR.LE.EPDT) THEN
           BETA=SNWEPR/max(1.e-8,EPDT)
           SNWE=0.
        ENDIF

!******************************************************************************
!       COEFFICIENTS FOR THOMAS ALGORITHM FOR TSO
!******************************************************************************

        cotso(1)=0.
        rhtso(1)=TSO(NZS)
        DO K=1,NZS2
          KN=NZS-K
          K1=2*KN-3
          X1=DTDZS(K1)*THDIFICE(KN-1)
          X2=DTDZS(K1+1)*THDIFICE(KN)
          FT=TSO(KN)+X1*(TSO(KN-1)-TSO(KN))                           &
             -X2*(TSO(KN)-TSO(KN+1))
          DENOM=1.+X1+X2-X2*cotso(K)
          cotso(K+1)=X1/DENOM
          rhtso(K+1)=(FT+X2*rhtso(K))/DENOM
       ENDDO
       !--- THE NZS element in COTSO and RHTSO will be for snow
       !--- There will be 2 layers in snow if it is deeper than DELTSN+SNTH
       IF(SNHEI.GE.SNTH) then
         if(snhei.le.DELTSN+SNTH) then
         !-- 1-layer snow model
           ilnb=1
           snprim=max(snth,snhei)
           soilt1=tso(1)
           tsob=tso(1)
           XSN = DELT/2./(zshalf(2)+0.5*SNPRIM)
           DDZSN = XSN / SNPRIM
           X1SN = DDZSN * thdifsn
           X2 = DTDZS(1)*THDIFICE(1)
           FT = TSO(1)+X1SN*(SOILT-TSO(1))                              &
                -X2*(TSO(1)-TSO(2))
           DENOM = 1. + X1SN + X2 -X2*cotso(NZS1)
           cotso(NZS)=X1SN/DENOM
           rhtso(NZS)=(FT+X2*rhtso(NZS1))/DENOM
           cotsn=cotso(NZS)
           rhtsn=rhtso(NZS)
           !*** Average temperature of snow pack (C)
           tsnav=0.5*(soilt+tso(1))                                     &
                     -273.15

         else
         !-- 2 layers in snow, SOILT1 is temperasture at DELTSN depth
           ilnb=2
           snprim=deltsn
           tsob=soilt1
           XSN = DELT/2./(0.5*SNHEI)
           XSN1= DELT/2./(zshalf(2)+0.5*(SNHEI-DELTSN))
           DDZSN = XSN / DELTSN
           DDZSN1 = XSN1 / (SNHEI-DELTSN)
           X1SN = DDZSN * thdifsn
           X1SN1 = DDZSN1 * thdifsn
           X2 = DTDZS(1)*THDIFICE(1)
           FT = TSO(1)+X1SN1*(SOILT1-TSO(1))                            &
                -X2*(TSO(1)-TSO(2))
           DENOM = 1. + X1SN1 + X2 - X2*cotso(NZS1)
           cotso(nzs)=x1sn1/denom
           rhtso(nzs)=(ft+x2*rhtso(nzs1))/denom
           ftsnow = soilt1+x1sn*(soilt-soilt1)                          &
                   -x1sn1*(soilt1-tso(1))
           denomsn = 1. + X1SN + X1SN1 - X1SN1*cotso(NZS)
           cotsn=x1sn/denomsn
           rhtsn=(ftsnow+X1SN1*rhtso(NZS))/denomsn
           !*** Average temperature of snow pack (C)
           tsnav=0.5/snhei*((soilt+soilt1)*deltsn                       &
                       +(soilt1+tso(1))*(SNHEI-DELTSN))                 &
                       -273.15
         endif
       ENDIF

       IF(SNHEI.LT.SNTH.AND.SNHEI.GT.0.) then
       !--- snow is too thin to be treated separately, therefore it
       !--- is combined with the first sea ice layer.
         snprim=SNHEI+zsmain(2)
         fsn=SNHEI/snprim
         fso=1.-fsn
         soilt1=tso(1)
         tsob=tso(2)
         XSN = DELT/2./((zshalf(3)-zsmain(2))+0.5*snprim)
         DDZSN = XSN /snprim
         X1SN = DDZSN * (fsn*thdifsn+fso*thdifice(1))
         X2=DTDZS(2)*THDIFICE(2)
         FT=TSO(2)+X1SN*(SOILT-TSO(2))-                              &
                       X2*(TSO(2)-TSO(3))
         denom = 1. + x1sn + x2 - x2*cotso(nzs-2)
         cotso(nzs1) = x1sn/denom
         rhtso(nzs1)=(FT+X2*rhtso(NZS-2))/denom
         tsnav=0.5*(soilt+tso(1))                                    &
                     -273.15
         cotso(nzs)=cotso(NZS1)
         rhtso(nzs)=rhtso(nzs1)
         cotsn=cotso(NZS)
         rhtsn=rhtso(NZS)
       ENDIF

!************************************************************************
!--- THE HEAT BALANCE EQUATION 
        !18apr08 nmelt is the flag for melting, and SNOH 
        ! is heat of snow phase changes
        nmelt=0
        SNOH=0.

        EPOT=-QKMS*(QVATM-QSG)
        RHCS=CAPICE(1)
        H=1.
        FKT=TKMS
        D1=cotso(NZS1)
        D2=rhtso(NZS1)
        TN=SOILT
        D9=THDIFICE(1)*RHCS*dzstop
        D10=TKMS*CP*RHO
        R211=.5*CONFLX/DELT
        R21=R211*CP*RHO
        R22=.5/(THDIFICE(1)*DELT*dzstop**2)
        R6=EMISS *STBOLT*.5*TN**4
        R7=R6/TN
        D11=RNET+R6

      IF(SNHEI.GE.SNTH) THEN 
        if(snhei.le.DELTSN+SNTH) then
        !--- 1-layer snow
          D1SN = cotso(NZS)
          D2SN = rhtso(NZS)
        else
        !--- 2-layer snow
          D1SN = cotsn
          D2SN = rhtsn
        endif
        D9SN= THDIFSN*RHOCSN / SNPRIM
        R22SN = SNPRIM*SNPRIM*0.5/(THDIFSN*DELT)
      ENDIF

      IF(SNHEI.LT.SNTH.AND.SNHEI.GT.0.) then
      !--- thin snow is combined with sea ice
        D1SN = D1
        D2SN = D2
        D9SN = (fsn*THDIFSN*RHOCSN+fso*THDIFICE(1)*RHCS)/        &
               snprim
        R22SN = snprim*snprim*0.5                                &
                /((fsn*THDIFSN+fso*THDIFICE(1))*delt)
      ENDIF
 

      IF(SNHEI.eq.0.)then
      !--- all snow is sublimated
        D9SN = D9
        R22SN = R22
        D1SN = D1
        D2SN = D2
      ENDIF


      !---- TDENOM for snow
      TDENOM = D9SN*(1.-D1SN +R22SN)+D10+R21+R7                    &
              +RAINF*CVW*PRCPMS                                    &
              +RHOnewCSN*NEWSNOW/DELT

      FKQ=QKMS*RHO
      R210=R211*RHO
      AA=XLVM*(BETA*FKQ+R210)/TDENOM
      BB=(D10*TABS+R21*TN+XLVM*(QVATM*                             &
         (BETA*FKQ)                                                &
        +R210*QVG)+D11+D9SN*(D2SN+R22SN*TN)                        &
        +RAINF*CVW*PRCPMS*max(273.15,TABS)                         &
        + RHOnewCSN*NEWSNOW/DELT*min(273.15,TABS)                  &
         )/TDENOM
        AA1=AA
        PP=PATM*1.E3
        AA1=AA1/PP
        !18apr08  - the iteration start point
 212    continue
        BB=BB-SNOH/TDENOM
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        print *,'VILKA-SNOW on SEAICE'
        print *,'tn,aa1,bb,pp,fkq,r210',                             &
                 tn,aa1,bb,pp,fkq,r210
        print *,'TABS,QVATM,TN,QVG=',TABS,QVATM,TN,QVG
    ENDIF

        CALL VILKA(TN,AA1,BB,PP,QS1,TS1,TBQ,KTAU,i,j)
        !--- it is saturation over snow
        QVG=QS1
        QSG=QS1
        QCG=0.

        !--- SOILT - skin temperature
        SOILT=TS1

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        print *,' AFTER VILKA-SNOW on SEAICE'
        print *,' TS1,QS1: ', ts1,qs1
    ENDIF
       ! Solution for temperature at 7.5 cm depth and snow-seaice interface
       IF(SNHEI.GE.SNTH) THEN
         if(snhei.gt.DELTSN+SNTH) then
         !-- 2-layer snow model
           SOILT1=min(273.15,rhtsn+cotsn*SOILT)
           TSO(1)=min(271.4,(rhtso(NZS)+cotso(NZS)*SOILT1))
           tsob=soilt1
         else
         !-- 1 layer in snow
           TSO(1)=min(271.4,(rhtso(NZS)+cotso(NZS)*SOILT))
           SOILT1=TSO(1)
           tsob=tso(1)
         endif
       ELSEIF  (SNHEI > 0. .and. SNHEI < SNTH) THEN
       ! blended
         TSO(2)=min(271.4,(rhtso(NZS1)+cotso(NZS1)*SOILT))
         tso(1)=min(271.4,(tso(2)+(soilt-tso(2))*fso))
         SOILT1=TSO(1)
         tsob=TSO(2)
       ELSE
       ! snow is melted
         TSO(1)=min(271.4,SOILT)
         SOILT1=min(271.4,SOILT)
         tsob=tso(1)
       ENDIF
       !---- Final solution for TSO in sea ice
       IF (SNHEI > 0. .and. SNHEI < SNTH) THEN
       ! blended or snow is melted
         DO K=3,NZS
           KK=NZS-K+1
           TSO(K)=min(271.4,rhtso(KK)+cotso(KK)*TSO(K-1))
         END DO
       ELSE
         DO K=2,NZS
           KK=NZS-K+1
           TSO(K)=min(271.4,rhtso(KK)+cotso(KK)*TSO(K-1))
         END DO
       ENDIF
       !--- For thin snow layer combined with the top soil layer
       !--- TSO(i,j,1) is computed by linear interpolation between SOILT
       !--- and TSO(i,j,2)
       !if(SNHEI.LT.SNTH.AND.SNHEI.GT.0.)then
       !  tso(1)=min(271.4,tso(2)+(soilt-tso(2))*fso)
       !  soilt1=tso(1)
       !  tsob = tso(2)
       !endif

       if(nmelt.eq.1) go to 220

       !--- IF SOILT > 273.15 F then melting of snow can happen
       !   IF(SOILT.GT.273.15.AND.SNWE.GT.0.) THEN
       ! if all snow can evaporate, then there is nothing to melt
       IF(SOILT.GT.273.15.AND.SNWEPR-BETA*EPOT*RAS*DELT.GT.0..AND.SNHEI.GT.0.) THEN
       !
         nmelt = 1
         !soiltfrac=273.15
         soiltfrac=snowfrac*273.15+(1.-snowfrac)*min(271.4,SOILT)

         QSG= QSN(soiltfrac,TBQ)/PP
         T3      = STBOLT*TNold*TNold*TNold
         UPFLUX  = T3 * 0.5*(TNold+SOILTfrac)
         XINET   = EMISS*(GLW-UPFLUX)
         EPOT = -QKMS*(QVATM-QSG)
         Q1=EPOT*RAS

         IF (Q1.LE.0.) THEN
         ! ---  condensation
           DEW=-EPOT

           QFX= XLVM*RHO*DEW
           EETA=QFX/XLVM
         ELSE
         ! ---  evaporation
           EETA = Q1 * BETA *1.E3
           ! to convert from kg m-2 s-1 to m s-1: 1/rho water=1.e-3************
           QFX= - XLVM * EETA
         ENDIF

         HFX=D10*(TABS-soiltfrac)

        IF(SNHEI.GE.SNTH)then
          SOH=thdifsn*RHOCSN*(soiltfrac-TSOB)/SNPRIM
          SNFLX=SOH
        ELSE
          SOH=(fsn*thdifsn*rhocsn+fso*thdifice(1)*rhcs)*                &
               (soiltfrac-TSOB)/snprim
          SNFLX=SOH
        ENDIF
          X= (R21+D9SN*R22SN)*(soiltfrac-TNOLD) +                        &
              XLVM*R210*(QSG-QGOLD)
          !-- SNOH is energy flux of snow phase change
          SNOH=RNET+QFX +HFX                                             &
               +RHOnewCSN*NEWSNOW/DELT*(min(273.15,TABS)-soiltfrac)      &
               -SOH-X+RAINF*CVW*PRCPMS*                                  &
               (max(273.15,TABS)-soiltfrac)

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
     print *,'SNOWSEAICE melt I,J,SNOH,RNET,QFX,HFX,SOH,X',i,j,SNOH,RNET,QFX,HFX,SOH,X
     print *,'RHOnewCSN*NEWSNOW/DELT*(min(273.15,TABS)-soiltfrac)',     &
              RHOnewCSN*NEWSNOW/DELT*(min(273.15,TABS)-soiltfrac)
     print *,'RAINF*CVW*PRCPMS*(max(273.15,TABS)-soiltfrac)',           &
              RAINF*CVW*PRCPMS*(max(273.15,TABS)-soiltfrac)
    ENDIF
          SNOH=AMAX1(0.,SNOH)
          !-- SMELT is speed of melting in M/S
          SMELT= SNOH /XLMELT*1.E-3
          SMELT=AMIN1(SMELT,SNWEPR/DELT-BETA*EPOT*RAS)
          SMELT=AMAX1(0.,SMELT)

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'1-SMELT i,j',smelt,i,j
    ENDIF
          !18apr08 - Egglston limit
          SMELT= amin1 (smelt, 5.6E-8*max(1.,(soilt-273.15)))
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'2-SMELT i,j',smelt,i,j
    ENDIF

          ! rr - potential melting
          rr=SNWEPR/delt-BETA*EPOT*RAS
          SMELT=min(SMELT,rr)
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
      print *,'3- SMELT i,j,smelt,rr',i,j,smelt,rr
    ENDIF
          SNOHGNEW=SMELT*XLMELT*1.E3
          SNODIF=AMAX1(0.,(SNOH-SNOHGNEW))

          SNOH=SNOHGNEW

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print*,'soiltfrac,soilt,SNOHGNEW,SNODIF=', &
            i,j,soiltfrac,soilt,snohgnew,snodif
       print *,'SNOH,SNODIF',SNOH,SNODIF
    ENDIF

          !*** From Koren et al. (1999) 13% of snow melt stays in the snow pack
          rsmfrac=min(0.18,(max(0.08,snwepr/0.10*0.13)))
          if(snhei > 0.01) then
            rsm=rsmfrac*smelt*delt
          else
          ! do not keep melted water if snow depth is less that 1 cm
            rsm=0.
          endif
          !18apr08 rsm is part of melted water that stays in snow as liquid
          SMELT=AMAX1(0.,SMELT-rsm/delt)
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'4-SMELT i,j,smelt,rsm,snwepr,rsmfrac', &
                    i,j,smelt,rsm,snwepr,rsmfrac
    ENDIF

          !-- update liquid equivalent of snow depth
          !-- for evaporation and snow melt
          SNWE = AMAX1(0.,(SNWEPR-                                 &
                 (SMELT+BETA*EPOT*RAS)*DELT                        &
                 !(SMELT+BETA*EPOT*RAS)*DELT*snowfrac               &
                                         ) )
          soilt=soiltfrac
          !--- If there is no snow melting then just evaporation
          !--- or condensation changes SNWE
        ELSE
          if(snhei.ne.0.) then
            EPOT=-QKMS*(QVATM-QSG)
            SNWE = AMAX1(0.,(SNWEPR-                               &
                   BETA*EPOT*RAS*DELT))
                   !BETA*EPOT*RAS*DELT*snowfrac))
          endif

        ENDIF

        ! no iteration for snow on sea ice, because it will produce
        ! skin temperature higher than it is possible with snow on sea ice
        !      if(nmelt.eq.1) goto 212  ! second iteration
 220    continue

        if(smelt > 0..and.  rsm > 0.) then
          if(snwe.le.rsm) then
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
     print *,'SEAICE SNWE<RSM snwe,rsm,smelt*delt,epot*ras*delt,beta', &
                              snwe,rsm,smelt*delt,epot*ras*delt,beta
    ENDIF
          else
          !*** Update snow density on effect of snow melt, melted
          !*** from the top of the snow. 13% of melted water
          !*** remains in the pack and changes its density.
          !*** Eq. 9 (with my correction) in Koren et al. (1999)

            xsn=(rhosn*(snwe-rsm)+1.e3*rsm)/                            &
                 snwe
            rhosn=MIN(MAX(58.8,XSN),500.)

            RHOCSN=2090.* RHOSN
            thdifsn = 0.265/RHOCSN
          endif
        endif

        snweprint=snwe
        snheiprint=snweprint*1.E3 / RHOSN

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
      print *, 'snweprint : ',snweprint
      print *, 'D9SN,SOILT,TSOB : ', D9SN,SOILT,TSOB
    ENDIF

        IF(SNHEI.GT.0.) THEN
          if(ilnb.gt.1) then
            tsnav=0.5/snhei*((soilt+soilt1)*deltsn                     &
                      +(soilt1+tso(1))*(SNHEI-DELTSN))                 &
                      -273.15
          else
            tsnav=0.5*(soilt+tso(1)) - 273.15
          endif
        ENDIF
        !--- RECALCULATION OF DEW USING NEW VALUE OF QSG
        DEW=0.
        PP=PATM*1.E3
        QSG= QSN(SOILT,TBQ)/PP
        EPOT = -FQ*(QVATM-QSG)
        IF(EPOT.LT.0.) THEN
        ! Sublimation
          DEW=-EPOT
        ENDIF

        SNOM=SNOM+SMELT*DELT*1.e3

        !--- THE DIAGNOSTICS OF SURFACE FLUXES

        T3      = STBOLT*TNold*TNold*TNold
        UPFLUX  = T3 *0.5*(SOILT+TNold)
        XINET   = EMISS*(GLW-UPFLUX)
        !RNET    = GSW + XINET
        HFT=-TKMS*CP*RHO*(TABS-SOILT)
        HFX=-TKMS*CP*RHO*(TABS-SOILT)                        &
            *(P1000mb*0.00001/Patm)**ROVCP
        Q1 = - FQ*RAS* (QVATM - QSG)
        IF (Q1.LT.0.) THEN
        ! ---  condensation
          if(myj) then
          !-- moisture flux for coupling with MYJ PBL
            EETA=-QKMS*RAS*(QVATM/(1.+QVATM) - QSG/(1.+QSG))*1.E3
          else ! myj
          !-- actual moisture flux from RUC LSM
            DEW=QKMS*(QVATM-QSG)
            EETA= - RHO*DEW
          endif ! myj
          QFX= XLVm*EETA
          EETA= - RHO*DEW
          sublim = EETA
        ELSE
        ! ---  evaporation
          if(myj) then
          !-- moisture flux for coupling with MYJ PBL
            EETA=-QKMS*RAS*BETA*(QVATM/(1.+QVATM) - QVG/(1.+QVG))*1.E3
          else ! myj
          ! to convert from m s-1 to kg m-2 s-1: *rho water=1.e3************
          !-- actual moisture flux from RUC LSM
            EETA = Q1*BETA*1.E3
          endif ! myj
          QFX= XLVm * EETA
          EETA = Q1*BETA*1.E3
          sublim = EETA
        ENDIF

        icemelt=0.
        IF(SNHEI.GE.SNTH)then
          S=thdifsn*RHOCSN*(soilt-TSOB)/SNPRIM
          SNFLX=S
        ELSEIF(SNHEI.lt.SNTH.and.SNHEI.GT.0.) then
          S=(fsn*thdifsn*rhocsn+fso*thdifice(1)*rhcs)*                &
            (soilt-TSOB)/snprim
          SNFLX=S
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
      print *,'SNOW is thin, snflx',i,j,snflx
    ENDIF
        ELSE 
          SNFLX=D9SN*(SOILT-TSOB)
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
      print *,'SNOW is GONE, snflx',i,j,snflx
    ENDIF
        ENDIF

        SNHEI=SNWE *1.E3 / RHOSN

    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'SNHEI,SNOH',i,j,SNHEI,SNOH
    ENDIF
        X= (R21+D9SN*R22SN)*(soilt-TNOLD) +              &
           XLVM*R210*(QSG-QGOLD)
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
     print *,'SNOWSEAICE storage ',i,j,x
     print *,'R21,D9sn,r22sn,soiltfrac,tnold,qsg,qgold,snprim', &
              R21,D9sn,r22sn,soiltfrac,tnold,qsg,qgold,snprim
    ENDIF
        X=X &
         -RHOnewCSN*NEWSNOW/DELT*(min(273.15,TABS)-SOILT)        &
         -RAINF*CVW*PRCPMS*(max(273.15,TABS)-SOILT)

        ! -- excess energy is spent on ice melt
        icemelt = RNET-HFT-XLVm*EETA-S-SNOH-X
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
        print *,'SNOWSEAICE icemelt=',icemelt
    ENDIF

        FLTOT=RNET-HFT-XLVm*EETA-S-SNOH-x-icemelt
    IF ( wrf_at_debug_level(LSMRUC_DBG_LVL) ) THEN
       print *,'i,j,snhei,qsg,soilt,soilt1,tso,TABS,QVATM', &
                i,j,snhei,qsg,soilt,soilt1,tso,tabs,qvatm
       print *,'SNOWSEAICE - FLTOT,RNET,HFT,QFX,S,SNOH,icemelt,snodif,X,SOILT=' &
                      ,FLTOT,RNET,HFT,XLVm*EETA,s,SNOH,icemelt,snodif,X,SOILT
    ENDIF
        !-- Restore sea-ice parameters if snow is less than threshold
        IF(SNHEI.EQ.0.)  then
          tsnav=soilt-273.15
          emiss=0.98
          znt=0.011
          alb=0.65
        ENDIF

!------------------------------------------------------------------------
   END SUBROUTINE SNOWSEAICE
!------------------------------------------------------------------------

       SUBROUTINE VILKA(TN,D1,D2,PP,QS,TS,TT,NSTEP,ii,j)
!--------------------------------------------------------------
!--- VILKA finds the solution of energy budget at the surface
!--- using table T,QS computed from Clausius-Klapeiron
!--------------------------------------------------------------
   REAL,     DIMENSION(1:5001),  INTENT(IN   )   ::  TT
   REAL,     INTENT(IN  )   ::  TN,D1,D2,PP
   INTEGER,  INTENT(IN  )   ::  NSTEP,ii,j

   REAL,     INTENT(OUT  )  ::  QS, TS

   REAL    ::  F1,T1,T2,RN
   INTEGER ::  I,I1
     
       I=(TN-1.7315E2)/.05+1
       T1=173.1+FLOAT(I)*.05
       F1=T1+D1*TT(I)-D2
       I1=I-F1/(.05+D1*(TT(I+1)-TT(I)))
       I=I1
       IF(I.GT.5000.OR.I.LT.1) GOTO 1
  10   I1=I
       T1=173.1+FLOAT(I)*.05
       F1=T1+D1*TT(I)-D2
       RN=F1/(.05+D1*(TT(I+1)-TT(I)))
       I=I-INT(RN)                      
       IF(I.GT.5000.OR.I.LT.1) GOTO 1
       IF(I1.NE.I) GOTO 10
       TS=T1-.05*RN
       QS=(TT(I)+(TT(I)-TT(I+1))*RN)/PP
       GOTO 20
   1   PRINT *,' SEA-ICE, AVOST IN VILKA   Table index= ',I
       print *,'I,J=',ii,j,'Psfc[hPa] = ',pp, 'Tsfc = ',tn
       CALL wrf_error_fatal ('  Crash in surface energy budget  ' )
   20  CONTINUE
!-----------------------------------------------------------------------
   END SUBROUTINE VILKA
!-----------------------------------------------------------------------

       FUNCTION QSN(TN,T)
!-----------------------------------------------------------------------
   REAL,     DIMENSION(1:5001),  INTENT(IN   )   ::  T
   REAL,     INTENT(IN  )   ::  TN

      REAL    QSN, R,R1,R2
      INTEGER I

       R=(TN-173.15)/.05+1.
       I=INT(R)
       IF(I.GE.1) goto 10
       I=1
       R=1.
  10   IF(I.LE.5000) GOTO 20
       I=5000
       R=5001.
  20   R1=T(I)
       R2=R-I
       QSN=(T(I+1)-R1)*R2 + R1
!       print *,' in QSN, I,R,R1,R2,T(I+1),TN, QSN', I,R,r1,r2,t(i+1),tn,QSN
!-----------------------------------------------------------------------
  END FUNCTION QSN
!-----------------------------------------------------------------

END MODULE module_sf_ruclsm_ice
