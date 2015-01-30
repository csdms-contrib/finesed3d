c  FLTIMER = whether to do timing
c  FLAVER  = whether to do axial average
c  FLTHRM  = whether to include temperature
c  FLPART  = whether to compute particles velocity
c  FLTWC   = whether to include 2-way coupling
c  FLOTAU  = whether settling is order tau (otherwise O(1))
c  FLNBSQ  = whether non-Boussinesq term are computed
c  FLMDRHS = whether mat der is computed using the RHS
c  FLGVAR  = whether gravity vector changes with y or z
c  FLHINDER = whether to include hindered settling
c  FLEro = wether to include erosional boundary

      LOGICAL FLTIMER,FLAVER,FLTHRM,FLPART,FLSETTLE,
     *        FLINERTIA,FLTWC,FLOTAU,FLHINDER,
     *        FLNBSQ,FLMDRHS,FLSV_YES,FLSV_NO,FLGVAR,
     *        FLEro
      COMMON /flags/ FLTIMER,FLAVER,FLTHRM,FLPART,
     *               FLSETTLE,FLINERTIA,FLTWC,FLHINDER,
     *               FLOTAU,FLNBSQ,FLMDRHS,FLSV_YES,FLSV_NO,
     *               FLGVAR,FLEro  

