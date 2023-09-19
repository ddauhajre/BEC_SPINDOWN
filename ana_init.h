
      ! Everything after the implicit none

      ! This initial condition assumes linera eof
      ! with Tcoef=2.0e-4

      integer :: istr,iend,jstr,jend, i,j,k, itrc, ierr
# ifdef BIOLOGY_BEC2
      real :: itrc_ind, ibio_strt

      !Parameters for vertical shape function of idealized BEC I.C.
      real :: bgc_prof
      real :: bec_alpha = 0.1
# endif

 
      real :: Sconst          = 32       ! Const. salinity
      real :: h_sbl       = 60       ! Surface BL depth


      ! necessary for FIRST_TIME_STEP flag to work, as per set_global_definitions.h
      ! this had no value for analytical examples previously.
#ifdef EXACT_RESTART
      forw_start=ntstart
#endif

      
      !---------------------- Iniitalize temp,salinity ---------------------------  
      
      ! Calculate temperature from buoyancy:
      ! Make sure that you are using linear equation of state

      do j=-1,ny+2
        do i=-1,nx+2
          do k=1,nz
	     ! 		FILL IN HERE WITH OTHER RECIPES FOR T,S
             t(i,j,k,1,itemp) = 15 !set to constant for now
          enddo
        enddo
      enddo

      do k=1,nz
         do j=-1,ny+2
	    do i=-1,nx+2
	       t(i,j,k,2,itemp) = t(i,j,k,1,itemp)
# ifdef SALINITY
            t(i,j,k,1,isalt)=Sconst
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
# endif
	    enddo
	 enddo
      enddo
      !-------------------------------------------------------------------
 
      !-------------------- Initialize no flow -----------------------
      do j=-1,ny+2
         do i=-1,nx+2
	      zeta(i,j,1)=0
	      zeta(i,j,2)=zeta(i,j,1)
	      ubar(i,j,1)=0
	      ubar(i,j,2)=ubar(i,j,1)
	      vbar(i,j,1)=0
	      vbar(i,j,2)=vbar(i,j,1)
	 enddo
      enddo
      do k=1,nz 
         do j=-1,ny+2
            do i=-1,nx+2
              u(i,j,k,1)=0
              u(i,j,k,2)=u(i,j,k,1)
              v(i,j,k,1)=0
              v(i,j,k,2)=v(i,j,k,1) 
	    enddo
	 enddo
      enddo  
      !---------------------------------------------------------------

      !-------------------- Initialize BEC tracers -----------------------
#ifdef BIOLOGY_BEC2
#ifdef SALINITY
      ibio_strt=isalt+1
# else
      ibio_strt=itemp+1
#endif
      !print *, 'ibio_str', ibio_strt
      !print *, 'nt', nt
      do itrc_ind=ibio_strt,nt
         do k=1,nz
            do j=-1,ny+2
              do i=-1,nx+2
	        
	        !Construct vertical shape function
		bgc_prof = 0.5 * (1 + t_anagamma(itrc_ind)* tanh(bec_alpha*(z_r(i,j,k) +h_sbl)))
	        
	        !Normalize to min/max of tracer
		t(i,j,k,1,itrc_ind)=t_anamin(itrc_ind) + (t_anamax(itrc_ind) - t_anamin(itrc_ind)) * bgc_prof
		
		t(i,j,k,2,itrc_ind)=t(i,j,k,1,itrc_ind)
	     enddo
           enddo
	 enddo
      enddo
#endif
      !-------------------------------------------------------------------

