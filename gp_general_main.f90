module gp_general_main_mod

use prec_mod
use constants_mod
use verbosity_mod
use lbfgs_mod
use lbfgsb_mod
use positions_mod
use clustering_mod
use gp_basic_dat_mod
use gp_basic_mod
use gp_general_dat_mod

implicit none

contains

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_core

  implicit none

  call gp_general_main_find_pos_min_max
  call gp_general_main_pred_pos
  call gp_general_main_sparsification
  call gp_general_main_calc
  call gp_general_main_find_function_min
  call gp_general_main_find_path

end subroutine gp_general_main_core

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_find_pos_min_max

  implicit none

  integer :: i

  if (.not. (inval .or. inder) ) return

  if (.not. genrange) then
     if(verbose) call verbosity_start_task('Finding min and max positions for general range')
     call positions_find_min_max(positions, abs(periodicities), min_positions, max_positions, range_positions)
     if(verbose) call verbosity_finish_task('Finding min and max positions for general range')
  else
     do i=1,ndim
        range_positions(i) = max_positions(i) - min_positions(i)
        if (range_positions(i) < 0) range_positions(i) = range_positions(i) + abs(periodicities(i))
     end do
  end if

end subroutine gp_general_main_find_pos_min_max

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_pred_pos

  implicit none

  integer :: i, j, k, n
  integer :: ios
  integer :: row
  integer :: alloc_failed
  logical :: file_exist
  integer :: multigrids

  if (.not. (outval .or. outvaldev .or. outder .or. outderdev .or. outhess .or. outint .or. outintdev .or. outmin) ) return

  ! prediction points
  ! from file
  if(.not. predgriddata) then

    if(verbose) call verbosity_start_task('Reading prediction positions from file')

    inquire(file=predposfile, exist=file_exist)
    if(.not. file_exist) then
        write(6,*) '>>> ERROR: PREDPOSFILE("', trim(adjustl(predposfile)), '") does not exist!'
        stop
    end if

    open(10,file=predposfile)

    ios=0
    n=0
    do while (ios .eq. 0)
       read(10,*,iostat=ios)
       n=n+1
    end do
    n=n-1

    rewind(10)

    pred_tot_grids = n
    allocate(pred_grid_positions(ndim,pred_tot_grids), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        write(6,*) '>>> ERROR: [gp_general_main_grid] Unable to allocate memory for pred_grid_positions!'
        stop
    end if

    do i=1,pred_tot_grids
     if(predrownum) then
        read(10,*) row, pred_grid_positions(:,i)
     else
        read(10,*) pred_grid_positions(:,i)
     end if
    end do

    close(10)

    if(verbose) call verbosity_finish_task('Reading prediction positions from file')

  ! equidistant in the region
  else

    if (.not. predrange) then
       pred_min_positions = min_positions
       pred_max_positions = max_positions
       pred_range_positions = range_positions
    else
       do i=1,ndim
          pred_range_positions(i) = pred_max_positions(i) - pred_min_positions(i)
          if (pred_range_positions(i) < 0) pred_range_positions(i) = pred_range_positions(i) + abs(periodicities(i))
       end do
    end if

    if(verbose) call verbosity_start_task('Constructing equidistant grid for prediction')

    do i=1,ndim
       pred_widths(i) = pred_range_positions(i) / real( pred_grid(i) )
    end do

    ! calculate total number of bins
    pred_tot_grids=1
    do i=1,ndim
       pred_tot_grids = pred_tot_grids * pred_grid(i)
    end do

    ! calculate grid positions
    allocate(pred_grid_positions(ndim,pred_tot_grids), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        write(6,*) '>>> ERROR: [gp_general_main_grid] Unable to allocate memory for gridpositions!'
        stop
    end if

    do n=1,pred_tot_grids
       do i=1,ndim
          multigrids = 1
          do j=i+1,ndim
             multigrids = multigrids * pred_grid(j)
          end do

          k = mod(((n-1)/multigrids),pred_grid(i))

          pred_grid_positions(i,n) = pred_min_positions(i) + real(k)*pred_widths(i) + pred_widths(i) / 2.0_dp
          if(periodicities(i) .ne. 0) pred_grid_positions(i,n) = pred_grid_positions(i,n) &
                                       - nint(pred_grid_positions(i,n)/abs(periodicities(i)))*abs(periodicities(i))
       end do
    end do

    if(verbose) call verbosity_finish_task('Constructing equidistant grid for prediction')

  end if

10 format(100A)
20 format(I2,1X,I5,1X,F12.5,1X,F12.5,1X,F12.5,1X,F12.5)

end subroutine gp_general_main_pred_pos

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_sparsification

  implicit none

  select case(sparsification)
     case('none')
     case('clustering')
       call gp_general_main_sparsification_clustering
     case('grid')
       call gp_general_main_sparsification_grid
     case('file')
       call gp_general_main_sparsification_file
     case('gpfile')
       call gp_general_main_sparsification_gpfile
  end select

end subroutine gp_general_main_sparsification

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_sparsification_clustering

  implicit none

  integer :: i
  integer :: alloc_failed

  if ( ntot .le. sparsepoints ) then
       if(verbose) then
          write(6,*) '>>>INFO: ntot(',ntot,') is less or equal than sparsepoints(',sparsepoints,')'
          write(6,*) '>>>INFO: clustering is skipped'
       end if

       clustering = 'every'
       sparsepoints = ntot
       sparsefreq = 1
  end if

  allocate(index_sparse(sparsepoints), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_sparsification_clustering] Unable to allocate memory for index_sparse!'
      stop
  end if

  if ( sparsepoints .gt. (ntot / sparsefreq) ) then
       if(verbose) then
          write(6,*) '>>>INFO: ntot/sparsefreq(',ntot/sparsefreq,') is less than sparsepoints(',sparsepoints,')'
          write(6,*) '>>>INFO: sparsefreq is adjusted accordingly (= ntot/sparsepoints)'
       end if
       sparsefreq = ntot / sparsepoints
  end if

  select case(clustering)
     case('every')
       if(verbose) call verbosity_start_task('Every clustering')
       index_sparse = (/ (i, i=1, sparsepoints) /)
       if(verbose) call verbosity_finish_task('Every clustering')
     case('random')
       if(verbose) call verbosity_start_task('Random clustering')
       call random_clustering_pick(positions(:,sparsefreq:ntot:sparsefreq),index_sparse)
       if(verbose) call verbosity_finish_task('Random clustering')
     case('kmeans')
       if(verbose) call verbosity_start_task('Kmeans clustering')
       call k_means_clustering_pick(positions(:,sparsefreq:ntot:sparsefreq),abs(periodicities),index_sparse, &
                                              thresholds=kmeansthresholds,verbose=verbose)
       if(verbose) call verbosity_finish_task('Kmeans clustering')
     case('kmeans++')
       if(verbose) call verbosity_start_task('Kmeans++ clustering')
       call k_meanspp_clustering_pick(positions(:,sparsefreq:ntot:sparsefreq),abs(periodicities),index_sparse)
       if(verbose) call verbosity_finish_task('Kmeans++ clustering')

       if(verbose) call verbosity_start_task('Kmeans clustering')
       call k_means_clustering_pick(positions(:,sparsefreq:ntot:sparsefreq),abs(periodicities),index_sparse, &
                                    thresholds=kmeansthresholds,verbose=verbose)
       if(verbose) call verbosity_finish_task('Kmeans clustering')
  end select

  index_sparse = sparsefreq * index_sparse

  allocate(pseudo_sparse(ndim,sparsepoints))
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_sparsification_clustering] Unable to allocate memory for pseudo_sparse!'
      stop
  end if

  pseudo_sparse = positions(:,index_sparse)

end subroutine gp_general_main_sparsification_clustering

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_sparsification_grid

  implicit none

  integer :: i, j, k, n 
  integer :: alloc_failed
  integer :: multigrids

  if(verbose) call verbosity_start_task('Sparse points on grid')

  if (.not. predrange) then
     sparse_min_positions = min_positions
     sparse_max_positions = max_positions
     sparse_range_positions = range_positions
  else
     do i=1,ndim
        sparse_range_positions(i) = sparse_max_positions(i) - sparse_min_positions(i)
        if (sparse_range_positions(i) < 0) sparse_range_positions(i) = sparse_range_positions(i) + abs(periodicities(i))
     end do
  end if

  do i=1,ndim
     sparse_widths(i) = sparse_range_positions(i) / real( sparse_grid(i) )
  end do

  sparsepoints = 1
  do i=1,ndim
     sparsepoints = sparsepoints * sparse_grid(i)
  end do

  allocate(pseudo_sparse(ndim,sparsepoints), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then             
      write(6,*) '>>> ERROR: [gp_general_main_sparsification_grid] Unable to allocate memory for pseudo_sparse!'
      stop                                   
  end if

  do n=1,sparsepoints
     do i=1,ndim
        multigrids = 1
        do j=i+1,ndim
           multigrids = multigrids * sparse_grid(j)
        end do

        k = mod(((n-1)/multigrids),sparse_grid(i))

        pseudo_sparse(i,n) = sparse_min_positions(i) + real(k)*sparse_widths(i) + sparse_widths(i) / 2.0_dp
        if ( periodicities(i) .ne. 0 ) pseudo_sparse(i,n) = pseudo_sparse(i,n) &
                                          - nint(pseudo_sparse(i,n)/abs(periodicities(i)))*abs(periodicities(i))
     end do
  end do

  if(verbose) call verbosity_finish_task('Sparse points on grid')

end subroutine gp_general_main_sparsification_grid

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_sparsification_file

  implicit none

  integer :: i, n
  integer :: ios
  integer :: row
  integer :: alloc_failed
  logical :: file_exist

  if(verbose) call verbosity_start_task('Sparse points from file')

  inquire(file=sparseposfile, exist=file_exist)
  if(.not. file_exist) then
      write(6,*) '>>> ERROR: SPARSEPOSFILE("', trim(adjustl(sparseposfile)), '") does not exist!'
      stop
  end if

  open(11,file=sparseposfile)

  ios=0
  n=0
  do while (ios .eq. 0)
     read(11,*,iostat=ios)
     n=n+1
  end do
  n=n-1

  rewind(11)

  sparsepoints = n
  allocate(pseudo_sparse(ndim,sparsepoints), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_sparsification_file] Unable to allocate memory for pseudo_sparse!'
      stop
  end if

  do i=1,sparsepoints
   if(sparserownum) then
      read(11,*) row, pseudo_sparse(:,i)
   else
      read(11,*) pseudo_sparse(:,i)
   end if
  end do

  close(4)

  if(verbose) call verbosity_finish_task('Sparse points from file')

end subroutine gp_general_main_sparsification_file

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_sparsification_gpfile

  implicit none

  integer :: alloc_failed

  if(gp%m_f .ne. 0) then
    allocate(pseudo_sparse(gp%n_dof,gp%m_f), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        write(6,*) '>>> ERROR: [gp_general_main_sparsification_gpfile] Unable to allocate memory for pseudo_sparse!'
        stop
    end if
    pseudo_sparse = gp%f_r
  else
    allocate(pseudo_sparse(gp%n_dof,gp%m_g), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        write(6,*) '>>> ERROR: [gp_general_main_sparsification_gpfile] Unable to allocate memory for pseudo_sparse!'
        stop
    end if
    pseudo_sparse = gp%g_r
  end if

end subroutine gp_general_main_sparsification_gpfile

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_calc

  implicit none

  if( .not. outopt ) then
      ! change
      call gp_general_main_change
  else
      ! optimization
      call gp_general_main_hyper_opt
  end if

  call gp_general_main_teaching
  call gp_general_main_prediction

end subroutine gp_general_main_calc

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_hyper_opt

  implicit none

  integer :: i, j, n
  real(dp) :: rmsg, gradmax, fval, lastfval, eps
  integer :: compmax
  real(dp), allocatable :: posi(:), grad(:)
  integer :: alloc_failed
  integer :: found
  real(dp), allocatable :: lower_bound(:), upper_bound(:)
  ! lbfgs
  integer, allocatable :: bound(:)
  real(dp), allocatable :: wa(:)
  integer, allocatable :: iwa(:)
  character(len=60) :: task, csave
  logical, allocatable :: lsave(:)
  integer, allocatable :: isave(:)
  real(dp), allocatable :: dsave(:)
  integer :: iprint
  

  ! init variables
  lastfval = 0.0_dp
  eps = epsilon(0.0_dp)

  ! check how many variables are to optimize actually
  if( .not. inval ) optsigmaval = 0.0_dp
  if( .not. inder ) optsigmader(:) = 0.0_dp

  optdim = 0
  if( optdelta .ne. 0.0_dp ) then
      optdim = optdim + 1
      orig_delta = delta
  end if
  if( optsigmaval .ne. 0.0_dp ) then
      optdim = optdim + 1
      allocate(orig_sigmaval(ntot), stat = alloc_failed)
      if( alloc_failed .ne. 0 ) then
          write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for orig_sigmaval!'
          stop
      end if
      orig_sigmaval = sigmaval
  end if
  if( any(optsigmader .ne. 0.0_dp) ) then
      optdim = optdim + count(optsigmader .ne. 0.0_dp)
      allocate(orig_sigmader(ndim,ntot), stat = alloc_failed)
      if( alloc_failed .ne. 0 ) then
          write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for orig_sigmader!'
          stop
      end if
      orig_sigmader = sigmader
  end if
  if( any(opttheta .ne. 0.0_dp) ) then
      optdim = optdim + count(opttheta .ne. 0.0_dp)
      allocate(orig_thetas(ndim), stat = alloc_failed)
      if( alloc_failed .ne. 0 ) then
          write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for orig_thetas!'
          stop
      end if
      orig_thetas = thetas
  end if

  if( optdim .eq. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] dimension of variables to be optimized is 0!'
      stop
  end if

  ! initialize arrays
  allocate(optpos(optdim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for optpos!'
      stop
  end if

  allocate(optgrad(optdim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for optgrad!'
      stop
  end if

  allocate(posi(optdim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for posi!'
      stop
  end if

  allocate(grad(optdim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for grad!'
      stop
  end if

  allocate(lower_bound(optdim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for lower_bound!'
      stop
  end if

  allocate(upper_bound(optdim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for upper_bound!'
      stop
  end if

  select case(numdiff)
     case('one-sided')
       n = 0
       if( optdelta .ne. 0.0_dp ) then
           n = n + 1
           lower_bound(n) = optdelta
       end if
       if( optsigmaval .ne. 0.0_dp ) then
           n = n + 1
           lower_bound(n) = optsigmaval
       end if
       do j=1, ndim
          if( optsigmader(j) .ne. 0.0_dp ) then
              n = n + 1
              lower_bound(n) = optsigmader(j)
          end if
       end do
       do j=1, ndim
          if( opttheta(j) .ne. 0.0_dp ) then
              n = n + 1
              lower_bound(n) = opttheta(j)
          end if
       end do
     case('two-sided')
       n = 0
       if( optdelta .ne. 0.0_dp ) then
           n = n + 1
           lower_bound(n) = two*optdelta
       end if
       if( optsigmaval .ne. 0.0_dp ) then
           n = n + 1
           lower_bound(n) = two*optsigmaval
       end if
       do j=1, ndim
          if( optsigmader(j) .ne. 0.0_dp ) then
              n = n + 1
              lower_bound(n) = two*optsigmader(j)
          end if
       end do 
       do j=1, ndim
          if( opttheta(j) .ne. 0.0_dp ) then
              n = n + 1
              lower_bound(n) = two*opttheta(j)
          end if
       end do
  end select

  ! allocate working arrays
  select case(optimizer)
     case('lbfgs', 'lbfgsb')
       allocate(bound(optdim), stat = alloc_failed)
       if( alloc_failed .ne. 0 ) then
           write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for bound!'
           stop
       end if

       allocate(wa(2*ncorr*optdim+5*optdim+11*ncorr*ncorr+8*ncorr), stat = alloc_failed)
       if( alloc_failed .ne. 0 ) then
           write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for wa!'
           stop
       end if

       allocate(iwa(3*optdim), stat = alloc_failed)
       if( alloc_failed .ne. 0 ) then
           write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for iwa!'
           stop
       end if

       allocate(lsave(4), stat = alloc_failed)
       if( alloc_failed .ne. 0 ) then
           write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for lsave!'
           stop
       end if

       allocate(isave(44), stat = alloc_failed)
       if( alloc_failed .ne. 0 ) then
           write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for isave!'
           stop
       end if

       allocate(dsave(29), stat = alloc_failed)
       if( alloc_failed .ne. 0 ) then
           write(6,*) '>>> ERROR: [gp_general_main_hyper_opt] Unable to allocate memory for dsave!'
           stop
       end if

       bound(:) = 1
       task = 'START'
  end select

  if(verbose) call verbosity_start_task('Optimization')

  found = 0
  posi(:) = 1.0_dp

  do i=1,maxiter-1

     ! calculate logL and gradients
     select case(numdiff)
        case('one-sided')
          call gp_general_main_hyper_opt_numdiff_one_sided(posi, fval, grad)
        case('two-sided')
          call gp_general_main_hyper_opt_numdiff_two_sided(posi, fval, grad)
     end select

     ! calculate maxgrad and rmsg
     rmsg = grad_rmsg(grad,gradmax,compmax)
     ! check function value change
     if( (i .ne. 1) .and. (abs(fval-lastfval) .le. maxdval) ) then
         if(verbose) then
            write(6,'(a,f15.8)') '>>> INFO: RMS of gradient: ', rmsg
            write(6,'(a,f15.8)') '>>> INFO: RMS of gradient treshold: ', maxrmsgrad
            write(6,'(a,f15.8)') '>>> INFO: Max gradient component: ', abs(gradmax)
            write(6,'(a,f15.8)') '>>> INFO: Max gradient component treshold: ', maxgrad
            write(6,'(a,f15.8)') '>>> INFO: Last function value change: ', abs(fval-lastfval)
            write(6,'(a,f15.8)') '>>> INFO: Function value change treshold: ', maxdval
            write(6,'(a)') '>>> INFO: Function value change is below trashold! Optimization is stopped.'
         end if
         found=i
         exit
     end if
     ! check gradient
     if( (abs(gradmax) .le. maxgrad) .and. (rmsg .le. maxrmsgrad) ) then
         if(verbose) then
            write(6,'(a,f15.8)') '>>> INFO: RMS of gradient: ', rmsg
            write(6,'(a,f15.8)') '>>> INFO: RMS of gradient treshold: ', maxrmsgrad
            write(6,'(a,f15.8)') '>>> INFO: Max gradient component: ', abs(gradmax)
            write(6,'(a,f15.8)') '>>> INFO: Max gradient component treshold: ', maxgrad
            write(6,'(a,f15.8)') '>>> INFO: Last function value change: ', abs(fval-lastfval)
            write(6,'(a,f15.8)') '>>> INFO: Function value change treshold: ', maxdval
            write(6,'(a)') '>>> INFO: Gradient tresholds were satisfied! Optimization is stopped.'
         end if
         found=i
         exit
     end if
     ! print intermediate results
     if(verbose) then
        write(6,'(a,f20.8)') '>>> INFO: RMS of gradient: ', rmsg
        write(6,'(a,f20.8,a,i8)') '>>> INFO: Max gradient component: ', abs(gradmax), ' component: ', compmax
        write(6,'(a,f20.8)') '>>> INFO: Function value: ', fval
        if(i .ne. 1) write(6,'(a,f15.8)') '>>> INFO: function value change: ', abs(fval-lastfval)
        write(6,'(a)', advance='no') '>>> INFO: Position of components (factor): '
        do j=1, optdim
           write(6,'(f15.8)', advance='no') posi(j)
        end do
        write(6,*)
     end if

     select case(optimizer)
        case('lbfgs', 'lbfgsb')
          if( (task(1:5) /= 'START') .and. (task(1:2) /= 'FG') .and. (task(1:5) /= 'NEW_X') ) exit
          do while(.true.)
             call lbfgsb(optdim,ncorr,posi,lower_bound,upper_bound,bound,-fval,-grad, &
                         0.0_dp,0.0_dp,wa,iwa,task,-1,csave,lsave,isave,dsave)
             if(task(1:5) /= 'NEW_X' .and. task(1:8) /= 'FG_START') exit
          end do
     end select

     lastfval = fval
  end do
  ! store results
  optfound = found
  optpos(:) = posi(:)
  optfval = lastfval
  optgrad(:) = grad(:)

  call gp_general_main_set_hyper_actual(posi)
  call finalise(self=gp)

  if(verbose) call verbosity_finish_task('Optimization')

end subroutine gp_general_main_hyper_opt

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_hyper_opt_numdiff_one_sided(posi, logL, grad)

  implicit none

  real(dp), intent(in) :: posi(:)
  real(dp), intent(out) :: logL
  real(dp), intent(out) :: grad(:)

  integer :: n, i

  ! first: calculate logL at the actual position
  call gp_general_main_set_hyper_actual(posi)

  call finalise(self=gp)
  call gp_general_main_teaching

  logL = gp%logL

  ! second: calculate logL at the corresponding shifted positions
  n = 0
  if( optdelta .ne. 0.0_dp ) then
      n = n + 1
      call gp_general_main_set_hyper_actual(posi)
      delta = orig_delta * (posi(n) + optdelta)

      call finalise(self=gp)
      call gp_general_main_teaching

      grad(n) = (gp%logL - logL) / optdelta
  end if
  if( optsigmaval .ne. 0.0_dp ) then
      n = n + 1
      call gp_general_main_set_hyper_actual(posi)
      sigmaval = orig_sigmaval * (posi(n) + optsigmaval)

      call finalise(self=gp)
      call gp_general_main_teaching

      grad(n) = (gp%logL - logL) / optsigmaval
  end if
  do i=1,ndim
     if( optsigmader(i) .ne. 0.0_dp ) then
         n = n + 1
         call gp_general_main_set_hyper_actual(posi)
         sigmader(:,i) = orig_sigmader(:,i) * (posi(n) + optsigmader(i))

         call finalise(self=gp)
         call gp_general_main_teaching

         grad(n) = (gp%logL - logL) / optsigmader(i)
     end if
  end do
  do i=1,ndim
     if( opttheta(i) .ne. 0.0_dp ) then
         n = n + 1
         call gp_general_main_set_hyper_actual(posi)
         thetas(i) = orig_thetas(i) * (posi(n) + opttheta(i))

         call finalise(self=gp)
         call gp_general_main_teaching

         grad(n) = (gp%logL - logL) / opttheta(i)
     end if
  end do

  call finalise(self=gp)

end subroutine gp_general_main_hyper_opt_numdiff_one_sided

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_hyper_opt_numdiff_two_sided(posi, logL, grad)

  implicit none

  real(dp), intent(in) :: posi(:)
  real(dp), intent(out) :: logL
  real(dp), intent(out) :: grad(:)

  integer :: n, i

  ! first: calculate logL at the corresponding shifted positions
  n = 0
  if( optdelta .ne. 0.0_dp ) then
      n = n + 1
      call gp_general_main_set_hyper_actual(posi)
      delta = orig_delta * (posi(n) + optdelta)

      call finalise(self=gp)
      call gp_general_main_teaching

      logL = gp%logL

      call gp_general_main_set_hyper_actual(posi)
      delta = orig_delta * (posi(n) - optdelta)

      call finalise(self=gp)
      call gp_general_main_teaching

      grad(n) = (logL - gp%logL) / optdelta / 2.0_dp
  end if
  if( optsigmaval .ne. 0.0_dp ) then
      n = n + 1
      call gp_general_main_set_hyper_actual(posi)
      sigmaval = orig_sigmaval * (posi(n) + optsigmaval)

      call finalise(self=gp)
      call gp_general_main_teaching

      logL = gp%logL

      call gp_general_main_set_hyper_actual(posi)
      sigmaval = orig_sigmaval * (posi(n) - optsigmaval)

      call finalise(self=gp)
      call gp_general_main_teaching

      grad(n) = (logL - gp%logL) / optsigmaval / 2.0_dp
  end if
  do i=1,ndim
     if( optsigmader(i) .ne. 0.0_dp ) then
         n = n + 1
         call gp_general_main_set_hyper_actual(posi)
         sigmader(:,i) = orig_sigmader(:,i) * (posi(n) + optsigmader(i))

         call finalise(self=gp)
         call gp_general_main_teaching

         logL = gp%logL

         call gp_general_main_set_hyper_actual(posi)
         sigmader(:,i) = orig_sigmader(:,i) * (posi(n) - optsigmader(i))

         call finalise(self=gp)
         call gp_general_main_teaching

         grad(n) = (logL - gp%logL) / optsigmader(i) / 2.0_dp
     end if
  end do
  do i=1,ndim
     if( opttheta(i) .ne. 0.0_dp ) then
         n = n + 1
         call gp_general_main_set_hyper_actual(posi)
         thetas(i) = orig_thetas(i) * (posi(n) + opttheta(i))

         call finalise(self=gp)
         call gp_general_main_teaching

         logL = gp%logL

         call gp_general_main_set_hyper_actual(posi)
         thetas(i) = orig_thetas(i) * (posi(n) - opttheta(i))

         call finalise(self=gp)
         call gp_general_main_teaching

         grad(n) = (logL - gp%logL) / opttheta(i) / 2.0_dp
     end if
  end do

  ! second: calculate logL at the actual position
  call gp_general_main_set_hyper_actual(posi)

  call finalise(self=gp)
  call gp_general_main_teaching

  logL = gp%logL

end subroutine gp_general_main_hyper_opt_numdiff_two_sided

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_set_hyper_actual(factor)

  implicit none

  real(dp), intent(in) :: factor(:)

  integer :: n, i

  n = 0

  if( optdelta .ne. 0.0_dp ) then
      n = n + 1
      delta = orig_delta * factor(n)
  end if

  if( optsigmaval .ne. 0.0_dp ) then
      n = n + 1
      sigmaval = orig_sigmaval * factor(n)
  end if

  if( any(optsigmader .ne. 0.0_dp) ) then
      do i=1,ndim
         if( optsigmader(i) .ne. 0.0_dp ) then
             n = n + 1
             sigmader(i,:) = orig_sigmader(i,:) * factor(n)
         end if
      end do
  end if

  if( any(opttheta .ne. 0.0_dp) ) then
      do i=1,ndim
         if( opttheta(i) .ne. 0.0_dp ) then
             n = n + 1
             thetas(i) = orig_thetas(i) * factor(n)
         end if
      end do
  end if

end subroutine gp_general_main_set_hyper_actual

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_change

  implicit none 

  if (.not. dochange) return
  select case(sparsetype)
     case('all')
       call gp_basic_change(self=gp, f_set=pseudo_sparse, g_set=pseudo_sparse, verbose=verbose)
     case('val')
       call gp_basic_change(self=gp, f_set=pseudo_sparse, verbose=verbose)
     case('der')
       call gp_basic_change(self=gp, g_set=pseudo_sparse, verbose=verbose)
  end select

end subroutine gp_general_main_change

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_teaching

  implicit none

  integer :: istart, istop

  if( .not. (inval .or. inder) ) return

  if(verbose) call verbosity_start_task('GP covariance calculation')

  ! calculate covariance matrix
  ! conventional GP
  if ( sparsification .eq. 'none' ) then

     ! do GP covariance matrix calculation
     if(inval .and. inder) then
        call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                            f_r=positions,f_v=values,f_s=sigmaval,g_r=positions,g_v=derivatives,g_s=sigmader, &
                            kernel_string=kernel_string,var_names=var_names, &
                            delta=delta,theta=thetas,periodicity=periodicities, &
                            jitter=jitter,factorization=factorization,verbose=verbose)
     else if(inval) then
        call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                            f_r=positions,f_v=values,f_s=sigmaval, &
                            kernel_string=kernel_string,var_names=var_names, &
                            delta=delta,theta=thetas,periodicity=periodicities, &
                            jitter=jitter,factorization=factorization,verbose=verbose)
     else if(inder) then
        call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                            g_r=positions,g_v=derivatives,g_s=sigmader, &
                            kernel_string=kernel_string,var_names=var_names, &
                            delta=delta,theta=thetas,periodicity=periodicities, &
                            jitter=jitter,factorization=factorization,verbose=verbose)
     end if

  ! sparsified GP with multiple iteration if required 
  else

     ! do GP covariance matrix calculation
     istart = 1
     istop = sparsemaxn / itot

     do while(istart .le. ntot)
         if(verbose) write(*,*) "istart = ", istart, " istop = ", istop

         if(inval .and. inder) then
            select case(sparsetype)
               case('all')
                 call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                                     f_r=positions(:,istart:istop),f_v=values(istart:istop),f_s=sigmaval(istart:istop), &
                                     g_r=positions(:,istart:istop),g_v=derivatives(:,istart:istop), &
                                     g_s=sigmader(:,istart:istop), &
                                     f_sparse_r=pseudo_sparse, g_sparse_r=pseudo_sparse, &
                                     kernel_string=kernel_string,var_names=var_names, &
                                     delta=delta,theta=thetas,periodicity=periodicities, &
                                     sparsemethod=sparsemethod,jitter=jitter,partial=.true., &
                                     factorization=factorization,verbose=verbose)
               case('val')
                 call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                                     f_r=positions(:,istart:istop),f_v=values(istart:istop),f_s=sigmaval(istart:istop), &
                                     g_r=positions(:,istart:istop),g_v=derivatives(:,istart:istop), &
                                     g_s=sigmader(:,istart:istop), &
                                     f_sparse_r=pseudo_sparse, &
                                     kernel_string=kernel_string,var_names=var_names, &
                                     delta=delta,theta=thetas,periodicity=periodicities, &
                                     sparsemethod=sparsemethod,jitter=jitter,partial=.true., &
                                     factorization=factorization,verbose=verbose)
               case('der')
                 call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                                     f_r=positions(:,istart:istop),f_v=values(istart:istop),f_s=sigmaval(istart:istop), &
                                     g_r=positions(:,istart:istop),g_v=derivatives(:,istart:istop), &
                                     g_s=sigmader(:,istart:istop), &
                                     f_sparse_r=pseudo_sparse, &
                                     kernel_string=kernel_string,var_names=var_names, &
                                     delta=delta,theta=thetas,periodicity=periodicities, &
                                     sparsemethod=sparsemethod,jitter=jitter,partial=.true., &
                                     factorization=factorization,verbose=verbose)
            end select
         else if(inval) then
            select case(sparsetype)
               case('all')
                 call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                                     f_r=positions(:,istart:istop),f_v=values(istart:istop),f_s=sigmaval(istart:istop), &
                                     f_sparse_r=pseudo_sparse, g_sparse_r=pseudo_sparse, &
                                     kernel_string=kernel_string,var_names=var_names, &
                                     delta=delta,theta=thetas,periodicity=periodicities, &
                                     sparsemethod=sparsemethod,jitter=jitter,partial=.true., &
                                     factorization=factorization,verbose=verbose)
               case('val')
                 call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                                     f_r=positions(:,istart:istop),f_v=values(istart:istop),f_s=sigmaval(istart:istop), &
                                     f_sparse_r=pseudo_sparse, &
                                     kernel_string=kernel_string,var_names=var_names, &
                                     delta=delta,theta=thetas,periodicity=periodicities, &
                                     sparsemethod=sparsemethod,jitter=jitter,partial=.true., &
                                     factorization=factorization,verbose=verbose)
               case('der')
                 call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                                     f_r=positions(:,istart:istop),f_v=values(istart:istop),f_s=sigmaval(istart:istop), &
                                     g_sparse_r=pseudo_sparse, &
                                     kernel_string=kernel_string,var_names=var_names, &
                                     delta=delta,theta=thetas,periodicity=periodicities, &
                                     sparsemethod=sparsemethod,jitter=jitter,partial=.true., &
                                     factorization=factorization,verbose=verbose)
            end select
         else if(inder) then
            select case(sparsetype)
               case('all')
                 call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                                     g_r=positions(:,istart:istop),g_v=derivatives(:,istart:istop), &
                                     g_s=sigmader(:,istart:istop), &
                                     f_sparse_r=pseudo_sparse, g_sparse_r=pseudo_sparse, &
                                     kernel_string=kernel_string,var_names=var_names, &
                                     delta=delta,theta=thetas,periodicity=periodicities, &
                                     sparsemethod=sparsemethod,jitter=jitter,partial=.true., &
                                     factorization=factorization,verbose=verbose)
               case('val')
                 call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                                     g_r=positions(:,istart:istop),g_v=derivatives(:,istart:istop), &
                                     g_s=sigmader(:,istart:istop), &
                                     f_sparse_r=pseudo_sparse, &
                                     kernel_string=kernel_string,var_names=var_names, &
                                     delta=delta,theta=thetas,periodicity=periodicities, &
                                     sparsemethod=sparsemethod,jitter=jitter,partial=.true., &
                                     factorization=factorization,verbose=verbose)
               case('der')
                 call gp_basic_teach(self=gp,kernel_type=kernel_type, &
                                     g_r=positions(:,istart:istop),g_v=derivatives(:,istart:istop), &
                                     g_s=sigmader(:,istart:istop), &
                                     g_sparse_r=pseudo_sparse, &
                                     kernel_string=kernel_string,var_names=var_names, &
                                     delta=delta,theta=thetas,periodicity=periodicities, &
                                     sparsemethod=sparsemethod,jitter=jitter,partial=.true., &
                                     factorization=factorization,verbose=verbose)
            end select
         end if

         istart = istop + 1
         istop = istop + sparsemaxn / itot
         if (istop .gt. ntot) istop = ntot 
     end do

     if(.not. partial) call gp_basic_complete(gp,factorization,verbose)

  end if

  if(verbose) call verbosity_finish_task('GP covariance calculation') 

end subroutine gp_general_main_teaching

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_prediction

  implicit none

  integer :: i, n
  integer :: alloc_failed

  if( .not. (outval .or. outvaldev .or. outder .or. outderdev .or. outhess .or. outint .or. outintdev) ) return

  ! complete gp_basic if necessary
  call gp_basic_complete(gp, factorization, verbose)

  if(verbose) call verbosity_start_task('Prediction of requested quantities') 

  if(outval) then
     allocate(val(pred_tot_grids), stat = alloc_failed)
     if( alloc_failed .ne. 0 ) then
         write(6,*) '>>> ERROR: [gp_general_main_prediction] Unable to allocate memory for val!'
         stop
     end if
  end if
  if(outvaldev) then
     allocate(valdev(pred_tot_grids), stat = alloc_failed)
     if( alloc_failed .ne. 0 ) then
         write(6,*) '>>> ERROR: [gp_general_main_prediction] Unable to allocate memory for valdev!'
         stop
     end if
  end if
  if(outder) then
     allocate(der(ndim,pred_tot_grids), stat = alloc_failed)
     if( alloc_failed .ne. 0 ) then
         write(6,*) '>>> ERROR: [gp_general_main_prediction] Unable to allocate memory for der!'
         stop
     end if
  end if
  if(outderdev) then
     allocate(derdev(ndim,pred_tot_grids), stat = alloc_failed)
     if( alloc_failed .ne. 0 ) then
         write(6,*) '>>> ERROR: [gp_general_main_prediction] Unable to allocate memory for derdev!'
         stop
     end if
  end if
  if(outhess) then
     allocate(hess(ndim,ndim,pred_tot_grids), stat = alloc_failed)
     if( alloc_failed .ne. 0 ) then
         write(6,*) '>>> ERROR: [gp_general_main_prediction] Unable to allocate memory for hess!'
         stop
     end if
  end if
  if(outint) then
     allocate(intg(pred_tot_grids), stat = alloc_failed)
     if( alloc_failed .ne. 0 ) then
         write(6,*) '>>> ERROR: [gp_general_main_prediction] Unable to allocate memory for intg!'
         stop
     end if
  end if
  if(outintdev) then
     allocate(intgdev(pred_tot_grids), stat = alloc_failed)
     if( alloc_failed .ne. 0 ) then
         write(6,*) '>>> ERROR: [gp_general_main_prediction] Unable to allocate memory for intgdev!'
         stop
     end if
  end if

  do n=1,pred_tot_grids
     if(outval) then
        val(n) = f_predict(gp,pred_grid_positions(:,n)) + aveval
     end if
     if(outvaldev) then
        valdev(n) = f_predict_var(gp,pred_grid_positions(:,n))
        if(valdev(n) .lt. 0.0_dp) then
           valdev(n) = 0.0_dp
        else
           valdev(n) = sqrt(valdev(n))
        end if
     end if
     if(outder) then
        der(:,n) = f_predict_grad(gp,pred_grid_positions(:,n))
     end if
     if(outderdev) then
        derdev(:,n) = f_predict_grad_var(gp,pred_grid_positions(:,n))
        do i=1,ndim
           if(derdev(i,n) .lt. 0.0_dp) then
              derdev(i,n) = 0.0_dp
           else
              derdev(i,n) = sqrt(derdev(i,n))
           end if
        end do
     end if
     if(outhess) then
        hess(:,:,n) = f_predict_hess(gp,pred_grid_positions(:,n))
     end if
     if(outint) then
        intg(n) = f_predict_int(gp,pred_grid_positions(:,n),interval)
     end if
     if(outintdev) then
        intgdev(n) = f_predict_int_var(gp,pred_grid_positions(:,n),interval)
        if(intgdev(n) .lt. 0.0_dp) then
           intgdev(n) = 0.0_dp
        else
           intgdev(n) = sqrt(intgdev(n))
        end if
     end if
  end do

  if(verbose) call verbosity_finish_task('Prediction of requested quantities') 

  ! if shift then find the minimum
  shiftminval = 0.0_dp
  if(outval .and. shift) then
     shiftminval=val(1)
     do n=2,pred_tot_grids
        if(val(n) < shiftminval) shiftminval = val(n)
     end do
  end if

  shiftminint = 0.0_dp
  if(outint .and. shift) then
     shiftminint=intg(1)
     do n=2,pred_tot_grids
        if(intg(n) < shiftminint) shiftminint = intg(n)
     end do
  end if

end subroutine gp_general_main_prediction

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_find_function_min

  implicit none

  integer :: i, j, k, n
  integer :: ios
  integer :: row
  integer :: alloc_failed
  logical :: file_exist
  integer :: multigrids

  if ( .not. outmin ) return

  ! prediction points
  ! from file
  if(.not. mingriddata) then

    if(verbose) call verbosity_start_task('Reading starting minimization positions from file')

    inquire(file=minposfile, exist=file_exist)
    if(.not. file_exist) then
        write(6,*) '>>> ERROR: MINPOSFILE("', trim(adjustl(minposfile)), '") does not exist!'
        stop
    end if

    open(12,file=minposfile)

    ios=0
    n=0
    do while (ios .eq. 0)
       read(12,*,iostat=ios)
       n=n+1
    end do
    n=n-1

    rewind(12)

    min_tot_grids = n
    allocate(min_grid_positions(ndim,min_tot_grids), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        write(6,*) '>>> ERROR: [gp_general_main_find_function_min] Unable to allocate memory for min_grid_positions!'
        stop
    end if

    do i=1,min_tot_grids
     if(minrownum) then
        read(12,*) row, min_grid_positions(:,i)
     else
        read(12,*) min_grid_positions(:,i)
     end if
    end do

    close(12)

    if(verbose) call verbosity_finish_task('Reading starting minimization positions from file')

  ! equidistant in the region
  else

    if (.not. minrange) then
       min_min_positions = min_positions
       min_max_positions = max_positions
       min_range_positions = range_positions
    else
       do i=1,ndim
          min_range_positions(i) = min_max_positions(i) - min_min_positions(i)
          if (min_range_positions(i) < 0) min_range_positions(i) = min_range_positions(i) + abs(periodicities(i))
       end do
    end if

    if(verbose) call verbosity_start_task('Constructing equidistant grid for minimization')

    do i=1,ndim
       min_widths(i) = min_range_positions(i) / real( min_grid(i) )
    end do

    ! calculate total number of bins
    min_tot_grids=1
    do i=1,ndim
       min_tot_grids = min_tot_grids * min_grid(i)
    end do

    ! calculate grid positions
    allocate(min_grid_positions(ndim,min_tot_grids), stat = alloc_failed)
    if( alloc_failed .ne. 0 ) then
        write(6,*) '>>> ERROR: [gp_general_main_find_function_min] Unable to allocate memory for min_grid_positions!'
        stop
    end if

    do n=1,min_tot_grids
       do i=1,ndim
          multigrids = 1
          do j=i+1,ndim
             multigrids = multigrids * min_grid(j)
          end do

          k = mod(((n-1)/multigrids),min_grid(i))

          min_grid_positions(i,n) = min_min_positions(i) + real(k)*min_widths(i) + min_widths(i) / 2.0_dp
          if(periodicities(i) .ne. 0) min_grid_positions(i,n) = min_grid_positions(i,n) &
                                       - nint(min_grid_positions(i,n)/abs(periodicities(i)))*abs(periodicities(i))
       end do
    end do

    if(verbose) call verbosity_finish_task('Constructing equidistant grid for minimization')

  end if
  

  ! initialize arrays
  allocate(minfound(min_tot_grids), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_find_function_min] Unable to allocate memory for minfound!'
      stop
  end if

  allocate(minpos(ndim,min_tot_grids), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_find_function_min] Unable to allocate memory for minpos!'
      stop                   
  end if

  allocate(minfval(min_tot_grids), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_find_function_min] Unable to allocate memory for minfval!'
      stop
  end if

  allocate(mingrad(ndim,min_tot_grids), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_find_function_min] Unable to allocate memory for mingrad!'
      stop
  end if

  select case(optimizer)
      case('lbfgs')
          call gp_general_main_min_lbfgs
  end select

end subroutine gp_general_main_find_function_min

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_min_lbfgs

  implicit none

  integer :: n, i
  real(dp) :: rmsg, gradmax, fval, lastfval, eps, xtol
  integer :: compmax
  integer :: iprint(2),iflag
  real(dp), allocatable :: work(:)
  real(dp), allocatable :: posi(:), diag(:), grad(:)
  integer :: alloc_failed
  integer :: found

  ! init variables
  lastfval = 0.0_dp
  iprint(1) = -1
  iprint(2) = 1
  eps = epsilon(0.0_dp)
  xtol = epsilon(0.0_dp)
 
  allocate(work(ndim*(2*ncorr+1)+2*ncorr), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_min_lbfgs] Unable to allocate memory for work!'
      stop
  end if

  allocate(posi(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_min_lbfgs] Unable to allocate memory for posi!'
      stop
  end if

  allocate(diag(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_min_lbfgs] Unable to allocate memory for diag!'
      stop
  end if

  allocate(grad(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_min_lbfgs] Unable to allocate memory for grad!'
      stop
  end if

  if(verbose) call verbosity_start_task('Minimization')
  do n=1,min_tot_grids
     if(verbose) write(6,'(a,i8,a)') 'Minimization from grid(',n,') has started...'
     iflag = 0
     found = 0
     posi(:) = min_grid_positions(:,n)
     diag(:) = 1.0_dp
     do i=1,maxiter-1
        ! calculate value
        fval=f_predict(gp,posi(:)) + aveval
        ! calculate gradient
        grad(:)=f_predict_grad(gp,posi(:))
        ! calculate maxgrad and rmsg
        rmsg = grad_rmsg(grad,gradmax,compmax) 
        ! check function value change
        if( (i .ne. 1) .and. (abs(fval-lastfval) .le. maxdval) ) then
            if(verbose) then
               write(6,'(a,f15.8)') '>>> INFO: RMS of gradient: ', rmsg
               write(6,'(a,f15.8)') '>>> INFO: RMS of gradient treshold: ', maxrmsgrad
               write(6,'(a,f15.8)') '>>> INFO: Max gradient component: ', abs(gradmax)
               write(6,'(a,f15.8)') '>>> INFO: Max gradient component treshold: ', maxgrad
               write(6,'(a,f15.8)') '>>> INFO: Last function value change: ', abs(fval-lastfval)
               write(6,'(a,f15.8)') '>>> INFO: Function value change treshold: ', maxdval
               write(6,'(a)') '>>> INFO: Function value change is below trashold! Minimization is stopped.'
            end if
            found=i
            exit
        end if
        ! check gradient
        if( (abs(gradmax) .le. maxgrad) .and. (rmsg .le. maxrmsgrad) ) then
            if(verbose) then
               write(6,'(a,f15.8)') '>>> INFO: RMS of gradient: ', rmsg
               write(6,'(a,f15.8)') '>>> INFO: RMS of gradient treshold: ', maxrmsgrad
               write(6,'(a,f15.8)') '>>> INFO: Max gradient component: ', abs(gradmax)
               write(6,'(a,f15.8)') '>>> INFO: Max gradient component treshold: ', maxgrad
               write(6,'(a,f15.8)') '>>> INFO: Last function value change: ', abs(fval-lastfval)
               write(6,'(a,f15.8)') '>>> INFO: Function value change treshold: ', maxdval
               write(6,'(a)') '>>> INFO: Gradient tresholds were satisfied! Minimization is stopped.'
            end if
            found=i
            exit
        end if
        ! print intermediate results
        if(verbose) then
           write(6,'(a,f15.8)') '>>> INFO: RMS of gradient: ', rmsg
           write(6,'(a,f15.8,a,i8)') '>>> INFO: Max gradient component: ', abs(gradmax), ' component: ', compmax
           write(6,'(a,f15.8)') '>>> INFO: Function value: ', fval
        end if
        ! lbfgs step
        call lbfgs(ndim,ncorr,posi,fval,grad,.false.,diag,iprint,eps,xtol,work,iflag)
        ! wrap if variables are periodic
        call positions_wrap(posi,abs(periodicities))
        if( iflag .eq. 0) exit
        if( iflag .le. 0) then
            write(6,*) '>>> ERROR: Internal L-BFGS driver error! Code = ', iflag
            found=iflag
            exit
        end if
        lastfval = fval
     end do
     ! store results
     minfound(n) = found
     minpos(:,n) = posi(:)
     minfval(n) = lastfval
     mingrad(:,n) = grad(:)
     if(verbose) write(6,'(a,i8,a)') 'Minimization from grid(',n,') has finished...'
  end do

  if(verbose) call verbosity_finish_task('Minimization')

end subroutine gp_general_main_min_lbfgs

!--------------------------------------------------------------------------------------------------

subroutine gp_general_main_find_path

  implicit none

  integer :: alloc_failed

  if ( .not. outpath ) return

  ! allocate initial paths
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_main_opt_lbfgs] Unable to allocate memory for work!'
      stop
  end if

end subroutine gp_general_main_find_path

!--------------------------------------------------------------------------------------------------

real(dp) function grad_rmsg(grad,gradmax,compmax)

  implicit none

  real(dp) :: grad(:)
  real(dp) :: gradmax
  integer :: compmax

  integer :: i, nd
  real(dp) :: norm

  gradmax = 0.0
  grad_rmsg = 0.0
  compmax = 0

  nd = size(grad)

  if(nd .eq. 0) return

  do i=1, nd
     norm = grad(i)**2
     if( abs(grad(i)) > abs(gradmax) ) then
         gradmax = grad(i)
         compmax = i
     end if
     grad_rmsg = grad_rmsg + norm 
  end do

  grad_rmsg = sqrt(grad_rmsg/real(nd))

  return

end function grad_rmsg

!--------------------------------------------------------------------------------------------------

end module gp_general_main_mod
