module gp_general_read_mod

use verbosity_mod
use fparser_mod
use gp_basic_dat_mod
use gp_basic_mod
use gp_general_dat_mod

implicit none

contains

!--------------------------------------------------------------------------------------------------

subroutine gp_general_read_core

  implicit none

  call gp_general_read_help
  call gp_general_read_args_type
  call gp_general_read_gp_files
  call gp_general_read_args_data
  call gp_general_read_alloc_var
  call gp_general_read_args_rest

end subroutine gp_general_read_core

!--------------------------------------------------------------------------------------------------

subroutine gp_general_read_help

  implicit none

  logical :: help
  logical :: help_input, help_output, help_basic, help_sparse, help_hyperopt, help_minim, help_optim, help_format
  integer :: i
  character(len=length) :: argi

  help = .false.
  help_input = .false.
  help_output = .false.
  help_basic = .false.
  help_sparse = .false.
  help_hyperopt = .false.
  help_minim = .false.
  help_optim = .false.
  help_format = .false.

  argcount = command_argument_count()

  if ( argcount .eq. 0 ) then
     help=.true.
  else
     do i=1,argcount
        call get_command_argument(i,argi)

        if (index(uppercase(argi),'-H ') .ne. 0) then
           help=.true.
           exit
        else if(index(uppercase(argi),'-HELP ') .ne. 0) then
           help=.true.
           exit
        else if(index(uppercase(argi),'--HELP ') .ne. 0) then
           help=.true.
           exit
        end if
     end do
  end if

  if ( help ) then
     write(6,*)
     write(6,10) 'GP general program'
     write(6,10) '   Available topics:'
     write(6,10) '     -help_input     => description of available input options'
     write(6,10) '     -help_output    => description of available output options'
     write(6,10) '     -help_basic     => description of basic GP options'
     write(6,10) '     -help_sparse    => description of sparse GP options'
     write(6,10) '     -help_hyperopt  => description of options of hyper parameter optimization'
     write(6,10) '     -help_minim     => description of minimization related options'
     write(6,10) '     -help_optim     => description of optimizer related options'
     write(6,10) '     -help_format    => description of file formats'
     write(6,10) ''
     write(6,10) '     -help_all       => description of all options'
     write(6,10) ''

     stop
  end if

  do i=1,argcount
     call get_command_argument(i,argi)

     if (index(uppercase(argi),'-HELP_ALL ') .ne. 0) then
        help_input=.true.
        help_output=.true.
        help_basic=.true.
        help_sparse=.true.
        help_hyperopt=.true.
        help_minim=.true.
        help_optim=.true.
        help_format=.true.
        exit
     end if

     if (index(uppercase(argi),'-HELP_INPUT ') .ne. 0) then
        help_input=.true.
     end if

     if(index(uppercase(argi),'-HELP_OUTPUT ') .ne. 0) then
        help_output=.true.
     end if

     if(index(uppercase(argi),'-HELP_BASIC ') .ne. 0) then
        help_basic=.true.
     end if

     if (index(uppercase(argi),'-HELP_SPARSE ') .ne. 0) then
        help_sparse=.true.
     end if

     if(index(uppercase(argi),'-HELP_HYPEROPT ') .ne. 0) then
        help_hyperopt=.true.
     end if

     if(index(uppercase(argi),'-HELP_MINIM ') .ne. 0) then
        help_minim=.true.
     end if

     if(index(uppercase(argi),'-HELP_OPTIM ') .ne. 0) then
        help_optim=.true.
     end if

     if(index(uppercase(argi),'-HELP_FORMAT ') .ne. 0) then
        help_format=.true.
     end if
  end do

  if ( help_input ) then
     write(6,10) ''
     write(6,10) 'Input options:'
     write(6,10) '--------------'
     write(6,10) '  -INVAL (function values are available, default=.false.)'
     write(6,10) '  -INDER (function derivatives are available, default=.false.)'
     write(6,10) '  -INGP (one or more gp files are available for reading informations, default=.false.)'
!    write(6,10) '  -DOCHANGE (change position of teaching points of existing teaching, default=.false.)'
     write(6,10) ''
  end if

  if ( help_output ) then
     write(6,10) ''
     write(6,10) 'Output options:'
     write(6,10) '---------------'
     write(6,10) '  -OUTVAL (calculate values, default=.false.)'
     write(6,10) '  -OUTVALDEV (calculate standard deviation of values, default=.false.)'
     write(6,10) '  -OUTDER (calculate derivatives, default=.false.)'
     write(6,10) '  -OUTDERDEV (calculate standard deviation of derivatives, default=.false.)'
     write(6,10) '  -OUTINT (calculate integrals, default=.false.)'
     write(6,10) '  -OUTINTDEV (calculate standard deviation of integrals, default=.false.)'
     write(6,10) '  -OUTHESS (calculate hessian, default=.false.)'
     write(6,10) '  -OUTGP (gp file is generated, default=.false.)'
     write(6,10) '  -OUTOPT (numerical optimization of hyperparameters, default=.false.)'
     write(6,10) '  -OUTMIN (find minima, default=.false.)'
     write(6,10) ''
  end if

  if ( help_basic ) then
     write(6,10) ''
     write(6,10) 'Basic options:'
     write(6,10) '--------------'
     write(6,10) '  -NDIM 3 (number of dimensions default=1)'
     write(6,10) '  -PERIODICITY 6.28318530718 0 -2.0 ... (define periodicities, =0: nonperiodic, <0: naive, default=0. 0. ...)'
     write(6,10) '  -TRAINFILE trainfile (name of training/input file, default=training.txt)'
     write(6,10) '  -PREDFILE predfile (name of prediction/output file, default=prediction.txt)'
     write(6,10) '  -LSTART 100 (starting line of data to be used, default=1)'
     write(6,10) '  -LSTOP 2000 (final line of data to be used, default=0)'
     write(6,10) '  -LINC 10 (increment number of line of data to be used, default=1)'
     write(6,10) '  -ROWNUM (first column is row numbers, default=.false.)'
     write(6,10) '  -AVEVAL -10.11 (average value of the function, default=0.0)'
     write(6,10) '  -KERNEL (type of kernel: se, user, default=se)'
     write(6,10) '  -DELTA 7.0 (delta value, default=1.0)'
     write(6,10) '  -THETA 1.0 3.14 ... (theta values)'
     write(6,10) '  -FKERNEL "..." (kernel function)'
     write(6,10) '  -PREFIX1 prefix (prefix of first variables in kernel function, default=x)'
     write(6,10) '  -PREFIX2 prefix (prefix of second variables in kernel function, default=y)'
     write(6,10) '  -SIGMAVAL 4.0 (standard deviation of function values, default reads it from file)'
     write(6,10) '  -SIGMADER 25.0 25.0 ... (standard deviation of derivatives, default reads it from file)'
     write(6,10) '  -JITTER jitter (jitter for covariance matrix inversion, default=1.0e-14)'
     write(6,10) '  -RANGE -3.1415926536 3.1415926536 -3. 5. 1. 10. ... (define universal region, default=0. 0. ...)'
     write(6,10) '  -PREDGRID 50 50 ... (grid sizes for prediciton, ie. prediction at these points, default=10 10 ...)'
     write(6,10) '  -PREDRANGE -3.1415926536 3.1415926536 -3. 5. 1. 10. ... (define grid region for prediction, default=0. 0. ...)'
     write(6,10) '  -PREDPOSFILE predposfile (read prediction positions from file)'
     write(6,10) '  -PREDROWNUM (first column is row number in predposfile, default=.false.)'
     write(6,10) '  -GPINFILE /path/to/gpinfile1 -GPINFILE /path/to/gpinfile2 ... (specification of input gp files)'
     write(6,10) '  -SHIFT (shift value profile starting from 0, default=.false.)'
     write(6,10) '  -GNUPLOT (gnuplot  format of the output, default=.false.)'
     write(6,10) '  -VERBOSE (increase information in output, default=.false.)'
     write(6,10) '  -GPFORMAT (gp print out mode: full, partial, compact, default=full)'
     write(6,10) '  -GPOUTFILE gpoutfile (name of gp file including all necessary information, default=gp.txt)'
     write(6,10) '  -INTERVAL -3.1415926536 3.1415926536 0 0 ... (define intervals for integration, default=0. 0. ...)'
     write(6,10) '  -FACTORIZATION qr (type of factorization method: ch, qr, bk, default=ch)'
     write(6,10) ''
  end if

  if ( help_sparse ) then
     write(6,10) ''
     write(6,10) 'Sparsification options:'
     write(6,10) '-----------------------'
     write(6,10) '  -PARTIAL (do a partial sparsification type calculation, default=.false.)'
     write(6,10) '  -SPARSIFICATION clustering (sparsification type: none, clustering, grid and file; default=none)'
     write(6,10) '  -SPARSEMETHOD fitc (sparsification method: dic, dtc, fitc; default=dtc)'
     write(6,10) '  -SPARSETYPE der (type of sparse points in sparsification: val, der, all; default=val)'
     write(6,10) '  -CLUSTERING kmeans (clustering method for sparsification: every, kmeans, kmeans++; default=every)'
     write(6,10) '  -SPARSEPOINTS 1000 (number of sparse points for clustering sparsification, default=100)'
     write(6,10) '  -SPARSEFREQ 100 (sparse frequency from data for clustering sparsification, default=1)'
     write(6,10) '  -SPARSERANGE -3.1415926536 3.1415926536 -3. 5. ... (define grid region for sparsification, default=0. 0. ...)'
     write(6,10) '  -SPARSEGRID 10 20 ... (grid sizes for grid sparsification, ie. sparse points are on a grid, default=10 10 ...)'
     write(6,10) '  -SPARSEPOSFILE sparseposfile (read sparse positions from file)'
     write(6,10) '  -SPARSEROWNUM (first column is row numbers in sparseposfile, default=.false.)'
     write(6,10) '  -SPARSEMAXN (maximum number of data points for iterative sparse calculation, default=0)'
     write(6,10) '  -KMEANSTHRESHOLD 0.1 0.01 0. ... (define thresholds for kmeans, default=0. 0. ...)'
     write(6,10) ''
  end if

  if ( help_hyperopt ) then
     write(6,10) ''
     write(6,10) 'Hyper parameter optimization options:'
     write(6,10) '-------------------------------------'
     write(6,10) '  -OPTDELTA 0.01 (relative stepsize for delta optimization, =0: no optim, default=0)'
     write(6,10) '  -OPTSIGMAVAL 0.001 (relative stepsize for value standard deviation optimization, =0: no optim, default=0)'
     write(6,10) '  -OPTSIGMADER 0.02 0.01 0 ... (rel. stepsizes for der st. deviation optim, =0: no optim, default=0 0 0 ...)'
     write(6,10) '  -OPTTHETA 0.1 0. 0.001  ... (rel. stepsizes for theta values optim, =0: no optim, default=0 0 0 ...)'
     write(6,10) '  -OPTFILE optfile (name of opt file including hypers optimum, default=opt.txt)'
     write(6,10) '  -NUMDIFF two-sided (specifies numerical differentiation: one-sided, two-sided, default=one-sided)'
     write(6,10) ''
  end if

  if ( help_minim ) then
     write(6,10) ''
     write(6,10) 'Function minimization options:'
     write(6,10) '------------------------------'
     write(6,10) '  -MINFILE minfile (name of min file including minima, default=min.txt)'
     write(6,10) '  -MINRANGE -3.1415926536 3.1415926536 -3. 5. 1. 10. ... (define grid region for minimization, default=0. 0. ...)'
     write(6,10) '  -MINGRID 10 20 ... (grid sizes for minimization, ie. starting minims from these points, default=10 10 ...)'
     write(6,10) '  -MINPOSFILE minposfile (read minimization starting positions from file)'
     write(6,10) '  -MINROWNUM (first column is row numbers in minposfile, default=.false.)'
     write(6,10) ''
  end if

  if ( help_optim ) then
     write(6,10) ''
     write(6,10) 'Optimizer options:'
     write(6,10) '---------------------'
     write(6,10) '  -SEED 12345 (define random seed, default=3141592)'
     write(6,10) '  -OPTIMIZER lbfgs (specifies optimizer for minimization/optimization, default=lbfgs)'
     write(6,10) '  -NCORR 5 (number of corrections for lbfgs optimizer, default=4)'
     write(6,10) '  -MAXITER 200 (maximum number of iterations for minimization, default=100)'
     write(6,10) '  -PRINTFREQ 10 (frequency of printing out results during minimization, default=1)'
     write(6,10) '  -MAXDVAL 0.1 (max change in value for minimization, default=0.000001)'
     write(6,10) '  -MAXGRAD 0.5 (max gradient component for minimization, default=0.0001)'
     write(6,10) '  -MAXRMSGRAD 1.0 (max rms of gradients for minimization, default=0.001)'
     write(6,10) ''
  end if

  if ( help_format ) then
     write(6,10) ''
     write(6,10) 'File formats:'
     write(6,10) '-------------'
     write(6,10) '  -> trainfile format: '
     write(6,10) '     [ROWNUM] POS(1) ... POS(ndim) [VAL] [SIGMAVAL] [DER(1) ... DER(ndim)] [SIGMADER(1) ... SIGMADER(ndim)]'
     write(6,10) ''
     write(6,10) '  -> predposfile/sparseposfile/minposfile format:'
     write(6,10) '     [ROWNUM] POS(1) ... POS(ndim)'  
     write(6,10) ''
     write(6,10) '  -> predfile format: '
     write(6,10) '     POS(1) ... POS(ndim) [VAL] [VALDEV] [DER(1) ... DER(ndim)] [DERDEV(1) ... DERDEV(ndim)]'
     write(6,10) '     [HESS(1,1) ... HESS(1,ndim) ... HESS(ndim,1) ... HESS(ndim,ndim)]'
     write(6,10) ''
  end if

  if (help_input .or. help_output .or. help_basic .or. help_sparse .or. &
      help_hyperopt .or. help_minim .or. help_optim .or. help_format) stop

  10 format(100A)

end subroutine gp_general_read_help

!--------------------------------------------------------------------------------------------------

subroutine gp_general_read_args_type

  implicit none
  
  integer :: i
  character(len=length) :: argi

  ! read type(s) of analysis
  do i=1,argcount
     call get_command_argument(i,argi)

     if (index(uppercase(argi),'-INVAL ') .ne. 0) then
        inval=.true.
     end if

     if (index(uppercase(argi),'-INDER ') .ne. 0) then
        inder=.true.
     end if

     if (index(uppercase(argi),'-INGP ') .ne. 0) then
        ingp=.true.
     end if
  end do

  if ( .not. (inval .or. inder .or. ingp) ) then
       write(6,*)
       write(6,*) '>>> ERROR:'
       write(6,*) 'At least one of the following commands must be specified:'
       write(6,*) '-INVAL -INDER'
       write(6,*) 'or:'
       write(6,*) '-INGP'
       write(6,*)
       stop
  end if

  do i=1,argcount
     call get_command_argument(i,argi)

     if (index(uppercase(argi),'-OUTVAL ') .ne. 0) then
        outval=.true.
     end if

     if (index(uppercase(argi),'-OUTVALDEV ') .ne. 0) then
        outvaldev=.true.
     end if

     if (index(uppercase(argi),'-OUTDER ') .ne. 0) then
        outder=.true.
     end if

     if (index(uppercase(argi),'-OUTDERDEV ') .ne. 0) then
        outderdev=.true.
     end if

     if (index(uppercase(argi),'-OUTHESS ') .ne. 0) then
        outhess=.true.
     end if

     if (index(uppercase(argi),'-OUTINT ') .ne. 0) then
        outint=.true.
     end if

     if (index(uppercase(argi),'-OUTINTDEV ') .ne. 0) then
        outintdev=.true.
     end if

     if (index(uppercase(argi),'-OUTGP ') .ne. 0) then
        outgp=.true.
     end if

     if (index(uppercase(argi),'-OUTOPT ') .ne. 0) then
        outopt=.true.
     end if

     if (index(uppercase(argi),'-OUTMIN ') .ne. 0) then
        outmin=.true.
     end if

     if (index(uppercase(argi),'-OUTPATH ') .ne. 0) then
        outpath=.true.
     end if
  end do

  if ( .not. (outval .or. &
              outvaldev .or. &
              outder .or. &
              outderdev .or. &
              outhess .or. &
              outint .or. &
              outintdev .or. &
              outgp .or. &
              outopt .or. &
              outmin .or. &
              outpath) ) then
       write(6,*)
       write(6,*) '>>> ERROR:'
       write(6,*) 'At least one of the following commands must be specified:'
       write(6,*) '-OUTVAL -OUTVALDEV -OUTDER -OUTDERDEV -OUTHESS -OUTINT -OUTINTDEV -OUTGP -OUTOPT -OUTMIN -OUTPATH'
       write(6,*)
       stop
  end if

  if ( ingp .and. outopt ) then
       write(6,*)
       write(6,*) '>>> ERROR:'
       write(6,*) '-OUTOPT cannot be specified with -INGP at the same time!'
       write(6,*)
       stop
  end if

  if ( outmin .and. outopt ) then
       write(6,*)
       write(6,*) '>>> ERROR:'
       write(6,*) '-OUTMIN cannot be specified with -OUTOPT at the same time!'
       write(6,*)
       stop
  end if

  if ( outpath .and. outopt ) then
       write(6,*)
       write(6,*) '>>> ERROR:'
       write(6,*) '-OUTPATH cannot be specified with -OUTOPT at the same time!'
       write(6,*) 
       stop
  end if

  if ( outmin .and. outpath ) then
       write(6,*)
       write(6,*) '>>> ERROR:'
       write(6,*) 'Either -OUTMIN or -OUTPATH can be specified at the same time!'
       write(6,*)
       stop
  end if

  do i=1,argcount
     call get_command_argument(i,argi)

     if (index(uppercase(argi),'-DOCHANGE ') .ne. 0) then
        dochange=.true.
     end if
  end do

  if ( dochange ) then
     if ( .not. ingp ) then
        write(6,*) '>>> ERROR:'
        write(6,*) '-DOCHANGE requires the option -INGP!'
        stop
     end if
  end if 

  do i=1,argcount
     call get_command_argument(i,argi)

     if (index(uppercase(argi),'-KERNEL ') .ne. 0) then
        call get_command_argument(i+1,argi)
        select case(trim(adjustl(lowercase(argi))))
           case('user')
              kernel_type = 0
           case('se')
              kernel_type = 1
           case default
              write(6,*) '>>> ERROR: KERNEL("', trim(adjustl(argi)), &
                          '") must be one of the followings: se, user!'
              stop
        end select
     end if

     if (index(uppercase(argi),'-AVEVAL ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*) aveval
     end if

     if (index(uppercase(argi),'-GPFORMAT ') .ne. 0) then
        call get_command_argument(i+1,argi)
        argi=trim(adjustl(lowercase(argi)))
        read(argi,*) gpformat
        if( (trim(adjustl(gpformat)) .ne. 'full') .and. &
            (trim(adjustl(gpformat)) .ne. 'partial') .and. &
            (trim(adjustl(gpformat)) .ne. 'compact') ) then
            write(6,*) '>>> ERROR: GPFORMAT("', trim(adjustl(gpformat)), &
                          '") must be one of the followings: full, partial, compact!'
            stop
        end if
     end if

     if (index(uppercase(argi),'-VERBOSE ') .ne. 0) then
        verbose=.true.
     end if
  end do

end subroutine gp_general_read_args_type

!--------------------------------------------------------------------------------------------------

subroutine gp_general_read_gp_files

  implicit none

  integer :: i, j
  character(len=length) :: argi
  integer :: alloc_failed
  type(gp_basic) :: gp_file
  logical :: file_exist

  if( .not. ingp) return

  do i=1,argcount
     call get_command_argument(i,argi)

     if (index(uppercase(argi),'-GPINFILE ') .ne. 0) then
        ngpinfiles = ngpinfiles + 1
     end if
  end do

  if ( ngpinfiles .eq. 0 ) then
     write(6,*) ">>> ERROR: -INGP was defined but no GPINFILE was specified!"
     stop
  end if

  allocate(gpinfiles(ngpinfiles), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_gp_files] Unable to allocate memory for gpinfiles!'
      stop
  end if

  j = 0
  do i=1,argcount
     call get_command_argument(i,argi)

     if (index(uppercase(argi),'-GPINFILE ') .ne. 0) then
        j = j + 1
        call get_command_argument(i+1,gpinfiles(j))
        inquire(file=gpinfiles(j), exist=file_exist)
        if(.not. file_exist) then
           write(6,*) '>>> ERROR: GPINFILE("', trim(adjustl(gpinfiles(j))), '") does not exist!'
           stop
        end if
     end if
  end do

  if(verbose) call verbosity_start_task('Reading gpinfile(s)')

  ! read first gp file => gp
  if(verbose) write(6,*) 'reading gpinfile(',1,')=',trim(adjustl(gpinfiles(1))),' file has started'
  call gp_basic_read(self=gp,fileid=5,filename=gpinfiles(1))
  if(verbose) write(6,*) 'reading gpinfile(',1,')=',trim(adjustl(gpinfiles(1))),' file has finished'

  ! read the other gp files (if any) and merge them with gp
  do i=2,ngpinfiles
     if(verbose) write(6,*) 'reading gpinfile(',i,')=',trim(gpinfiles(i)),' file has started'
     call gp_basic_read(self=gp_file,fileid=5,filename=gpinfiles(i))
     if(verbose) write(6,*) 'reading gpinfile(',i,')=',trim(gpinfiles(i)),' file has finished'
     if(verbose) write(6,*) 'merging gpinfile(',i,')=',trim(gpinfiles(i)),' file has started'
     call gp_basic_merge(gp_file, gp)
     if(verbose) write(6,*) 'merging gpinfile(',i,')=',trim(gpinfiles(i)),' file has finished'
  end do

  if(verbose) call verbosity_finish_task('Reading gpinfile(s)')

end subroutine gp_general_read_gp_files

!--------------------------------------------------------------------------------------------------

subroutine gp_general_read_args_data

  implicit none

  integer :: i, j
  character(len=length) :: argi, argj
  integer :: ios
  integer :: row
  integer :: alloc_failed
  logical :: file_exist

  if( .not. (inval .or. inder) ) return

  do i=1,argcount
     call get_command_argument(i,argi)

     if (index(uppercase(argi),'-TRAINFILE ') .ne. 0) then
        call get_command_argument(i+1,trainfile)
        inquire(file=trainfile, exist=file_exist)
        if(.not. file_exist) then
           write(6,*) '>>> ERROR: TRAINFILE("', trim(adjustl(trainfile)), '") does not exist!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-NDIM ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) ndim
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: NDIM("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(ndim .le. 0) then
           write(6,*) '>>> ERROR: NDIM(', ndim, ') must be larger than 0!'
           stop
        end if

        if(ingp) then
          if(ndim .ne. gp%n_dof) then
             write(6,*) '>>> ERROR: ndim(', ndim, ') must be equal to n_dof(', gp%n_dof, ')!'
             stop
          end if
        end if
     end if

     if (index(uppercase(argi),'-LSTART ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*, iostat=ios) lstart
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: LSTART("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(lstart .le. 0) then
           write(6,*) '>>> ERROR: LSTART(', lstart, ') must be larger than 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-LSTOP ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) lstop
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: LSTOP("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(lstop .le. 0) then
           write(6,*) '>>> ERROR: LSTOP(', lstop, ') must be larger than 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-LINC ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) linc
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: LINC("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(linc .le. 0) then
           write(6,*) '>>> ERROR: LINC(', linc, ') must be larger than 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-ROWNUM ') .ne. 0) then
        rownum=.true.
     end if
  end do

  if(verbose) call verbosity_start_task('Reading data file')

  inquire(file=trainfile, exist=file_exist)
  if(.not. file_exist) then
      write(6,*) '>>> ERROR: TRAINFILE("', trim(adjustl(trainfile)), '") does not exist!'
      stop
  end if

  open(5,file=trainfile)

  ! read data file
  ios=0
  nline=0
  do while (ios .eq. 0)
     read(5,*,iostat=ios)
     nline=nline+1
  end do
  nline=nline-1

  rewind(5)

  if ( (lstop .eq. 0) .or. (lstop .gt. nline) ) lstop=nline
  if ( lstart .gt. lstop ) then
       write(6,*) '>>> ERROR: LSTART(', lstart, ') must be smaller than LSTOP(', lstop, ')!'
       stop
  end if
  if ( linc .gt. (lstop-lstart+1) ) then
       write(6,*) '>>> ERROR: LINC(', linc, ') must be smaller than LSTOP-LSTART+1(', lstop-lstart+1, ')!'
       stop
  end if
  ntot=(lstop-lstart)/linc+1

  if( inval .and. inder ) then
      itot = ndim+1
      stot = ntot*(ndim+1)
  else if( inval ) then
      itot = 1
      stot = ntot
  else if( inder ) then
      itot = ndim
      stot = ntot*ndim
  end if

  allocate(positions(ndim,ntot), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_args_data] Unable to allocate memory for positions!'
      stop
  end if

  if( inval ) then
      allocate(values(ntot), stat = alloc_failed)
      if( alloc_failed .ne. 0 ) then
          write(6,*) '>>> ERROR: [gp_general_read_args_data] Unable to allocate memory for values!'
          stop
      end if
      allocate(sigmaval(ntot), stat = alloc_failed)
      if( alloc_failed .ne. 0 ) then
          write(6,*) '>>> ERROR: [gp_general_read_args_data] Unable to allocate memory for sigmaval!'
          stop
      end if
  end if

  if( inder ) then
      allocate(derivatives(ndim,ntot), stat = alloc_failed)
      if( alloc_failed .ne. 0 ) then
          write(6,*) '>>> ERROR: [gp_general_read_args_data] Unable to allocate memory for derivatives!'
          stop
      end if
      allocate(sigmader(ndim,ntot), stat = alloc_failed)
      if( alloc_failed .ne. 0 ) then
          write(6,*) '>>> ERROR: [gp_general_read_args_data] Unable to allocate memory for sigmader!'
          stop
      end if
  end if

  ! if variances are defined in command line then take them
  do i=1,argcount
     call get_command_argument(i,argi)

     if (index(uppercase(argi),'-SIGMAVAL ') .ne. 0) then
        commandsigmaval=.true.
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) sigmaval(1)
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: SIGMAVAL("', trim(adjustl(argi)), '") has a wrong fromat!'
           stop
        end if
        if(sigmaval(1) .lt. 0) then
           write(6,*) '>>> ERROR: SIGMAVAL(', sigmaval(1), ') must be nonnegative!'
           stop
        end if
        sigmaval(2:ntot) = sigmaval(1)
     end if

     if (index(uppercase(argi),'-SIGMADER ') .ne. 0) then
        commandsigmader=.true.
        do j=1,ndim
           call get_command_argument(i+j,argj)
           read(argj,*,iostat=ios) sigmader(j,1)
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: ', j, 'th component of SIGMADER("', trim(adjustl(argj)), '") has a wrong format!'
           end if
           if(sigmader(j,1) .lt. 0) then
              write(6,*) '>>> ERROR: ', j, 'th component of SIGMADER(', sigmader(j,1), ') must be nonnegative!'
              stop
           end if
           sigmader(j,2:ntot) = sigmader(j,1)
        end do
     end if
  end do

  j=0
  do i=1,lstop
     if(i .lt. lstart) then
        read(5,*)
        cycle
     end if
     if(mod((i-lstart),linc) .ne. 0) then
        read(5,*)
        cycle
     end if 
     j=j+1
     if(inval .and. inder) then
        if(rownum) then
           if(commandsigmaval .and. commandsigmader) then
             read(5,*) row, positions(:,j), values(j), derivatives(:,j)
           else if(commandsigmaval) then
             read(5,*) row, positions(:,j), values(j), sigmaval(j), derivatives(:,j)
           else if(commandsigmader) then
             read(5,*) row, positions(:,j), values(j), derivatives(:,j), sigmader(:,j)
           else
             read(5,*) row, positions(:,j), values(j), sigmaval(j), derivatives(:,j), sigmader(:,j)
           end if
        else
           if(commandsigmaval .and. commandsigmader) then
             read(5,*) positions(:,j), values(j), derivatives(:,j)
           else if(commandsigmaval) then
             read(5,*) positions(:,j), values(j), sigmaval(j), derivatives(:,j)
           else if(commandsigmader) then
             read(5,*) positions(:,j), values(j), derivatives(:,j), sigmader(:,j)
           else
             read(5,*) positions(:,j), values(j), sigmaval(j), derivatives(:,j), sigmader(:,j)
           end if
        end if
     else if(inval) then
        if(rownum) then
           if(commandsigmaval) then
             read(5,*) row, positions(:,j), values(j)
           else
             read(5,*) row, positions(:,j), values(j), sigmaval(j)
           end if
        else
           if(commandsigmaval) then
             read(5,*) positions(:,j), values(j)
           else
             read(5,*) positions(:,j), values(j), sigmaval(j)
           end if
        end if
     else if(inder) then
        if(rownum) then
           if(commandsigmader) then
             read(5,*) row, positions(:,j), derivatives(:,j)
           else
             read(5,*) row, positions(:,j), derivatives(:,j), sigmader(:,j)
           end if
        else
           if(commandsigmader) then
             read(5,*) positions(:,j), derivatives(:,j)
           else
             read(5,*) positions(:,j), derivatives(:,j), sigmader(:,j)
           end if
        end if
     end if
     if(j .eq. ntot) exit
  end do

  close(5)

  if(verbose) call verbosity_finish_task('Reading data file')

  ! shift values by aveval
  if(allocated(values)) then
     values = values - aveval
  end if

end subroutine gp_general_read_args_data

!--------------------------------------------------------------------------------------------------

subroutine gp_general_read_alloc_var

  implicit none

  integer :: alloc_failed

  if ( ingp ) then
       ndim = gp%n_dof
       kernel_type = gp%kernel_type
  end if

  select case (kernel_type)
      case(USER_kernel)
         allocate(var_names(2*ndim), stat = alloc_failed)
         if( alloc_failed .ne. 0 ) then
             write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for var_names!'
             stop
         end if

      case(SE_kernel)
         allocate(thetas(ndim), stat = alloc_failed)
         if( alloc_failed .ne. 0 ) then
             write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for thetas!'
             stop
         end if

         thetas(:)=0.0_dp

         allocate(optsigmader(ndim), stat = alloc_failed)
         if( alloc_failed .ne. 0 ) then
             write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for optsigmader!'
             stop
         end if

         optsigmader(:)=0.0_dp

         allocate(opttheta(ndim), stat = alloc_failed)
         if( alloc_failed .ne. 0 ) then
             write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for opttheta!'
             stop
         end if

         opttheta(:)=0.0_dp

  end select

  allocate(periodicities(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for periodicities!'
      stop
  end if

  periodicities(:)=0.0_dp

  allocate(min_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for min_positions!'
      stop
  end if

  allocate(max_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for max_positions!'
      stop
  end if

  allocate(range_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for range_positions!'
      stop
  end if

  allocate(pred_grid(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for pred_grid!'
      stop
  end if

  pred_grid(:)=10

  allocate(pred_min_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for pred_min_positions!'
      stop
  end if

  allocate(pred_max_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for pred_max_positions!'
      stop
  end if

  allocate(pred_range_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for pred_range_positions!'
      stop
  end if

  allocate(pred_widths(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for pred_widths!'
      stop
  end if

  allocate(sparse_grid(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for sparse_grid!'
      stop
  end if

  sparse_grid(:)=10

  allocate(sparse_min_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for sparse_min_positions!'
      stop
  end if

  allocate(sparse_max_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for sparse_max_positions!'
      stop
  end if

  allocate(sparse_range_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for sparse_range_positions!'
      stop
  end if

  allocate(sparse_widths(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for sparse_widths!'
      stop
  end if

  allocate(min_grid(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for min_grid!'
      stop
  end if

  min_grid(:)=10

  allocate(min_min_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for min_min_positions!'
      stop
  end if

  allocate(min_max_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for min_max_positions!'
      stop
  end if

  allocate(min_range_positions(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for min_range_positions!'
      stop
  end if

  allocate(min_widths(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for min_widths!'
      stop
  end if

  allocate(kmeansthresholds(ndim), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_alloc_var] Unable to allocate memory for kmeansthresholds!'
      stop
  end if

  kmeansthresholds(:)=0.0_dp

  if ( ingp ) then
     select case(kernel_type)
        case(USER_kernel)

        case(SE_kernel)
          thetas = sqrt(gp%theta_sq)
          periodicities = gp%periodicity
          delta = sqrt(gp%delta_sq)
     end select
     if (gp%sparsified) then
        sparsification = 'gpfile'
        sparsemethod = gp%sparsemethod
        if ( (gp%m_f .ne. 0) .and. (gp%m_g .ne. 0) ) then
           sparsetype = 'all'
        else if ( gp%m_f .ne. 0 ) then
           sparsetype = 'val'
        else if ( gp%m_g .ne. 0 ) then
           sparsetype = 'der'
        else
           write(6,*) '>>> ERROR: [gp_general_read_alloc_var] cannot define sparsetype!'
           stop
        end if
     end if
  end if

  allocate(interval(ndim,2), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_read_alloc_var] Unable to allocate memory for interval!'
      stop
  end if

  interval(:,:) = 0

end subroutine gp_general_read_alloc_var

!--------------------------------------------------------------------------------------------------

subroutine gp_general_read_args_rest

  implicit none

  integer :: i, j
  integer :: ios
  character(len=length) :: argi, argj
  logical :: file_exist
  integer :: alloc_failed
  integer :: isize
  integer, allocatable :: iseed(:)
  character(len=length) :: c

  if ( .not. ingp ) then
     select case(kernel_type)
        case(USER_kernel)
           do i=1,argcount
              call get_command_argument(i,argi)

              if (index(uppercase(argi),'-FKERNEL ') .ne. 0) then
                 call get_command_argument(i+1,kernel_string)
              end if

              if (index(uppercase(argi),'-PREFIX1 ') .ne. 0) then
                 call get_command_argument(i+1,prefix1)
              end if

              if (index(uppercase(argi),'-PREFIX2 ') .ne. 0) then
                 call get_command_argument(i+1,prefix2)
              end if

           end do

           call removeallspaces(kernel_string)

           do i=1,ndim
              write(c,*) i
              var_names(i) = trim(adjustl(prefix1))//trim(adjustl(c))
           end do
           do i=1,ndim
              write(c,*) i
              var_names(ndim+i) = trim(adjustl(prefix2))//trim(adjustl(c))
           end do

           call checkvariables(var_names)

        case(SE_kernel)
           do i=1,argcount
              call get_command_argument(i,argi)

              if (index(uppercase(argi),'-DELTA ') .ne. 0) then
                 call get_command_argument(i+1,argi)
                 read(argi,*,iostat=ios) delta
                 if(ios .ne. 0) then
                    write(6,*) '>>> ERROR: DELTA("', trim(adjustl(argi)), '") has a wrong fromat!'
                    stop
                 end if
                 if(delta .le. 0) then
                    write(6,*) '>>> ERROR: DELTA(', delta, ') must be larger than 0!'
                    stop
                 end if
              end if

              if (index(uppercase(argi),'-THETA ') .ne. 0) then
                 do j=1,ndim
                    call get_command_argument(i+j,argj)
                    read(argj,*,iostat=ios) thetas(j)
                    if(ios .ne. 0) then
                       write(6,*) '>>> ERROR: ', j, 'th component of THETA("', trim(adjustl(argj)), '") has a wrong fromat!'
                       stop
                    end if
                    if(thetas(j) .le. 0) then
                       write(6,*) '>>> ERROR: ', j, 'th component of THETA(', thetas(j), ') must be larger than 0!'
                       stop
                    end if
                 end do
              end if

           end do
     end select
  end if

  do i=1,argcount
     call get_command_argument(i,argi)

     if (index(uppercase(argi),'-PERIODICITY ') .ne. 0) then
        do j=1,ndim
           call get_command_argument(i+j,argj)
           read(argj,*,iostat=ios) periodicities(j)
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: ', j, 'th component of PERIODICITY("', trim(adjustl(argj)), '") has a wrong fromat!'
              stop
           end if
        end do
     end if

     if (index(uppercase(argi),'-RANGE ') .ne. 0) then
        genrange=.true.
        do j=1,ndim
           call get_command_argument(i+2*(j-1)+1,argj)
           read(argj,*,iostat=ios) min_positions(j)
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: ', 2*(j-1)+1, 'th component of RANGE("', trim(adjustl(argj)), '") has a wrong fromat!'
              stop
           end if
           call get_command_argument(i+2*(j-1)+2,argj)
           read(argj,*,iostat=ios) max_positions(j)
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: ', 2*(j-1)+2, 'th component of RANGE("', trim(adjustl(argj)), '") has a wrong fromat!'
              stop
           end if
        end do
     end if

     if (index(uppercase(argi),'-FACTORIZATION ') .ne. 0) then
        call get_command_argument(i+1,argi)
        select case(trim(adjustl(lowercase(argi))))
           case('ch')
              factorization = CH
           case('qr')
              factorization = QR
           case('bk')
              factorization = BK
           case default
              write(6,*) '>>> ERROR: FACTORIZATION("', trim(adjustl(argi)), &
                          '") must be one of the followings: ch, qr, bk!'
              stop
        end select
     end if

     if (index(uppercase(argi),'-INTERVAL ') .ne. 0) then
        do j=1,ndim
           call get_command_argument(i+2*(j-1)+1,argj)
           read(argj,*,iostat=ios) interval(j,1)
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: ', 2*(j-1)+1, 'th component of INTERVAL("', trim(adjustl(argj)), '") has a wrong fromat!'
              stop
           end if
           call get_command_argument(i+2*(j-1)+2,argj)
           read(argj,*,iostat=ios) interval(j,2)
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: ', 2*(j-1)+2, 'th component of INTERVAL("', trim(adjustl(argj)), '") has a wrong fromat!'
              stop
           end if
        end do
     end if

     if (index(uppercase(argi),'-SPARSEMAXN ') .ne. 0) then
         call get_command_argument(i+1,argi)
         read(argi,*,iostat=ios) sparsemaxn
         if(ios .ne. 0) then
            write(6,*) '>>> ERROR: SPARSEMAXN("', trim(adjustl(argi)), '") has a wrong format!'
            stop
         end if
         if(sparsemaxn .le. 0) then
            write(6,*) '>>> ERROR: SPARSEMAXN(', sparsemaxn, ') must be larger than 0!'
            stop
         end if
     end if

     if (index(uppercase(argi),'-PARTIAL ') .ne. 0) then
        partial=.true.
        if(outval .or. outvaldev .or. outder .or. outderdev) then
           write(6,*) '>>> ERROR: PARTIAL cannot be combined with any of the followings: OUTVAL, OUTVALDEV, OUTDER, OUTDERDEV!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-JITTER ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) jitter
        if(ios .ne. 0) then
            write(6,*) '>>> ERROR: JITTER("', trim(adjustl(argi)), '") has a wrong format!'
            stop
         end if
        if(jitter .lt. 0) then
           write(6,*) '>>> ERROR: JITTER(', jitter, ') must be larger or equal to 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-PREDGRID ') .ne. 0) then
        do j=1,ndim
           call get_command_argument(i+j,argj)
           read(argj,*,iostat=ios) pred_grid(j)
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: ', j, 'th component of PREDGRID("', trim(adjustl(argj)), '") has a wrong fromat!'
              stop
           end if
           if(pred_grid(j) .le. 0) then
              write(6,*) '>>> ERROR: ', j, 'the component of PREDGRID(', pred_grid(j), ') must be larger than 0!'
              stop
           end if
        end do
     end if

     if (index(uppercase(argi),'-PREDPOSFILE ') .ne. 0) then
        predgriddata=.false.
        call get_command_argument(i+1,predposfile)
        inquire(file=predposfile, exist=file_exist)
        if(.not. file_exist) then
           write(6,*) '>>> ERROR: PREDPOSFILE("', trim(adjustl(predposfile)), '") does not exist!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-PREDROWNUM ') .ne. 0) then
        predrownum = .true.
     end if

     if (index(uppercase(argi),'-PREDRANGE ') .ne. 0) then
        predrange=.true.
        do j=1,ndim
           call get_command_argument(i+2*(j-1)+1,argj)
           read(argj,*,iostat=ios) pred_min_positions(j)
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: ', 2*(j-1)+1, 'th component of PREDRANGE("', trim(adjustl(argj)), '") has a wrong fromat!'
              stop
           end if
           call get_command_argument(i+2*(j-1)+2,argj)
           read(argj,*,iostat=ios) pred_max_positions(j)
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: ', 2*(j-1)+2, 'th component of PREDRANGE("', trim(adjustl(argj)), '") has a wrong fromat!'
              stop
           end if
        end do
     end if

     if (index(uppercase(argi),'-OPTSIGMAVAL ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) optsigmaval
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: OPTSIGMAVAL("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(optsigmaval .lt. 0) then
           write(6,*) '>>> ERROR: OPTSIGMAVAL(', optsigmaval, ') must be greater or equal to 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-OPTSIGMADER ') .ne. 0) then
        do j=1,ndim
           call get_command_argument(i+j,argj)
           read(argj,*,iostat=ios) optsigmader(j)
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: ', j, 'th component of OPTSIGMADER("', trim(adjustl(argj)), '") has a wrong fromat!'
              stop
           end if
           if(optsigmader(j) .lt. 0) then
              write(6,*) '>>> ERROR: ', j, 'th component of OPTSIGMADER(', optsigmader(j), ') must be greater or equal to 0!'
              stop
           end if
        end do
     end if

     select case(kernel_type)
        case(USER_kernel)

        case(SE_kernel)
           if (index(uppercase(argi),'-OPTDELTA ') .ne. 0) then
              call get_command_argument(i+1,argi)
              read(argi,*,iostat=ios) optdelta
              if(ios .ne. 0) then
                 write(6,*) '>>> ERROR: OPTDELTA("', trim(adjustl(argi)), '") has a wrong format!'
                 stop
              end if
              if(optdelta .lt. 0) then
                 write(6,*) '>>> ERROR: OPTDELTA(', optdelta, ') must be greater or equal to 0!'
                 stop
              end if
           end if

           if (index(uppercase(argi),'-OPTTHETA ') .ne. 0) then
              do j=1,ndim
                 call get_command_argument(i+j,argj)
                 read(argj,*,iostat=ios) opttheta(j)
                 if(ios .ne. 0) then
                    write(6,*) '>>> ERROR: ', j, 'th component of OPTTHETA("', trim(adjustl(argj)), '") has a wrong fromat!'
                    stop
                 end if
                 if(opttheta(j) .lt. 0) then
                    write(6,*) '>>> ERROR: ', j, 'th component of OPTTHETA(', opttheta(j), ') must be greater or equal to 0!'
                    stop
                 end if
              end do
           end if

     end select

     if (index(uppercase(argi),'-NUMDIFF ') .ne. 0) then
        call get_command_argument(i+1,argi)
        argi=trim(adjustl(lowercase(argi)))
        read(argi,*) numdiff
        if( (trim(adjustl(numdiff)) .ne. "one-sided") .and. &
            (trim(adjustl(numdiff)) .ne. "two-sided") ) then
           write(6,*) '>>> ERROR: NUMDIFF("', trim(adjustl(numdiff)), '") must be one of the followings: one-sided, two-sided!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-OPTIMIZER ') .ne. 0) then
        call get_command_argument(i+1,argi)
        argi=trim(adjustl(lowercase(argi)))
        read(argi,*) optimizer
        if( (trim(adjustl(optimizer)) .ne. "lbfgs") ) then
           write(6,*) '>>> ERROR: OPTIMIZER("', optimizer, '") must be one of the followings: lbfgs!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-NCORR ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) ncorr
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: NCORR("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(ncorr .le. 0) then
           write(6,*) '>>> ERROR: NCORR(', ncorr, ') must be larger than 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-MAXITER ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) maxiter
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: MAXITER("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(maxiter .le. 0) then
           write(6,*) '>>> ERROR: MAXITER(', maxiter, ') must be larger than 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-PRINTFREQ ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) printfreq
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: PRINTFREQ("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(printfreq .le. 0) then
           write(6,*) '>>> ERROR: PRINTFREQ(', printfreq, ') must be larger than 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-MAXDVAL ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) maxdval
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: MAXDVAL("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(maxdval .le. 0) then
           write(6,*) '>>> ERROR: MAXDVAL(', maxdval, ') must be larger than 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-MAXGRAD ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*) maxgrad
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: MAXGRAD("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(maxgrad .le. 0) then
           write(6,*) '>>> ERROR: MAXGRAD(', maxgrad, ') must be larger than 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-MAXRMSGRAD ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) maxrmsgrad
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: MAXRMSGRAD("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(maxrmsgrad .le. 0) then
           write(6,*) '>>> ERROR: MAXRMSGRAD(', maxrmsgrad, ') must be larger than 0!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-PATHSEARCH ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*) pathsearch
        if(pathsearch .ne. "neb") then
           write(6,*) '>>> ERROR: PATHSEARCH("', pathsearch, '") must be one of the followings: neb!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-NEBFORCE ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) nebforce
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: NEBFORCE("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
        if(nebforce .le. 0) then
           write(6,*) '>>> ERROR: NEBFORCE(', nebforce, ') must be larger than 0!'
           stop
        end if
     end if

     if (index(uppercase(argi), '-SEED ') .ne. 0) then
        call get_command_argument(i+1,argi)
        read(argi,*,iostat=ios) seed
        if(ios .ne. 0) then
           write(6,*) '>>> ERROR: SEED("', trim(adjustl(argi)), '") has a wrong format!'
           stop
        end if
     end if

     if (index(uppercase(argi),'-SHIFT ') .ne. 0) then
        shift=.true.
     end if

     if (index(uppercase(argi),'-GNUPLOT ') .ne. 0) then
        gnuplot=.true.
     end if

     if (index(uppercase(argi),'-PREDFILE ') .ne. 0) then
        call get_command_argument(i+1,predfile)
     end if

     if (index(uppercase(argi),'-GPOUTFILE ') .ne. 0) then
        call get_command_argument(i+1,gpoutfile)
     end if

     if (index(uppercase(argi),'-OPTFILE ') .ne. 0) then
        call get_command_argument(i+1,optfile)
     end if

     if (index(uppercase(argi),'-MINFILE ') .ne. 0) then
        call get_command_argument(i+1,minfile)
     end if
  end do

  if( sparsemaxn .gt. stot ) sparsemaxn = stot

  if ( (.not. ingp) .or. dochange ) then
     do i=1,argcount
        call get_command_argument(i,argi)

        if (index(uppercase(argi),'-SPARSIFICATION ') .ne. 0) then
           call get_command_argument(i+1,argi)
           argi=trim(adjustl(lowercase(argi)))
           read(argi,*) sparsification
           if( (trim(adjustl(sparsification)) .ne. 'none') .and. &
               (trim(adjustl(sparsification)) .ne. 'clustering') .and. &
               (trim(adjustl(sparsification)) .ne. 'grid') .and. &
               (trim(adjustl(sparsification)) .ne. 'file') ) then
               write(6,*) '>>> ERROR: SPARSIFICATION("', trim(adjustl(sparsification)), &
                          '") must be one of the followings: none, clustering, grid, file!'
               stop
           end if
        end if

        if (index(uppercase(argi),'-SPARSEMETHOD ') .ne. 0) then
           call get_command_argument(i+1,argi)
           argi=trim(adjustl(lowercase(argi)))
           read(argi,*) sparsemethod
           if( (trim(adjustl(sparsemethod)) .ne. 'dic') .and. &
               (trim(adjustl(sparsemethod)) .ne. 'dtc') .and. &
               (trim(adjustl(sparsemethod)) .ne. 'fitc') ) then
               write(6,*) '>>> ERROR: SPARSEMETHOD("', trim(adjustl(sparsemethod)), &
                          '") must be one of the followings: dic, dtc, fitc!'
               stop
           end if
        end if

        if (index(uppercase(argi),'-SPARSETYPE ') .ne. 0) then
           call get_command_argument(i+1,argi)
           argi=trim(adjustl(lowercase(argi)))
           read(argi,*) sparsetype
           if( (trim(adjustl(sparsetype)) .ne. 'val') .and. &
               (trim(adjustl(sparsetype)) .ne. 'der') .and. &
               (trim(adjustl(sparsetype)) .ne. 'all') ) then
               write(6,*) '>>> ERROR: SPARSETYPE("', trim(adjustl(sparsetype)), &
                          ' must be one of the followings: val, der, all!'
               stop
           end if
        end if

        if (index(uppercase(argi),'-CLUSTERING ') .ne. 0) then
           call get_command_argument(i+1,argi)
           argi=trim(adjustl(lowercase(argi)))
           read(argi,*) clustering
           if( (trim(adjustl(clustering)) .ne. 'every') .and. &
               (trim(adjustl(clustering)) .ne. 'random') .and. &
               (trim(adjustl(clustering)) .ne. 'kmeans') .and. &
               (trim(adjustl(clustering)) .ne. 'kmeans++') ) then
               write(6,*) '>>> ERROR: CLUSTERING("', trim(adjustl(clustering)), &
                          '") must be one of the followings: every, random, kmeans, kmeans++!'
               stop
           end if
        end if

        if (index(uppercase(argi),'-SPARSEPOINTS ') .ne. 0) then
           call get_command_argument(i+1,argi)
           read(argi,*,iostat=ios) sparsepoints
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: SPARSEPOINTS("', trim(adjustl(argi)), '") has a wrong format!'
              stop
           end if
           if(sparsepoints .le. 0) then
              write(6,*) '>>> ERROR: SPARSEPOINTS(', sparsepoints, ') must be larger than 0!'
              stop
           end if
        end if

        if (index(uppercase(argi),'-SPARSEFREQ ') .ne. 0) then
           call get_command_argument(i+1,argi)
           read(argi,*,iostat=ios) sparsefreq
           if(ios .ne. 0) then
              write(6,*) '>>> ERROR: SPARSEFREQ("', trim(adjustl(argi)), '") has a wrong format!'
              stop
           end if
           if(sparsefreq .le. 0) then
              write(6,*) '>>> ERROR: SPARSEFREQ(', sparsefreq, ') must be larger than 0!'
              stop
           end if
        end if

        if (index(uppercase(argi),'-SPARSEGRID ') .ne. 0) then
           do j=1,ndim
              call get_command_argument(i+j,argj)
              read(argj,*,iostat=ios) sparse_grid(j)
              if(ios .ne. 0) then
                 write(6,*) '>>> ERROR: ', j, 'th component of SPARSEGRID("', trim(adjustl(argj)), '") has a wrong format!'
                 stop
              end if
              if(sparse_grid(j) .le. 0) then
                 write(6,*) '>>> ERROR: ', j, 'th component of SPARSEGRID(', sparse_grid(j), ') must be larger than 0!'
                 stop
              end if
           end do
        end if

        if (index(uppercase(argi),'-SPARSEPOSFILE ') .ne. 0) then
           sparsegriddata=.false.
           call get_command_argument(i+1,sparseposfile)
           inquire(file=sparseposfile, exist=file_exist)
           if(.not. file_exist) then
              write(6,*) '>>> ERROR: SPARSEPOSFILE("', trim(adjustl(sparseposfile)), '") does not exist!'
              stop
           end if
        end if

        if (index(uppercase(argi),'-SPARSEROWNUM ') .ne. 0) then
           sparserownum=.true.
        end if

        if (index(uppercase(argi),'-SPARSERANGE ') .ne. 0) then
           sparserange=.true.
           do j=1,ndim
              call get_command_argument(i+2*(j-1)+1,argj)
              read(argj,*,iostat=ios) sparse_min_positions(j)
              if(ios .ne. 0) then
                 write(6,*) '>>> ERROR: ', 2*(j-1)+1, 'th component of SPARSERANGE("', trim(adjustl(argj)), '") has a wrong fromat!'
                 stop
              end if
              call get_command_argument(i+2*(j-1)+2,argj)
              read(argj,*,iostat=ios) sparse_max_positions(j)
              if(ios .ne. 0) then
                 write(6,*) '>>> ERROR: ', 2*(j-1)+2, 'th component of SPARSERANGE("', trim(adjustl(argj)), '") has a wrong fromat!'
                 stop
              end if
           end do
        end if

        if (index(uppercase(argi),'-KMEANSTHRESHOLD ') .ne. 0) then
           do j=1,ndim
              call get_command_argument(i+j,argj)
              read(argj,*,iostat=ios) kmeansthresholds(j)
              if(ios .ne. 0) then
                 write(6,*) '>>> ERROR: ', j, 'th component of KMEANSTHRESHOLD("', trim(adjustl(argj)), '") has a wrong format!'
                 stop
              end if
           end do
        end if

     end do
  end if

  select case (sparsification)
     case('cluster')
        if ( .not. (inval .or. inder) ) then
           write(6,*) '>>> ERROR: SPARSIFICATION was declared as cluster no input was given!'
           stop
        end if
     case('grid')
        if ( .not. (inval .or. inder .or. sparserange) ) then
           write(6,*) '>>> ERROR: SPARSIFICATION was declared as grid but neither input nor sparserange were given!'
           stop
        end if
     case('file')
        if (sparseposfile .eq. '' ) then
           write(6,*) '>>> ERROR: SPARSIFICATION was declared as file but no SPARSEPOSFILE was given!'
           stop
        end if
  end select

  if ( dochange ) then
     if (sparsification .eq. 'none') then
        write(6,*) '>>> ERROR: DOCHANGE was declared but no sparsification type was defined!'
        stop
     end if
  end if

  if ( outmin ) then
     do i=1,argcount
        call get_command_argument(i,argi)

        if (index(uppercase(argi),'-MINGRID ') .ne. 0) then
           do j=1,ndim
              call get_command_argument(i+j,argj)
              read(argj,*,iostat=ios) min_grid(j)
              if(ios .ne. 0) then
                 write(6,*) '>>> ERROR: ', j, 'th component of MINGRID("', trim(adjustl(argj)), '") has a wrong format!'
                 stop
              end if
              if(min_grid(j) .le. 0) then
                 write(6,*) '>>> ERROR: ', j, 'th component of MINGRID(', min_grid(j), ') must be larger than 0!'
                 stop
              end if
           end do
        end if

        if (index(uppercase(argi),'-MINPOSFILE ') .ne. 0) then
           mingriddata=.false.
           call get_command_argument(i+1,minposfile)
           inquire(file=minposfile, exist=file_exist)
           if(.not. file_exist) then
              write(6,*) '>>> ERROR: MINPOSFILE("', trim(adjustl(minposfile)), '") does not exist!'
              stop
           end if
        end if

        if (index(uppercase(argi),'-MINROWNUM ') .ne. 0) then
           minrownum=.true.
        end if

        if (index(uppercase(argi),'-MINRANGE ') .ne. 0) then
           minrange=.true.
           do j=1,ndim
              call get_command_argument(i+2*(j-1)+1,argj)
              read(argj,*,iostat=ios) min_min_positions(j)
              if(ios .ne. 0) then
                 write(6,*) '>>> ERROR: ', 2*(j-1)+1, 'th component of MINRANGE("', trim(adjustl(argj)), '") has a wrong fromat!'
                 stop
              end if
              call get_command_argument(i+2*(j-1)+2,argj)
              read(argj,*,iostat=ios) min_max_positions(j)
              if(ios .ne. 0) then
                 write(6,*) '>>> ERROR: ', 2*(j-1)+2, 'th component of MINRANGE("', trim(adjustl(argj)), '") has a wrong fromat!'
                 stop
              end if
           end do
        end if

     end do
  end if

  call random_seed(size=isize)
  allocate(iseed(isize), stat = alloc_failed)
  if( alloc_failed .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_read_args_rest] Unable to allocate memory for iseed!'
      stop
  end if
  do i=1, isize
     iseed(i) = seed + 37 * (i-1) 
  end do
  call random_seed(put=iseed)

  if(verbose) then
     write(6,*) 'GP calculation'
     write(6,*) 'Type of calculation:'
     if(ingp) then
        write(6,*) ' -> GP file(s) is (are) provided'
        if(dochange) then
           write(6,*) ' -> Change of sparse points positions is requested'
        end if
     end if
     if(inval) then
        if(inder) then
           write(6,*) ' -> Teaching from values and derivatives'
        else
           write(6,*) ' -> Teaching from values only'
        end if
     else if(inder) then
        write(6,*) ' -> Teaching from derivatives only'
     end if
     write(6,*) 'ndim = ', ndim
     write(6,*) 'ntot = ', ntot
     if(inval .or. inder .or. dochange) then
        write(6,*) 'sparsification = ', trim(sparsification)
        if(sparsification .eq. 'clustering') then
           write(6,*) 'clustering = ', trim(clustering)
        end if
        if(sparsification .ne. 'none') then
           write(6,*) 'sparsemethod = ', trim(sparsemethod)
           write(6,*) 'sparsetype = ', trim(sparsetype)
        end if
     end if
     write(6,*) 'seed = ', seed
  end if

end subroutine gp_general_read_args_rest

function uppercase(strin) result(strout)

  implicit none

  character(len=*) :: strin
  character(len=len(strin)) :: strout
  integer :: i, j
  character(len=26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  strout = strin

  do i=1, len_trim(strin)
     j = index(lower, strin(i:i))
     if(j>0) strout(i:i)=upper(j:j)
  end do

end function uppercase

function lowercase(strin) result(strout)

  implicit none

  character(len=*) :: strin
  character(len=len(strin)) :: strout
  integer :: i, j
  character(len=26), parameter :: lower = 'abcdefghijklmnopqrstuvwxyz'
  character(len=26), parameter :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  strout = strin

  do i=1, len_trim(strin)
     j = index(upper, strin(i:i))
     if(j>0) strout(i:i)=lower(j:j)
  end do

end function lowercase

!--------------------------------------------------------------------------------------------------

end module gp_general_read_mod
