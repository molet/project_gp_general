module gp_general_dat_mod

use prec_mod
use constants_mod
use gp_basic_dat_mod

implicit none

! variables for input reading
integer :: argcount
character(len=length) :: trainfile='training.txt'
character(len=length) :: predfile='prediction.txt'
character(len=length) :: gpoutfile='gp.txt'
character(len=length) :: optfile='opt.txt'
character(len=length) :: minfile='min.txt'
character(len=length), allocatable :: gpinfiles(:)
character(len=length) :: predposfile=''
character(len=length) :: sparseposfile=''
character(len=length) :: minposfile=''
character(len=length) :: gpformat='full'
integer :: kernel_type=1
character(len=length) :: kernel_string
character(len=length) :: prefix1='x'
character(len=length) :: prefix2='y'
character(len=length), allocatable :: var_names(:)
integer :: nline
integer :: lstart=1
integer :: lstop=0
integer :: linc=1
integer :: ntot
integer :: itot
integer :: stot=0
logical :: inval=.false.
logical :: inder=.false.
logical :: ingp=.false.
integer :: ngpinfiles=0
logical :: dochange=.false.
logical :: outval=.false.
logical :: outvaldev=.false.
logical :: outder=.false.
logical :: outderdev=.false.
logical :: outhess=.false.
logical :: outint=.false.
logical :: outintdev=.false.
logical :: outgp=.false.
logical :: outopt=.false.
logical :: outmin=.false.
logical :: outpath=.false.
logical :: commandsigmaval=.false.
logical :: commandsigmader=.false.
logical :: rownum=.false.
logical :: genrange=.false.
logical :: predrownum=.false.
logical :: predrange=.false.
logical :: predgriddata=.true.
logical :: sparserownum=.false.
logical :: sparserange=.false.
logical :: sparsegriddata=.true.
logical :: minrownum=.false.
logical :: minrange=.false.
logical :: mingriddata=.true.

logical :: shift=.false.
logical :: gnuplot=.false.
logical :: verbose=.false.
integer :: seed=3141592

! GP related variables
integer :: ndim=1
real(dp), allocatable :: positions(:,:)
real(dp), allocatable :: values(:), sigmaval(:)
real(dp), allocatable :: derivatives(:,:), sigmader(:,:)
real(dp), allocatable :: thetas(:)
real(dp), allocatable :: periodicities(:)
real(dp), allocatable :: interval(:,:)
integer, allocatable :: index_sparse(:)
real(dp), allocatable :: pseudo_sparse(:,:)
real(dp) :: aveval=0.0
real(dp) :: delta=1.0
real(dp) :: jitter=1.0e-14_dp
integer :: factorization=CH
logical :: partial=.false.
real(dp), allocatable :: min_positions(:), max_positions(:)
real(dp), allocatable :: range_positions(:)

! sparsification related variables
character(len=length) :: sparsification='none'
character(len=length) :: sparsemethod='dtc'
character(len=length) :: sparsetype='val'
integer :: sparsepoints=1000
character(len=length) :: clustering='every'
integer, allocatable :: sparse_grid(:)
real(dp), allocatable :: sparse_min_positions(:), sparse_max_positions(:)
real(dp), allocatable :: sparse_range_positions(:), sparse_widths(:)
integer :: sparsemaxn=100000
integer, allocatable :: sparse_iter(:)
integer :: sparsefreq=1
real(dp), allocatable :: kmeansthresholds(:)

! prediction related variables
integer, allocatable :: pred_grid(:)
real(dp), allocatable :: pred_min_positions(:), pred_max_positions(:)
real(dp), allocatable :: pred_range_positions(:), pred_widths(:)
integer :: pred_tot_grids
real(dp), allocatable :: pred_grid_positions(:,:)

! other variables
real(dp), allocatable :: val(:), valdev(:)
real(dp), allocatable :: der(:,:), derdev(:,:)
real(dp), allocatable :: hess(:,:,:)
real(dp), allocatable :: intg(:), intgdev(:)
real(dp) :: shiftminval, shiftminint
type(gp_basic) :: gp

! optimization related variables
character(len=length) :: optimizer="lbfgs"
integer :: maxiter=100
integer :: printfreq=1
real(dp) :: maxdval=0.000001_dp
real(dp) :: maxgrad=0.0001_dp
real(dp) :: maxrmsgrad=0.001_dp
integer :: ncorr=4

! hyper optimization related variables
character(len=length) :: numdiff="one-sided"
integer :: optdim = 0
real(dp) :: optdelta=0.0_dp
real(dp) :: optsigmaval=0.0_dp
real(dp), allocatable :: optsigmader(:)
real(dp), allocatable :: opttheta(:)
real(dp) :: orig_delta
real(dp), allocatable :: orig_sigmaval(:)
real(dp), allocatable :: orig_sigmader(:,:)
real(dp), allocatable :: orig_thetas(:)
integer :: optfound
real(dp), allocatable :: optpos(:)
real(dp) :: optfval
real(dp), allocatable :: optgrad(:)

! minimization related variables
integer, allocatable :: min_grid(:)
real(dp), allocatable :: min_min_positions(:), min_max_positions(:)
real(dp), allocatable :: min_range_positions(:), min_widths(:)
integer :: min_tot_grids
real(dp), allocatable :: min_grid_positions(:,:)
integer, allocatable :: minfound(:)
real(dp), allocatable :: minpos(:,:)
real(dp), allocatable :: minfval(:)
real(dp), allocatable :: mingrad(:,:)

! minimum energy pathway searching
character(len=length) :: pathsearch="neb"
real(dp) :: nebforce=100.0_dp

end module gp_general_dat_mod
