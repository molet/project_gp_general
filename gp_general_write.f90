module gp_general_write_mod

use prec_mod
use verbosity_mod
use gp_basic_dat_mod
use gp_basic_mod
use gp_general_dat_mod

implicit none

contains

!--------------------------------------------------------------------------------------------------

subroutine gp_general_write_core

  implicit none

  ! single point calculation
  call gp_general_write_single_results(predfile)
  ! hyper optimization result
  call gp_general_write_hyper_optimum(optfile)
  ! gp file
  call gp_general_write_gp_file(gpoutfile)
  ! function minima
  call gp_general_write_function_minima(minfile) 

end subroutine gp_general_write_core

!--------------------------------------------------------------------------------------------------

subroutine gp_general_write_single_results(filename)

  implicit none

  character(len=*), intent(in) :: filename
  integer :: i, j, n
  integer :: ierr

  if( .not. (outval .or. outvaldev .or. outder .or. outderdev .or. outhess .or. outint .or. outintdev) ) return

  if(verbose) call verbosity_start_task('Printing results')
  ! do the printing
  open(unit = 7, file = filename, status = 'UNKNOWN', position = 'REWIND', form = 'FORMATTED', iostat = ierr)
  
  if( ierr .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_write_single_results] Unable to open file: ', filename
      stop
  end if

  do n=1,pred_tot_grids

     do i=1,ndim
        write(7,10,advance='NO') pred_grid_positions(i,n)
     end do

     if(outval) then
        write(7,20,advance='NO') val(n) - shiftminval
     end if

     if(outvaldev) then
        write(7,20,advance='NO') valdev(n)
     end if

     if(outder) then
        do i=1,ndim
           write(7,20,advance='NO') der(i,n)
        end do
     end if

     if(outderdev) then
        do i=1,ndim
           write(7,20,advance='NO') derdev(i,n)
        end do
     end if

     if(outhess) then
        do i=1,ndim
           do j=1,ndim
              write(7,20,advance='NO') hess(i,j,n)
           end do
        end do
     end if

     if(outint) then
        write(7,20,advance='NO') intg(n) - shiftminint
     end if

     if(outintdev) then
        write(7,20,advance='NO') intgdev(n)
     end if

     write(7,*)

     if(gnuplot .and. (ndim > 1)) then
        if(.not. predgriddata) then
           do i=1,ndim
              if( mod(n,pred_grid(i)) .eq. 0 ) then
                  write(7,*)
                  exit
              end if
           end do
        else
           if(any(pred_grid_positions(1:ndim-1,n) .ne. pred_grid_positions(1:ndim-1,n+1))) then
              write(7,*)
           end if
        end if
     end if

  end do

  close(7)

  if(verbose) call verbosity_finish_task('Printing results')

10 format(1X,F15.8)
20 format(1X,F15.8)

end subroutine gp_general_write_single_results

!--------------------------------------------------------------------------------------------------

subroutine gp_general_write_hyper_optimum(filename)

  implicit none

  character(len=*), intent(in) :: filename

  integer :: i
  integer :: ierr

  if ( .not. outopt) return

  open(unit = 9, file = filename, status = 'UNKNOWN', position = 'REWIND', form = 'FORMATTED', iostat = ierr)

  if( ierr .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_write_hyper_optimum] Unable to open file: ', filename
      stop
  end if

  do i=1,optdim
     write(9,10,advance='NO') 1.0
  end do

  write(9,20,advance='NO') '->'

  do i=1,optdim
     write(9,10,advance='NO') optpos(i)
  end do

  write(9,30,advance='NO') optfound

  write(9,40,advance='NO') optfval

  do i=1,optdim
     write(9,10,advance='NO') optgrad(i)
  end do

  write(9,*)

  close(9)

10 format(1X,F15.8)
20 format(1X,2A)
30 format(1X,I10)
40 format(1X,F15.8)

end subroutine gp_general_write_hyper_optimum

!--------------------------------------------------------------------------------------------------

subroutine gp_general_write_gp_file(filename)

  implicit none

  character(len=*), intent(in) :: filename

  if ( .not. outgp) return

  ! print out GP if requested
  if(verbose) call verbosity_start_task('Printing GP parameters')
  call gp_basic_write(self=gp,fileid=8,filename=filename,mode=gpformat)
  if(verbose) call verbosity_finish_task('Printing GP parameters')

end subroutine gp_general_write_gp_file

!--------------------------------------------------------------------------------------------------

subroutine gp_general_write_function_minima(filename)

  implicit none

  character(len=*), intent(in) :: filename

  integer :: i, n
  integer :: ierr

  if ( .not. outmin) return

  ! print out minima if requested
  if(verbose) call verbosity_start_task('Printing out found minima')

  open(unit = 9, file = filename, status = 'UNKNOWN', position = 'REWIND', form = 'FORMATTED', iostat = ierr)

  if( ierr .ne. 0 ) then
      write(6,*) '>>> ERROR: [gp_general_write_function_minima] Unable to open file: ', filename
      stop                                                                    
  end if 

  do n=1,min_tot_grids 

     do i=1,ndim
        write(9,10,advance='NO') min_grid_positions(i,n)
     end do

     write(9,20,advance='NO') '->'

     do i=1,ndim
        write(9,10,advance='NO') minpos(i,n)
     end do

     write(9,30,advance='NO') minfound(n)

     write(9,40,advance='NO') minfval(n)

     do i=1,ndim
        write(9,10,advance='NO') mingrad(i,n)
     end do

     write(9,*)

     if(gnuplot .and. (.not. mingriddata) .and. (ndim > 1)) then
        do i=1,ndim
           if( mod(n,min_grid(i)) .eq. 0 ) then
               write(9,*)
               exit
           end if
        end do
     end if

  end do

  close(9)
  
  if(verbose) call verbosity_finish_task('Printing out found minima')

10 format(1X,F15.8)
20 format(1X,2A)
30 format(1X,I10)
40 format(1X,F15.8)

end subroutine gp_general_write_function_minima

!--------------------------------------------------------------------------------------------------

end module gp_general_write_mod
