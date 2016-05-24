program gp_general

use gp_general_read_mod
use gp_general_main_mod
use gp_general_write_mod

implicit none

! Read arguments and allocate variables
call gp_general_read_core

! Do the job
call gp_general_main_core

! Print out results
call gp_general_write_core

end program gp_general
