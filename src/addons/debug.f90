subroutine dbg_open_file
use modmain
implicit none
open(fdbgout,file=trim(adjustl(fdbgname)),form="FORMATTED",status="OLD",&
  position="APPEND")
return
end

subroutine dbg_close_file
use modmain
implicit none
close(fdbgout)
return
end

!==============================================================================
! Debugging subroutine
! To show filename and line number, define the preprocessor macro
!
! #define ASSERT(isok,msg,ival) \
!   CALL assert( isok, msg, ival, sub, __FILE__, __LINE__ )
!
! then use the macro to call this subroutine

SUBROUTINE assert( lisok, msg, ival, sub, fname, line )
  IMPLICIT NONE

  ! Arguments
  LOGICAL, INTENT(IN) :: lisok
  CHARACTER(LEN=*), INTENT(IN) :: msg
  INTEGER, INTENT(IN) :: ival
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: sub, fname
  INTEGER, INTENT(IN), OPTIONAL :: line

  IF( .NOT. lisok ) THEN

     IF( PRESENT(sub) ) THEN
        WRITE(*,*) 'Assert(', sub, '): ', msg, ival
     ELSE
        WRITE(*,*) msg, ival
     END IF ! sub

     IF( PRESENT(fname) .AND. PRESENT(line) ) THEN
        WRITE(*,*) TRIM(fname), ' line ', line
     END IF ! fname, line

     STOP '** ASSERTION ERROR **'

  END IF ! lisok
  RETURN
END SUBROUTINE assert
