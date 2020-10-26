#ifdef _USE_NVTX_

! Fortran bindings for a small subset of the NVIDIA Tools Extensions library
module nvtx
  use iso_c_binding
  public :: nvtxrangepusha, nvtxrangepop
  public :: nvtxrangepushaargb
  interface

    ! Annotate the timeline with a message
    ! Parameters:
    ! * string : the message in a string format
    subroutine nvtxrangepusha(string) bind(C, name="nvtxRangePushA")
      use iso_c_binding, only : c_char
      character(kind=c_char) :: string(*)
    end subroutine nvtxrangepusha

    ! Annotate the timeline with both a message and an ARGB color
    ! Parameters:
    ! * string : the message in a string format
    ! * argb   : the color in argb format (example: Z'FF880000'
    subroutine nvtxrangepushaargb(string,argb) bind(C, name="_nvtxRangePushAARGB")
      use iso_c_binding, only : c_char, c_int
      character(kind=c_char) :: string(*)
      integer(kind=c_int), value  :: argb
    end subroutine nvtxrangepushaargb

    ! Place a mark on the timeline with a message
    ! Parameters:
    ! * string : the message in a string format
    ! NOT YET EXPOSED
    subroutine nvtxMarkA(string) bind(C, name="nvtxMarkA")
      use iso_c_binding, only : c_char
      character(kind=c_char) :: string(*)
    end subroutine

    ! Name an OS thread
    ! NOT YET EXPOSED
    subroutine nvtxNameOsThread(tid, string) bind(C, name="nvtxNameOsThread")
      use iso_c_binding, only : c_int, c_char
      integer(kind=c_int) :: tid
      character(kind=c_char) :: string(*)
    end subroutine
  end interface
end module nvtx

!-------------------------------------------------------------------------------

MODULE mod_nvtx

  USE ISO_C_BINDING
  USE nvtx
  IMPLICIT NONE

  PUBLIC :: nvtxStartRange, nvtxEndRange, nvtxRangePush, nvtxRangePop

  INTERFACE
     SUBROUTINE nvtxStartRange( label, color ) &
          BIND(C, name="nvtxStartRange")
       USE ISO_C_BINDING, ONLY: C_CHAR, C_INT
       IMPLICIT NONE
       CHARACTER(KIND=C_CHAR), INTENT(IN), VALUE :: label(*)
       INTEGER(KIND=C_INT), INTENT(IN), VALUE :: color
     END SUBROUTINE nvtxStartRange
  END INTERFACE

  INTERFACE
     SUBROUTINE nvtxEndRange( label ) &
          BIND(C, name="nvtxEndRange")
       USE ISO_C_BINDING, ONLY: C_CHAR, C_INT
       IMPLICIT NONE
       CHARACTER(KIND=C_CHAR), INTENT(IN), VALUE :: label(*)
     END SUBROUTINE nvtxEndRange
  END INTERFACE

  INTERFACE
     SUBROUTINE nvtxRangePush() &
          BIND(C, name="nvtxRangePush")
       IMPLICIT NONE
     END SUBROUTINE nvtxRangePush
  END INTERFACE

!  INTERFACE
!     SUBROUTINE nvtxRangePop() &
!          BIND(C, name="nvtxRangePop")
!       IMPLICIT NONE
!     END SUBROUTINE nvtxRangePop
!  END INTERFACE

  TYPE, BIND(C) :: nvtxEventAttributes
     integer(C_INT16_T) :: version=1
     integer(C_INT16_T) :: size=48 !
     integer(C_INT) :: category=0
     integer(C_INT) :: colorType=1 ! NVTX_COLOR_ARGB = 1
     integer(C_INT) :: color
     integer(C_INT) :: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
     integer(C_INT) :: reserved0
     integer(C_INT64_T) :: payload ! union uint,int,double
     integer(C_INT) :: messageType=1 ! NVTX_MESSAGE_TYPE_ASCII = 1 
     type(C_PTR) :: message  ! ascii char
  END TYPE nvtxEventAttributes

END MODULE mod_nvtx

#endif /* _USE_NVTX_ */
