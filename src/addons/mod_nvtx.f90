#ifdef _USE_NVTX_

!-------------------------------------------------------------------------------
! Fortran bindings for a small subset of the NVIDIA Tools Extensions library
MODULE nvtx

  USE ISO_C_BINDING

  PUBLIC :: nvtxEventAttributes
  PUBLIC :: nvtxRangePushA, nvtxRangePushAArgb, nvtxMarkA, nvtxRangePushA
  PUBLIC :: nvtxRangePushEx, nvtxRangePop
  PUBLIC :: nvtxStartRange, nvtxEndRange

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
     type(C_PTR) :: message ! ascii char
  END TYPE nvtxEventAttributes

  PRIVATE TYPE(nvtxEventAttributes) :: event

  INTERFACE
    ! Annotate the timeline with a message
    ! Parameters:
    ! * string : the message in a string format
    subroutine nvtxRangePushA(string) bind(C, name="nvtxRangePushA")
      use iso_c_binding, only: c_char
      character(kind=c_char) :: string(*)
    end subroutine nvtxrangepusha
  END INTERFACE

  INTERFACE
    ! Annotate the timeline with both a message and an ARGB color
    ! Parameters:
    ! * string : the message in a string format
    ! * argb   : the color in argb format (example: Z'FF880000')
    subroutine nvtxRangePushAArgb(string,argb) bind(C, name="nvtxRangePushAARGB")
      use iso_c_binding, only : c_char, c_int
      character(kind=c_char) :: string(*)
      integer(kind=c_int), value  :: argb
    end subroutine nvtxRangePushAArgb
  END INTERFACE

  INTERFACE
    ! Place a mark on the timeline with a message
    ! Parameters:
    ! * string : the message in a string format
    subroutine nvtxMarkA(string) bind(C, name="nvtxMarkA")
      use iso_c_binding, only : c_char
      character(kind=c_char) :: string(*)
    end subroutine
  END INTERFACE

  INTERFACE
    ! Name an OS thread
    ! Parameters:
    ! * tid : the thread ID
    ! * string : the message in a string format
    subroutine nvtxNameOsThread(tid, string) bind(C, name="nvtxNameOsThread")
      use iso_c_binding, only : c_int, c_char
      integer(kind=c_int) :: tid
      character(kind=c_char) :: string(*)
    end subroutine
  END INTERFACE

  INTERFACE
     ! Push range with custom label and standard color
     ! Parameters:
     ! * string : the message in a string format
     ! * color : the color in argb format (example: Z'FF880000')
     SUBROUTINE nvtxRangePushA( label, color ) BIND(C, name="nvtxRangePushA")
       USE ISO_C_BINDING, ONLY: C_CHAR, C_INT
       IMPLICIT NONE
       CHARACTER(KIND=C_CHAR), INTENT(IN), VALUE :: label(*)
       INTEGER(KIND=C_INT), INTENT(IN), VALUE :: color
     END SUBROUTINE nvtxRangePushA
  END INTERFACE

  INTERFACE
     ! Push range with custom label and custom color
     ! Parameters:
     ! * event : the trace event
     SUBROUTINE nvtxRangePushEx( event ) BIND(C, name="nvtxRangePushEx")
       USE ISO_C_BINDING, ONLY: C_CHAR, C_INT
       USE nvtx, IMPORT: nvtxEventAttributes
       IMPLICIT NONE
       TYPE(nvtxEventAttributes) :: event
     END SUBROUTINE nvtxRangePushEx
  END INTERFACE

  INTERFACE
     ! Pop range
     SUBROUTINE nvtxRangePop() BIND(C, name="nvtxRangePop")
       IMPLICIT NONE
     END SUBROUTINE nvtxRangePop
  END INTERFACE

CONTAINS

  ! Start range with a predefined color id
  subroutine nvtxStartRange( name, id )
    character(kind=c_char,len=*) :: name
    integer, optional :: id
    type(nvtxEventAttributes) :: event

    tempName=trim(name)//c_null_char

    if ( .not. present(id) ) then
       call nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call nvtxRangePushEx(event)
    end if
  end subroutine nvtxStartRange

  ! End range
  subroutine nvtxEndRange
    call nvtxRangePop
  end subroutine nvtxEndRange

END MODULE nvtx

#endif /* _USE_NVTX_ */
