#ifdef _USE_NVTX_

MODULE mod_nvtx

  INTERFACE
     SUBROUTINE nvtxStartRange(label,color)
       BIND(C, name="nvtxStartRange")
       IMPLICIT NONE
       CHARACTER, INTENT(IN) :: label
       INTEGER, INTENT(IN) :: color
     END SUBROUTINE nvtxStartRange

     SUBROUTINE nvtxEndRange
       BIND(C, name="nvtxEndRange")
       IMPLICIT NONE
     END SUBROUTINE nvtxEndRange

     SUBROUTINE nvtxRangePush
       BIND(C, name="nvtxRangePush")
       IMPLICIT NONE
     END SUBROUTINE nvtxRangePush

     SUBROUTINE nvtxRangePop
       BIND(C, name="nvtxRangePop")
       IMPLICIT NONE
     END SUBROUTINE nvtxRangePop

  END INTERFACE

#endif /* _USE_NVTX_ */
