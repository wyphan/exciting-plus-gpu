      program main

      implicit none
      include 'silo.inc'

      integer ndim, err, dbid, opid, maxop 
      integer stat, i, j, n, lmname, lvname, itemp
      integer dims(3)
      real(4) dtemp
      real(4) orig(3), delta(3,3)
      real(4), allocatable :: coords(:,:)
      real(4), allocatable :: wf(:,:,:)
      character(256) fname, bname, mname, vname, pre, ext

      call getarg(1,fname)
      call readdxh(fname, dims, orig, delta)
      ndim = max(dims(1), dims(2), dims(3))
      allocate(coords(3,ndim))
      do i=1,3
        do j=1,dims(i)
          coords(i,j) = (j-1)*delta(i,i)+orig(i)
        enddo
      enddo
      allocate(wf(dims(1),dims(2),dims(3)))

      mname = 'quad_mesh'
      lmname = 9
      vname = 'wan_func'
      lvname = 8
      maxop = 1

      read(fname, '(3x,I3.3)')n
      pre = 'wf'
      ext = '.bin'
      write(bname, '(A2,I3.3,A4)')trim(pre),n,trim(ext)
      call readbin(bname, dims, wf)

      err = dbcreate('out.silo', 8, 0, DB_LOCAL, 'file info',
     &               9, DB_PDB, dbid)

      err = dbmkoptlist(maxop, opid)
c...Strangly, it seems that row major acually means column major
c...since wf is a fortran array it is stored in a column major
c...fashion but setting row major order yields the correct result
      err = dbaddiopt(opid, DBOPT_MAJORORDER, DB_ROWMAJOR)

      err = dbputqm(dbid, mname, lmname, 'x', 1, 'y', 1, 'z', 1,
     &              coords(1,:), coords(2,:), coords(3,:), dims, 3,
     &              DB_FLOAT, DB_COLLINEAR, opid, stat)

      err = dbputqv1(dbid, vname, lvname, mname, lmname, wf,
     &               dims, 3, DB_F77NULL, 0, DB_FLOAT, DB_NODECENT,
     &               opid, stat)

      err = dbclose(dbid)
      deallocate(wf,coords)

      end program
