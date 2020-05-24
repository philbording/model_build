       Program mdlb
!
!   Model Building Program
!   R. P. Bording
!   Copy Right  May 2020
!
       parameter (ndpts = 4096)
       parameter (ndl = 4096)
       parameter (ndx=2000,ndz=2000)

       character*80 line

       character*80 model_name

       dimension ip1v(ndl),ip2v(ndl)
       dimension vpv(ndl),vsv(ndl),dev(ndl)
       dimension iptype(ndl)
!        = 1 single line segment  == line
!        = 2 multi-line segment   == segment

       dimension xzp(2,ndpts)
       dimension kxzp(ndpts)
       dimension vp (ndx,ndz)
       dimension vs (ndx,ndz)
       dimension den(ndx,ndz)

         
       open ( 10, file="velmod.dat",form="formatted",
     x           status="unknown")

       read(10,"(a)")  model_name
       write(6,*) "                 "
       write(6,"(a)")  model_name

       read(10,*) nmx,nmz
       read(10,*) xmin,zmin
       read(10,*) dx,dz
       read(10,*) vpmin,vsmin,denmin
       write(6,*) "                 "
       write(6,*) " Grid sizes nx,nz = ",nmx,nmz
       write(6,*) " Grid Spacing dx,dz ",dx,dz
       xmax = (nmx-1)*dx - xmin
       zmax = (nmz-1)*dz - zmin
       write(6,*) " Xmin, Zmin ",xmin,zmin
       write(6,*)" Min velocity, Vp, Vs, Density ",vpmin,vsmin,denmin
       write(6,*) "                 "
       close(10)
       open ( 11, file="velline.dat",form="formatted",
     x           status="unknown")

       ipc = 0
         write(6,*) " line inputs "
         write(6,*) "   Point 1, Point 2, Vp, Vs, Density  "
 99    continue
c      read( 11,*,end=100) ip1,ip2,vpl,vsl,denl
c      read( 11,*,end=100) ip1,ip2,vpl
       read( 11,"(a)",end=100) line
       write( 6,"(a)") line
       if( line(1:4) .eq. "Line" )  then 
         itype = 1
        go to 99
       endif
       if( line(1:4) .eq. "line" )  then 
         itype = 1
        go to 99
       endif
       if( line(1:4) .eq. "Segm" )  then 
         itype = 2
        go to 99
       endif
       if( line(1:4) .eq. "segm" )  then 
         itype = 2
        go to 99
       endif
       read(line,*) ip1,ip2,vpl
!      if( ip1 .eq. 0 .and. ip2 .eq. 0 ) go to 100
       vsl = 1300.0
       denl = 2.2
       if( mod(ipc,20) .eq. 0 ) then
         write(6,*) " line inputs "
         write(6,*) "   Point 1, Point 2, Vp, Vs, Density  "
       endif
       write( 6, *) ip1,ip2,vpl,vsl,denl
        ipc = ipc + 1
        iptype(ipc) = itype
        ip1v(ipc) = ip1
        ip2v(ipc) = ip2
        vpv(ipc) = vpl
        vsv(ipc) = vsl
        dev(ipc) = denl

        go to 99
 100    continue
       npc = ipc
       close(11)
       write(6,*) "                               "
!
! End of line segment input
!
       do i=1,ndpts
        kxzp(i) = 0
       enddo
!
       write(6,*) "                               "
       open ( 11, file="velpts.dat",form="formatted",
     x           status="unknown")
       write(6,*) "  Scale of input points     "
       read( 11, *)  xp_scale,zp_scale
       write( 6, *)  " x,z scale ",xp_scale,zp_scale
       write(6,*) "                               "
 98    continue
       read( 11,*,end=101)  ixzp,xp,zp
       kxzp(ixzp) = 1
       xzp(1,ixzp) = xp*xp_scale
       xzp(2,ixzp) = zp*zp_scale
       if( xzp(1,ixzp) .lt. xmin ) xzp(1,ixzp) = xmin
       if( xzp(1,ixzp) .gt. xmax ) xzp(1,ixzp) = xmax
       if( xzp(2,ixzp) .lt. zmin ) xzp(2,ixzp) = zmin
       if( xzp(2,ixzp) .gt. zmax ) xzp(2,ixzp) = zmax
       write( 6, *) " grid point number, x,z ", ixzp,xp,zp,
     x    xzp(1,ixzp),xzp(2,ixzp)
       go to 98
 101   continue
       close(11)
       write(6,*) "                               "
!
! end of grid point inputs, these points make up the line segments
!
       call setmin(vp,vs,den,vpmin,vsmin,denmin)

       write(6,*) " Arrays are set to zero "
!
       call bvelmodel(ip1v,ip2v,kxzp,iptype,xzp,
     x             nmx,nmz,
     x             vpv,vsv,denv,vp,vs,den,
     x             npc,dx,dz,xmin,zmin)
!
       write(6,*) " Grid sizes nx,nz = ",nmx,nmz

       open ( 20, file="vp.dat",form="unformatted",
     x           status="unknown")
          do ix=1,nmx
           write(20) (vp(ix,jz),jz=1,nmz)
          enddo
       close(20)
          do ix=1,nmx
           write(20) (vp(ix,jz),jz=1,nmz)
          enddo
          write(21) ((vs(ix,jz),jz=1,nmz),ix=1,nmx)
          write(22) ((den(ix,jz),jz=1,nmz),ix=1,nmx)
            end 


       subroutine setmin(vp,vs,den,vpmin,vsmin,denmin)

       parameter (ndx=2000,ndz=2000)

       dimension vp (ndx,ndz)
       dimension vs (ndx,ndz)
       dimension den(ndx,ndz)

c
c  python f2py calling sequence directives
c
c
c

cf2py  intent (in) vpmin
cf2py  intent (in) vsmin
cf2py  intent (in) denmin

cf2py  intent (out) vp
cf2py  intent (out) vs
cf2py  intent (out) den

c
c  end of directives
c
       do ix=1,ndx
        do jz=1,ndz
         vp(ix,jz) = vpmin
         vs(ix,jz) = vsmin
         den(ix,jz) = denmin
        enddo
       enddo
       return
       end

       subroutine bvelmodel(ip1v,ip2v,kxzp,iptype,xzp,
     x             nmx,nmz,
     x             vpv,vsv,denv,vp,vs,den,
     x             npc,dx,dz,xmin,zmin)

       parameter (ndpts = 4096)
       parameter (ndl = 4096)
       parameter (ndx=2000,ndz=2000)

       dimension ip1v(ndl),ip2v(ndl)
       dimension vpv(ndl),vsv(ndl),dev(ndl)
       dimension iptype(ndl)

!        = 1 single line segment  == line
!        = 2 multi-line segment   == segment

       dimension xzp(2,ndpts)
       dimension kxzp(ndpts)
       dimension vp (ndx,ndz)
       dimension vs (ndx,ndz)
       dimension den(ndx,ndz)
c
c  python f2py calling sequence directives
c
c
c

cf2py  intent (in) ip1v
cf2py  intent (in) 1p2v
cf2py  intent (in) kxpz
cf2py  intent (in) iptype 
cf2py  intent (in) xzp
cf2py  intent (in) nmx
cf2py  intent (in) nmz

cf2py  intent (in) vpv
cf2py  intent (in) vsv
cf2py  intent (in) denv

cf2py  intent (in) npc
cf2py  intent (in) dx
cf2py  intent (in) dz
cf2py  intent (in) xmin
cf2py  intent (in) zmin

cf2py  intent (in,out) vp
cf2py  intent (in,out) vs
cf2py  intent (in,out) den

c
c  end of directives
c

       write(6,*) " npc ",npc

       do ipc=1,npc

!
!     process each line segment
!
!   they must be sequenced from top to bottom
!
       ip1a = ip1v(ipc)
       ip2a = ip2v(ipc)
       write(6,*) " ip1a,ip2a ",ip1a,ip2a
       if( ip1a .gt. ip2a ) then
           ittp = ip1a
           ip1a  = ip2a
           ip2a = ittp
       endif
!
       if( kxzp(ip1a) .eq. 1 ) then
       if( kxzp(ip2a) .eq. 1 ) then
       
       if( iptype(ipc) .eq. 1 ) then
        kip = ip2a-ip1a
       endif

       if( iptype(ipc) .eq. 2 ) then
        kip = 1
       endif

         write(6,*) "              "
           write(6,*) ip1a,ip2a,kip
         write(6,*) "              "
       do ip1=ip1a,ip2a,kip
       if( ip1 .ne. ip2a ) then
       if( iptype(ipc) .eq. 1 ) then
        ip2 = ip1+kip
       endif
       if( iptype(ipc) .eq. 2 ) then
        ip2 = ip1 + 1
       endif
!        point is ok
       xa = xzp(1,ip1)
       za = xzp(2,ip1)
       xb = xzp(1,ip2)
       zb = xzp(2,ip2)
       write(6,*) xa,za,xb,zb
!
       if( xb .lt. xa ) then
         xt = xa
         xa = xb
         xb = xt
         zt = za
         za = zb
         zb = zt
       endif
!
       ixa = (xa-dx)/dx
       ixb = (xb+dx)/dx
       if( ixa .lt. 1 )  ixa = 1
       if( ixb .gt. nmx ) ixb = nmx
       do ix = ixa,ixb
        xxp = (ix-1)*dx + xmin
        if(  xxp .ge. xa ) then
         if( xxp .le. xb ) then
          amu = (xxp-xa)/(xb-xa)
          zzp = za + amu*(zb-za)

          jza = (zzp-2.0*dz)/dz
          if( jza .lt. 1 ) jza = 1
          jzb = nmz
          do jz = jza,jzb
           zpz = (jz-1)*dz + zmin
           if( zpz .ge. zzp ) then
            vp(ix,jz) = vpv(ipc)
            vs(ix,jz) = vsv(ipc)
            den(ix,jz) = dev(ipc)
           endif
          enddo
          endif
          endif
         enddo

!
           endif
          enddo
       else
       write(6,*) " error   p1 or p2 is bad ",ip1,ip2
       endif
       endif
       enddo
           




       return
       end
