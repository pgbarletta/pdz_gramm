      program min_cost
      implicit none

      integer, parameter :: ikind=selected_real_kind(p=6)! p= 6, 15 o 18 (nro de figuras)
      integer, parameter :: jkind=selected_real_kind(p=15)! p= 6, 15 o 18 (nro de figuras)
      character*20 infile_1, infile_2, outfile_1, mc, nc 
      integer i, j, k, t, mc, nc, flag, ierr
      double precision mods(max_mod, max_mod)

      call getarg (1, infile_1)!       matriz de modos de referencia
      call getarg (2, infile_2)!       matriz de modos de referencia 
      call getarg (3, mc)!         
      call getarg (4, nc)!         
      call getarg (5, outfile_1)!      

      open (11, file=infile_1)
      open (12, file=infile_2)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!! lee archivos !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      do i=1, tresn+6!                          leo el archivo de modos
        read(11,91) (modos2(i,j),j=1, pnumber)
      enddo
      do i=1, tresn+6!                          leo el archivo de modos
        read(12,91) (modosp2(i,j),j=1, pnumber2)
      enddo



* Esta parte me ordena los modos que por la   *
* perturabacion pudieron saltar en frecuencia *
      do i=1,ndat-6
         do j=1,ndat-6
            scpr(i,j)=0.0d0
            do k=1,ndat
               scpr(i,j)=scpr(i,j)+modos(k,i)*modosp(k,j)
            enddo
         enddo
      enddo
      do i=1,ndat-6
        do j=1,ndat-6
           ascpr(i,j)=int(scpr(i,j)**2*1.d5)
        enddo
      enddo
      do i=1,ndat-6
         do j=1,ndat-6
            if((j.lt.(i-4)).
     $or.(j.gt.(i+4))) then
              ascpr(i,j)=-1*ifix(sngl(1.d5))
            endif
         enddo
      enddo
      do i=1,ndat-6
         do j=1,ndat-6
            ascpr(i,j)=-1*ascpr(i,j)
         enddo
      enddo


      call apc(ndat-6,ascpr,iorden,z)



95      format(1i3,1x,1x,1x,1i6)
96      format(1i3,1x,1x,1x,F9.6)
    end
