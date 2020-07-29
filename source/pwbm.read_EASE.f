
      subroutine write_ann_subset(fname,igrid_fst,igrid_lst,dataval)
c  Write the ALT data as text file to save space
      implicit none

      integer   icell,igrid_fst,igrid_lst
      integer   nvals
      parameter (nvals=1206129)
      real*8    dataval(nvals)
      character fname*100,str2*20

      str2='"CellID" "ALT(m)"'
 124  format(i10,f10.2)

      open(302,file=fname)
      write(302,'(a20)') trim(str2)
      do icell=igrid_fst,igrid_lst
        write(302,124) icell,dataval(icell)
      enddo
      close(302)
      end

      subroutine write_ann_subset_v1(fname,igrid_fst,igrid_lst,ngrid,
     $   numdim,dataval)
       implicit none

       integer   nvals,numdim,ngrid
       integer   igrid_fst,igrid_lst
       integer   id,igrid,icell
       parameter (nvals=1206129)
       real*8    dataval(nvals,numdim)
       real      tempval(ngrid,numdim)
       character*100 fname

       open(unit=304,
     $     file=trim(fname),
     $     form='unformatted',
     $     access='direct',
     $     recl=ngrid*numdim*4,
     $     status='replace')

       do id=1, numdim
          do igrid=1,ngrid
              icell=igrid+igrid_fst-1
              tempval(igrid,id)=dataval(icell,id)
          enddo
       enddo

       write(304,rec=1) tempval
       close(304)
       end
 
      subroutine write_soil_subset(fname,igrid_fst,igrid_lst,ngrid,
     $    nlayer,numdim,dataval)
c  Subroutine writes the output file with binary format in float32 type
       implicit none

       integer nlayer,numdim,ngrid
       integer igrid_fst,igrid_lst
       integer id,ilayer,igrid,icell
       real*8  dataval(ngrid,numdim,nlayer)
       real    tempval(ngrid,numdim,nlayer)
       character*100 fname

       open(unit=303,
     $     file=trim(fname),
     $     form='unformatted',
     $     access='direct',
     $     recl=ngrid*nlayer*numdim*4,
     $     status='replace')

       do ilayer=1,nlayer
          do id=1, numdim
          do igrid=1,ngrid
             tempval(igrid,id,ilayer)=dataval(igrid,id,ilayer)
          enddo
          enddo
       enddo

       write(303,rec=1) tempval
       close(303)
       end


       subroutine write_ann_data_v1(fname,row_vals,col_vals,dataval)
       implicit none
       
       integer nvals,nrow,ncol,ndim
       integer i,j,id,irow,icol,icell
       parameter (ndim=46,nvals=1206129)
       parameter (nrow=1650,ncol=1776)
       integer   row_vals(nvals),col_vals(nvals)
       real*8    dataval(nvals,ndim)
       real    tempval(ncol,nrow)
       character*100 fname

       open(unit=201,
     $     file=trim(fname),
     $     form='unformatted',
     $     access='direct',
     $     recl=nrow*ncol*4,
     $     status='replace')

       do id=1,ndim
          do i=1,nrow
          do j=1,ncol
             tempval(j,i)=-999.0d0
          enddo
          enddo

          do icell=1,nvals
            irow=row_vals(icell)
            icol=col_vals(icell)
            tempval(icol,irow)=dataval(icell,id)
          enddo
         
          write(201,rec=id) tempval
       enddo
       close(201)
       end

      subroutine write_soil_data_v2(fname,row_vals,col_vals,iyr,nlayer,
     $ dataval)
       implicit none
       
       integer nvals,nrow,ncol,nlayer,length,mdate
       integer i,j,iindex,irow,icol,iyr,id,ilayer,icell
       parameter (nvals=1206129)
       parameter (nrow=1650,ncol=1776)
       integer   row_vals(nvals),col_vals(nvals)
       real*8    dataval(nvals,46,15)
       real      tempval(ncol,nrow)
       character*100 fname, fullname
       character*6   datestr

       length = index(fname, ' ') - 1
       do id=1, 46
          mdate = iyr*100+id
          write(datestr,'(i6)') mdate
          fullname = fname(1:length)//datestr//'.0-3m.flt32'

          open(unit=201,
     $       file=trim(fullname),
     $       form='unformatted',
     $       access='direct',
     $       recl=nrow*ncol*4,
     $       status='replace')

          do ilayer=1, nlayer
             do i=1,nrow
             do j=1,ncol
               tempval(j,i)=-999.0d0
             enddo
             enddo

             do icell=1,nvals
                irow=row_vals(icell)
                icol=col_vals(icell)
                tempval(icol,irow)=dataval(icell,id,ilayer)
             enddo
         
             write(201,rec=ilayer) tempval
          enddo
          close(201)
      enddo
      end


      subroutine cal_statis(input_array,annval,ndim,flag)
      implicit none
      real*8  input_array(ndim)
      real*8  annval, fsum
      integer num, id, flag, ndim
   
      fsum = 0.0d0
      num = 0
      do id=1,ndim
         if(input_array(id) .gt. -999.0d0) then
            fsum = fsum + input_array(id)
            num = num + 1
         endif
      enddo 
      if(num .eq. ndim) then
         if(flag .eq. 1) then 
            annval = fsum
         else
            annval = fsum/dble(num)
         endif
      else
         annval = -9999.0d0
      endif
    
      end           


       subroutine read_data_ann(fpath,dataname,iyr,
     $            row_vals,col_vals,dataval)
c  Subroutine reads the input file with binary format in float32 type
c  inputs:  fpath    - Character string designating path and file prefix name for input
c           dataname - the name of the dataset
c           row_vals - the row index for each grid
c           col_vals - the col index for each grid
c  Outputs: dataval  - Array of data values returned to main program 
       implicit none
 
       integer nvals,iyr
       integer nrow,ncol,ndim
       integer id,icell,irow,icol
       integer i,j,num
       integer length,length1
       parameter (nvals=1206129)
       parameter (nrow=1650)
       parameter (ncol=1776)
       parameter (ndim=46)
       integer row_vals(nvals),col_vals(nvals)
       real*8  dataval(46,nvals)
       real  tempval(ncol,nrow)
       character*100 fpath,dataname,fname
       character yearstr*4
   
       length = index(fpath,' ') - 1
       length1 = index(dataname, ' ') - 1
       write(yearstr,'(i4)') iyr
       fname = fpath(1:length)//dataname(1:length1)//
     $ yearstr//'.flt32'
       print*, fname

       open(unit=101,
     $     file=trim(fname),
     $     form='unformatted',
     $     access='direct',
     $     recl=4*nrow*ncol)
       
       do id=1, ndim
c          print*,id,' reading '
          read(101,rec=id) tempval
          num = 0
          do icell=1, nvals
             irow = row_vals(icell)
             icol = col_vals(icell)
             dataval(id,icell)=tempval(icol,irow)
             if(dataval(id,icell) .le. -999.0d0) then
               num = num + 1
             endif
          enddo
          if(num .gt. 0) then
             print*,num,':number of missing data' 
             stop
          endif
       enddo      
       close(101)
      
       end
