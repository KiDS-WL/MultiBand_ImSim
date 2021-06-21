      implicit none

C Fortran code to apply a post-measurement catalogue correction for
C weight anisotropy
C LM April-May 2014

      integer j, nbin,sbin,i,n,bin(1000000),k,status,ngal,ich_len,
     > iargc, writeout
      real*8 data(23),num(25),sumxy(25),sumxx(25),smin,smax,
     > maxmoment,efuncmax, gvar, gmom, pmoment, theta, pi,
     > gcvar,slope(25),gtheta,deltatheta,gvarcorr,gwcorr,
     > xmean(25),sumx(25),sumy(25),e1,e2,sumw,e1c,e2c,sumc,
     > delta1, delta2, qval, dmod,parval(7),l9,l17,lgvar
      character*128 filename,param,string,outname

      if (iargc() .ne. 2) then
         write (*,'(2a)') 
     >    'isoweight <lensfit ascii table name> ',
     >    '<output ascii table name>'
         call exit(0)
      endif

      pi = 3.14159265

C get the prior parameters
      call getenv('PRIOR_PARAMETERS', param)
      if (param.eq."") then
         write (*,'(a)') 
     > ' error: PRIOR_PARAMETERS environment variable not specified'
         call exit(2)
      else
         write (*,'(2a)') 'PRIOR PARAMETERS ',param
      endif

      open (1,file=param,status='old',iostat=status)
      if (status.ne.0) then
         write (*,'(2a)') ' could not open PRIOR_PARAMETERS file ',param
         call exit(2)
      endif
      n=0
      status=0
      do while (status.eq.0 .and. n.lt.7)
         read(1,'(a)',iostat=status) string
         if (status.eq.0 .and. string(1:1).ne.'/') then
            n=n+1
            read(string,*) parval(n)
         endif
      enddo
      close (1)

      efuncmax = parval(4)
      write (*,*) ' prior efunc_emax = ',efuncmax
      if (efuncmax.lt.0.7 .or. efuncmax.gt.0.9) call exit(2)

      call getarg(1,filename)
      open (1,file=filename,status='old')

      call getarg(2,outname)
      open (3,file=outname,status='unknown')

      do j = 1,24
         read (1,'(a)') string
         write (3,'(a)') string(:ich_len(string))
      enddo

      nbin = 10

      do j = 1,nbin
         num(j)=0.
         sumxy(j) = 0.
         sumxx(j) = 0.
         sumy(j) = 0.
         sumx(j) = 0.
      enddo

      maxmoment = efuncmax*efuncmax/2.

C  size bin limits
      smin = log(1.)
      smax = log(6.)

C  snr bin limits
C      smin = log(8.)
C      smax = log(200.)

      ngal = 0
      status = 0
      writeout = 0
      j = 0
      do while (status.eq.0)
         read (1,*,iostat=status) data
         if (status.eq.0) then
            j = j + 1
            bin(j) = 0
            ngal = ngal + 1

C  size binning
            if (data(7).gt.0.) then
               sbin = 1 + int(nbin*(log(data(7))-smin)/(smax-smin))
               if (sbin .gt. nbin) sbin = nbin
               if (sbin .lt. 1) sbin = 1
               bin(j) = sbin
            endif

C  snr binning
C            if (data(11).gt.0.) then
C               sbin = 1 + int(nbin*(log(data(11))-smin)/(smax-smin))
C               if (sbin .gt. nbin) sbin = nbin
C               if (sbin .lt. 1) sbin = 1
C               bin(j) = sbin
C            endif


            if (data(5).gt.0. .and. data(17).gt.0. 
     >           .and. nint(data(6)).eq.0
     >           .and. data(9).gt.0. .and. data(7).gt.0.) then
               
C  extract offsets needed to convert variance into weight
               gvar = data(17)
               gmom = 2./data(5)
               gcvar = gvar*maxmoment/(maxmoment-gvar)
               pmoment = gmom - gcvar
               if (writeout.eq.0) then
                  write (*,*) ' prior sigma ',sqrt(pmoment/2.)
                  writeout = 1
               endif

C normalise the measurement variance by flux squared
               gvar = gvar*data(9)**2
               if (gvar.gt.0.) then
                  gvar = log(gvar)
               else
                  write (*,*)' gvar ',data(9),gvar
                  write (*,*) data
               endif
               if (gvar.gt.log(20.)) gvar = log(20.)

               l9 = -99.
               if (data(9).gt.0.) l9 = log10(data(9))
               l17 = -99.
               if (data(17).gt.0.) l17 = log10(data(17))
               lgvar = -99.
C               write (9,*) data(9),gvar,data(17),l9,l17,
C     >              log10(data(7))
               

               
C  psf orientation
               theta = atan2(data(13),data(12))
C  galaxy orientation
               gtheta = atan2(data(4),data(3))
C  mod difference
               deltatheta = cos(gtheta-theta)
            
C  correlate flux-normalised gvar against delta-theta in bins of size

               sumxy(sbin) = sumxy(sbin) + gvar*deltatheta
               sumxx(sbin) = sumxx(sbin) + deltatheta*deltatheta
               sumy(sbin) = sumy(sbin) + gvar
               sumx(sbin) = sumx(sbin) + deltatheta
               num(sbin) = num(sbin) + 1.
               
            endif
         endif
      enddo

      close(1)

C  aggregate bins with small numbers
      do j = 1,nbin
         i=j+1
         do while (num(j).lt.500 .and. i.le.nbin)
            num(j) = num(j) + num(i)
            sumxx(j) = sumxx(j) + sumxx(i)
            sumxy(j) = sumxy(j) + sumxy(i)
            sumx(j) = sumx(j) + sumx(i)
            sumy(j) = sumy(j) + sumy(i)
            sumx(i) =0.
            sumy(i) =0.
            sumxy(i) =0.
            sumxx(i) =0.
            num(i) = 0.
            i = i + 1
            do k = 1,ngal
               if (bin(k).eq.i) bin(k)=j
            enddo
         enddo
      enddo

      n = 0
      do j = 1,nbin
         if (num(j).gt.0.) n= n + 1
      enddo

      if (num(n).lt.500) then
         j = n-1
         i = n
            num(j) = num(j) + num(i)
            sumxx(j) = sumxx(j) + sumxx(i)
            sumxy(j) = sumxy(j) + sumxy(i)
            sumx(j) = sumx(j) + sumx(i)
            sumy(j) = sumy(j) + sumy(i)
            sumx(i) =0.
            sumy(i) =0.
            sumxy(i) =0.
            sumxx(i) =0.
            num(i) = 0.
          n = n - 1
            do k = 1,ngal
               if (bin(k).eq.i) bin(k)=j
            enddo
       endif


      nbin = n

      write (*,'(a)') 'bin number mean slope '
      do j = 1,nbin
         slope(j) = 0.
         if (num(j) .gt. 10.) then
            slope(j) = (sumxy(j)-sumx(j)*sumy(j)/num(j))/
     >           (sumxx(j)-sumx(j)*sumx(j)/num(j))
            xmean(j) = sumx(j)/num(j)
          endif
         write (*,*) j,num(j),xmean(j),slope(j)
      enddo

C  reopen file and worked out corrected values

      open (1,file=filename,status='old')

      do j = 1,24
         read (1,*)
      enddo

      e1=0.
      e2=0.
      sumw=0.
      e1c=0.
      e2c=0.
      sumc=0.
      n=0

      do j = 1,ngal
         read (1,*) data
         if (data(5).gt.0. .and. data(17).gt.0. 
C     >        .and. nint(data(6)).eq.0
     >        .and. data(9).gt.0. .and. data(7).gt.0.) then

C  extract offsets needed to convert variance into weight
            gvar = data(17)
            gmom = 2./data(5)
            if (maxmoment .gt. gvar) then

               gcvar = gvar*maxmoment/(maxmoment-gvar)
               pmoment = gmom - gcvar

               sbin = bin(j)
               if (sbin .le. 0) then
                  write (*,*) ' error bin = ',sbin
                  call exit(2)
               endif

C  psf orientation
               theta = atan2(data(13),data(12))
C     galaxy orientation
               gtheta = atan2(data(4),data(3))
C  mod difference
               deltatheta = cos(gtheta-theta)

C  correction
               gvarcorr = gvar*(1. - 
     >              slope(sbin)*(deltatheta-xmean(sbin)))

C  limit correction to 5 percent
               if (gvarcorr < 0.95*gvar) gvarcorr = 0.95*gvar
               if (gvarcorr > 1.05*gvar) gvarcorr = 1.05*gvar

C  new moment
               if (maxmoment > gvarcorr) then
                  gwcorr = 2./(pmoment + 
     >                 gvarcorr*maxmoment/(maxmoment-gvarcorr))
               else
                  gwcorr = 0.
               endif

               if (nint(data(6)).eq.0) then
                  e1c = e1c + gwcorr*data(3)
                  e2c = e2c + gwcorr*data(4)
                  e1 = e1 + data(5)*data(3)
                  e2 = e2 + data(5)*data(4)
                  sumw = sumw + data(5)
                  sumc = sumc + gwcorr
                  n = n + 1
               endif

C               write (*,*) sbin,slope(sbin),xmean(sbin),deltatheta,
C     >              data(5),gwcorr,data(17),gvarcorr

               data(5) = gwcorr
               data(17) = gvarcorr

            endif
         endif

         write (3,
     > '(f10.6,f11.6,3f8.4,i3,2f10.4,2f9.1,f12.6,9f8.4,f6.2,i4,i12)')
     >        (data(i),i=1,5),nint(data(6)),(data(i),i=7,21),
     >        nint(data(22)),nint(data(23))

      enddo

      write (*,*) n, ' zero-fitclass galaxies used with mean weight ',
     >     sumw/n
      write (*,'(a,2f10.6)') 
     > ' uncorrected weighted mean shear values ',e1/sumw,e2/sumw
      write (*,'(a,2f10.6)') 
     > ' corrected weighted mean shear values ',e1c/sumc,e2c/sumc

      end


C+
      INTEGER FUNCTION ICH_LEN(STRING)
C
C     I C H _ L E N
C
C     Returns the 'logical length' of a character string.
C
C     Parameter -      (">" input, "<" output)
C
C     (>) STRING    (Character) The string whose length is required.
C
C     Returns -
C
C     (<) ICH_LEN   (Integer) The number - starting from 1 - of the
C                   last non-blank character in STRING.  Under some
C                   circumstances this is the logical length of the
C                   string.  If the string is completely blank, ICH_LEN
C                   returns 0.
C
C                                         KS / UCL  8th June 1982
C+
      CHARACTER*(*) STRING
C
      ICH_LEN=0
      DO I=1,LEN(STRING)
         IF (STRING(I:I).NE.' ') THEN
            ICH_LEN=I
         END IF
      END DO
C
      END
C+

