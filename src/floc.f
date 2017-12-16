c*********SUBROUTINE

      subroutine floc(lw, l1, submat, lev, nf, nprobe, nsub, nlev, rcor)
c	  input dimensions
      integer nf, nprobe, nsub, nlev
c     input
	  double precision lw(nsub), l1(nsub,nf), submat(nsub,nprobe), 
     *                   lev(nlev)

c     output
      double precision rcor(nprobe,nf)
c     variables
      integer i, j, ii, jj, k
      double precision qual(nsub), poicla(nlev), z(nlev)

c     initialize rcor	
      do 10, ii = 1, nprobe
         do 11 jj = 1, nf
            rcor(ii,jj)=0
11       continue
10    continue

      do 100, ii = 1, nprobe

c       initialize qual
        do 101, j = 1, nsub
           qual(j)=submat(j, ii)
101     continue


        do 110, jj = 1, nf
c          initialize and compute poicla		
           do 112, k = 1, nlev 
              poicla(k)=0
112	       continue

           do 120, k = 1, nlev
              do 121, i = 1, nsub
                 if(qual(i).eq.lev(k)) then
                    poicla(k)=lw(i)+poicla(k)
                 end if
121           continue
120        continue

c          initialize and compute z
           do 113, k = 1, nlev
              z(k)=0
113        continue

           do 130, k = 1, nlev
              do 131, i = 1, nsub
                 if(qual(i).eq.lev(k)) then
                    z(k)=(l1(i, jj)*lw(i)+z(k))
                 end if
131	          continue
130        continue

           do 140, k = 1, nlev
              if(poicla(k).gt.1e-12)  then 
                 z(k)=z(k)/poicla(k)
              end if
140        continue

           do 150, k = 1, nlev
              rcor(ii,jj)=rcor(ii,jj)+poicla(k) * z(k) * z(k)
150        continue

110        continue
100   continue
      return
      end
	
c	program main

c   dimensions	
c	integer nf, nprobe, nsub, nlev
		
c	input
c	real lw(285), l1(285,2), submat(285,830), lev(3)

c	output
c	real rcor(830,2)

c	variables
c	integer i,j

c	define dimensions
c	nf=2 
c	nprobe=830
c	nsub=285
c	nlev=3

c	read input
c	open(8, file='lw.txt')
c		read(8,*) lw
c	close(8)
	
c	open(9, file='l1.txt')
c		read(9,*) l1
c	close(9)
	
c	open(10, file='submat.txt')
c		read(10,*) submat
c	close(10)
	
c	lev(1)=-1
c	lev(2)=0
c	lev(3)=1
	
	
c	call foc(lw, l1, submat, lev, nf, nprobe, nsub, nlev, rcor)
	
c	do 99, i = 1, nprobe
c		print *, (rcor(i,j),j=1,nf)
c99	continue

c	stop
c	end