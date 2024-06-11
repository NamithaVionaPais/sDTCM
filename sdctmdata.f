****************************************************************
*                            sDCTM                              *
*                                                               *
*****************************************************************

*****************************************************************
*nm : Number of documents					*
*nc : Number of classes associated to the response		*
*nt : Length of the time series					*
*nv : length of the vocabulary					*
*nj : Number of latent topics					*
*nem: Number of VEM iterations					*
*****************************************************************

*****************************************************************
*OBSERVED DATA							*
*y(nm,nc): Response variable					*
*w(nm,nt,nv): word time series  				*
*****************************************************************

*****************************************************************
*LATENT VARIABLES 				             	*
*delta(nm,nt,nj): latent topic structure 	                *
*Z(nm,nt,nj): latent topic indicator				*
*****************************************************************

*****************************************************************
*MODEL PARAMETERS						*
*phi(nj,nj)  							*
*sigma(nj,nj)						        *
*beta(nj,nv): matrix of word probabilities			*
*eta(nc,nj) vector of regression parameters			*
*****************************************************************

*****************************************************************
*VARIATIONAL PARAMETERS	                  			*
*dhat(nm,nt,nj)							*
*s2h (1,1)						        *
*gamma(nm,nt,nj)			                        *
*****************************************************************



*****************************************************************
*                Generative Process                             *
*****************************************************************
*  For each word w(id,it,:)     				*
*  Choose delta[t,:] ~ N_nj(phi%*%delta[t-1,:],sigma)           *
*  Choose Z[t,:] ~ MULT(1, f(delta[t,:]))                       *
*  f (.) is the Multinomial logit transformation                *
*  Choose word w(id,it,:) ~ MULT(1,beta_Z[t,:])                 *
*								*
*  Draw class label y(id,:) ~softmax(zbar,eta)                  *
*  zbar(1,nj) = mean(Z[t,:])                                    *
****************************************************************


*****************************************************************
*       routines called and purpose:                            *
*****************************************************************
* For E step						        *
* 						                *
* For M step						        *
* fbeta: update equation for beta 				*
* fphi : function to optimize to update phi			*
* fsig : function to optimize to update sigma			*
* feta : function to optimize to update eta			*						                
*****************************************************************



*****************************************************************
*****************************************************************
*                         MAIN PROGRAM                          *
*****************************************************************
*****************************************************************

	implicit real*8 (a-h,o-z)
	parameter(nj=3,nv=4,nm=86,nt=57,nc=2,nem=100,nr=500)
	real*8 tsigma(nj,nj),tphi(nj,nj),tbeta(nj,nv),teta(nc,nj)
	real*8 shp,sum,tmpgam(nv),tempr(nj)
	real*8 delprev(nj,1),delta(nm,nt,nj)
	real*8 delvec(nj,1),delvec1(1,nj)
	real*8 temp(nj,1),sumdelj,pidt(nm,nt,nj)
	real*8 izdt(nm,nt,nj),sumz(nj),zbar(nj,1)
	real*8 wdt(nm,nt,nv), sump, matp(nc,1)
	real*8 rank fSigma(nj,nj),y(nm,nc)
	real   pij(nj),pij1(nv), pijy(nc),ysum(nc)
	integer irmult(1,nj),irmult1(1,nv),irmultc(1,nc)
	integer ind1
	integer ntem,ns


	real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
	real*8 sgamma,gamma(nm,nt,nj)
	real*8 res1(1),s2h, acc1
	real*8 dhat(nm,nt,nj),dhatd(nt,nj),dhatd1(nt,nj)
	real*8 tempn(1)
	real*8 gammad(nt,nj)
	real*8 gammaup(nm,nt,nj)
	real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
	real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
	real*8 fmd(nt,nj,1),fvd(nt,nj,nj)
	real*8 smd(nt,nj,1),svd(nt,nj,nj)
	real*8 elbovec(nem)
	
	integer nmntnj,njnj,ncnj,nmnt,ntnj,nk
	real*8 dhatvec(nm*nt*nj),gammavec(nm*nt*nj)
	real*8 betau(nj,nv), gamu(nm,nt,nj)
	real*8 phivec(nj*nj),sigvec(nj+(nj*(nj-1)/2)),etavec(nc*nj)
	real*8 fdhr,fs2,fg,fp,fs, fe,elbo,elbon,elbop
	real*8 fp1, fs11

	real*8 dhatgvec(nm*nt*nj),xsc1(nm*nt*nj),g1(nm*nt*nj)
	real*8 dhatdvec(nt*nj),xsc11(nt*nj),dhatdgvec(nt*nj)
	real*8 dlb(nt*nj),dub(nt*nj)
	real*8 g11(nt*nj)
	real*8 s2hg,s2hlb,s2hub, plb(nj*nj),pub(nj*nj)
	real*8 phigvec(nj*nj),xsc2(nj*nj),g2(nj*nj)
	real*8 sigmagvec(nj+(nj*(nj-1)/2))
	real*8 etagvec(nc*nj),xsc3(nc*nj),g3(nc*nj)
	real*8 eub(nc*nj),elb(nc*nj)
	real*8 slb(nj*nj),sub(nj*nj)
	real*8 fdhatval,fs2hval,fgammaval
	real*8 fphival,fsigval,fetaval
	real*8 nact1,iact1((2*nm*nt*nj)+(nm*nt))
	real*8 alamda1(nm*nt*nj)
	real*8 ip(7),rp(7)
	real*8 ip1(7),rp1(7)
	real*8 phiout(nj,nj), fout
	real*8 sigout(nj,nj),fsout
	real*8 etaout(nc,nj),feout
	real*8 dhatdout(nt,nj),fdout
	real*8 p1
	integer ytr(nm)
	CHARACTER*30  prdata(nm,nt)

	real*8 treta(nj,nc),yp(nm,nc)
	real*8 gsum,matr(nc,1),gbar(nj,1)
	real*8 matr1(nc)

	real*8 iden(nj,nj)
	common /matrices/iden
	
	real*8 ftol,fvalue
        integer ibtype,maxfcn	
        CHARACTER*30 INPUT1,INPUT2

        character*40 outputfile1,outputfile2
        external DRNNOR,DRNGAM,RNSET
	external DCHFAC,DRNMVN,RNMTN
	external DBCPOL,DUMCGF,DLCONF
	external DEVLRG, DUMINF
	external DEVLSF
	external fdhat,fs2h,fgamma,fdhatd
	external fbeta,fphi,fsig,feta
	

c	 common/alldat/dat,ldat
c	 common/epsmom/mue,vare,skewe,kurte

	 common/vpar/gamma,dhat,s2h
	 common/mpar/phi,sigma,beta,eta
	 common/dat/wdt,y
         common/gd/gammad
 	 common/fphi1/fp1
	 common/fsig1/fs11

c OUTPUT
        outputfile1 ='sdctmd.out'
        OPEN (UNIT=6,FILE=outputfile1,ACCESS='SEQUENTIAL',
     &      STATUS='UNKNOWN')

        outputfile2 ='datasdctm.out'
        OPEN (UNIT=7,FILE=outputfile2,ACCESS='SEQUENTIAL',
     &      STATUS='UNKNOWN')


c set seed
	iseed=221196
	call RNSET(iseed)

c identity matrix
	call MAKEI(nj,iden)


*****************************************************************
*	       INPUT DATA					*
*****************************************************************

	INPUT1='traindata.txt'
        INPUT2='responsetrain.txt'
        OPEN(unit=3,file=INPUT1,status='old')
        OPEN(unit=4,file=INPUT2,status='old')
        do im=1,nm
	 read(3,*)(prdata(im,it),it=1,nt)
        enddo


	 write(6,*)(prdata(1,it),it=1,nt)


	do id=1,nm
	 do ic=1,nc
	 y(id,ic)=0d0
	 enddo 
	enddo

	do id=1,43
	 y(id,1)=1d0
	enddo

	do id=44,nm
	 y(id,2)=1d0
	enddo

	do id=1,nm
	 write(6,*)(y(id,ic),ic=1,nc)
	enddo


	do id=1,nm
	 do it=1,nt
	  do iv=1,nv
	   wdt(id,it,iv)=0d0
	  enddo
	 enddo
	enddo

	do id=1,nm
	 do it=1,nt
	  if(prdata(id,it) .eq. 'a') then
	   wdt(id,it,1)=1d0
	  endif
	  if(prdata(id,it) .eq. 'g') then
	   wdt(id,it,2)=1d0
	  endif
	  if(prdata(id,it) .eq. 'c') then
	   wdt(id,it,3)=1d0
	  endif
	  if(prdata(id,it) .eq. 't') then
	   wdt(id,it,4)=1d0
	  endif
	 enddo
	enddo

	do it=1,nt
	 write(6,*)(wdt(1,it,iv),iv=1,nv)
	enddo



c	go to 9999

*****************************************************************
c set values for gamma(nm,nt,nj) 
*****************************************************************

      do id=1,nm
       do it=1,nt
	   sgamma=0d0
        do j=1,nj
         gamma(id,it,j)=1d0
	    sgamma=sgamma+gamma(id,it,j) 
        enddo
	   do j=1,nj
         gamma(id,it,j)=gamma(id,it,j)/sgamma
	   enddo
       enddo
      enddo
     
*****************************************************************
c set values for sigma2hat-inverse gamma prior
*****************************************************************

	 call DRNGAM(1,1d0,res1)
         s2h=1d0/res1(1)
c	write(6,*)'s2h',s2h


*****************************************************************
c set values for dhat-normal prior
*****************************************************************

	do id=1,nm
	 do it=1,nt
	  do j=1,nj
	   call DRNNOR(1,tempn(1))
	   dhat(id,it,j)=tempn(1)
	  enddo
	 enddo
	enddo

*****************************************************************
*             Initialize model parameters                       *
*****************************************************************


*****************************************************************
c Initial Sigma : sigma(nj,nj)
*****************************************************************

	sigma(1,1)=0.3714d0
	sigma(2,2)=0.1305d0
	sigma(3,3)=0.3809d0
	sigma(1,2)=0.12491d0
	sigma(2,1)=sigma(1,2)
	sigma(1,3)=0.1229d0
	sigma(3,1)=sigma(1,3)
	sigma(2,3)=0.1327d0
	sigma(3,2)=sigma(2,3)

	write(6,*)'Initial Sigma'
	do i=1,nj
	 write(6,*)(sigma(i,j),j=1,nj)
	enddo

*****************************************************************
c Initial Phi : phi(nj,nj)
*****************************************************************

	phi(1,1)=0.1250d0
	phi(1,2)=0.1207d0
	phi(1,3)=0.1278d0
	phi(2,1)=0.1308d0
	phi(2,2)=0.1317d0
	phi(2,3)=0.1312d0
	phi(3,1)=0.1200d0
	phi(3,2)=0.1197d0
	phi(3,3)=0.1281d0

	write(6,*)'Initial Phi'
	do i=1,nj
	 write(6,*)(phi(i,j),j=1,nj)
	enddo



*****************************************************************
c Initial matrix of word probabilities beta: beta(nj,nv)
*****************************************************************
	shp=1d0
	do j=1,nj
  	 call DRNGAM(nv,shp,tmpgam)
	 sum=0d0
	 do iv=1,nv
	  sum=sum+tmpgam(iv)
	 enddo
	 do iv=1,nv
	  beta(j,iv)=tmpgam(iv)/sum
	 enddo
	enddo

	write(6,*)'Initial beta'
	do i=1,nj
	 write(6,*)(beta(i,iv),iv=1,nv)
	enddo		
		
c check whether each row sum in beta is 1
	do j=1,nj
	 sum=0d0
	 do iv=1,nv
          sum=sum+beta(j,iv)
	 enddo
c	 write(6,*)j,sum
	enddo

*****************************************************************
c Initial eta: eta(nc,nj),each element of eta(CxJ) is from N(0,1) 
*****************************************************************


	eta(1,1)=0.5d0
	eta(1,2)=1.5d0
	eta(1,3)=1.5d0
	eta(2,1)=1.5d0
	eta(2,2)=1.5d0
	eta(2,3)=0.5d0


c	do ic=1,nc
c	 do j=1,nj
c	  eta(ic,j)=DRNUNF()
c	 enddo
c	enddo

	write(6,*)'Initial eta'
	do ic=1,nc
	 write(6,*)(eta(ic,j),j=1,nj)
	enddo

c End of setting initial values

*****************************************************************
c Calling each function
*****************************************************************


*****************************************************************
*                         VEM                                   *
*****************************************************************
*  E Step: For each document, find optimizing values of the	*
*          variational parameters.    				*
*  M Step: Maximize the resulting lower bound on the LL wrt     *
*         the model parameters.					*
****************************************************************
c Start EM iteration

	nmntnj=nm*nt*nj
	nmnt=nm*nt
	njnj=nj*nj
	ncnj=nc*nj
	ntnj=nt*nj
	nk=nj+(nj*(nj-1)/2)

c arguments general 
	ftol=1d-3
	fvalue=0d0
	ibtype=0
	maxfcn=500

c arguments for dhat- unconstrained optimization
c For ntnj variables for each doc
	do it=1,nt
	 do j=1,nj
	  intnj=((it-1)*nj)+j
	  xsc11(intnj)=1d0
	  dlb(intnj)=-10d0
	  dub(intnj)=-10d0
	 enddo
	enddo


c arguments for s2h- subject to bounds
	s2hg=s2h
	s2hlb=0.001d0
	s2hub=10d0

c arguments for gamma- NONE, fixed point est

c arguments for phi- unconstrained
c NONE

c arguments for sigma-constrained
c positive definite condition defined within f
	call  smat2vec(sigma,nj,sigmagvec) 
	do i=1,nk
	  slb(i)=0.0001d0
	  sub(i)=10d0	  
	enddo


	ntem=5
	ns=10

c arguments for beta- NONE, fixed point est

c arguments for eta- unconstrained

c ELBO before iter

	call elboLL(phi,sigma,beta,eta,gamma, 
     & dhat,s2h,wdt,y,elbo)
	write(6,*)'elbo',elbo

      call TDATE(iday,imonth,iyear)
      call TIMDY(ihour,iminute,isec)
      time0=CPSEC()
      write(6,*)'Current day: ',iday,imonth,iyear
      write(6,*)'Begin time:',ihour,iminute,isec

*****************************************************************
       do iter=1,nem
c	go to 6666
*****************************************************************
c E-step
*****************************************************************
	call elboLL(phi,sigma,beta,eta,gamma, 
     & dhat,s2h,wdt,y,elbop)

	write(6,*)'iter number:',iter
	write(6,*)'before dhat'

c Update dhat

	do id=1,nm
	 do it=1,nt
	  do j=1,nj
	   dhatd(it,j)= dhat(id,it,j)
	   gammad(it,j)=gamma(id,it,j)
	  enddo
	  enddo


	  call updatedhatd(dhatd,ntem,ns,dhatdout,fdout)
	
	  do it=1,nt
	   do j=1,nj
	    dhatd(it,j)=dhatdout(it,j)
	   enddo
	  enddo

          do it=1,nt
	   do j=1,nj
	    dhat(id,it,j)=dhatd(it,j)
	   enddo
	  enddo
         enddo

	write(6,*)(dhat(1,1,j),j=1,nj)
	write(6,*)'after dhat'
      call TDATE(iday,imonth,iyear)
      call TIMDY(ihour,iminute,isec)
      time01=CPSEC()
      write(6,*)'Current day: ',iday,imonth,iyear
      write(6,*)'Begin time:',ihour,iminute,isec


*****************************************************************
c Update s2h
	write(6,*)'before s2h'
	s2hg=s2h
	call DBCPOL(fs2h,1,s2hg,ibtype,s2hlb,s2hub,ftol,
     &    maxfcn,s2h,fs2hval)
      write(6,*)'s2h update',s2h
	write(6,*)'after s2h'
      call TDATE(iday,imonth,iyear)
      call TIMDY(ihour,iminute,isec)
      time02=CPSEC()
      write(6,*)'Current day: ',iday,imonth,iyear
      write(6,*)'Begin time:',ihour,iminute,isec
***************************************************************** 
c Update gamma
	write(6,*)'before gamma'

	call fgamma(gammaup)
	write(6,*)'after gamma'
	do id=1,nm
	 do it=1,nt
	  do j=1,nj
	  gamma(id,it,j)=gammaup(id,it,j)
	  enddo
	 enddo
	enddo
	write(6,*)(gamma(1,1,j),j=1,nj)
	write(6,*)(gamma(2,1,j),j=1,nj)
	write(6,*)'after gamma'
      call TDATE(iday,imonth,iyear)
      call TIMDY(ihour,iminute,isec)
      time03=CPSEC()
      write(6,*)'Current day: ',iday,imonth,iyear
      write(6,*)'Begin time:',ihour,iminute,isec

*****************************************************************
c M-step
*****************************************************************

c Update phi
	 write(6,*)'before phi'

	call updatephi(phi,ntem,ns,phiout,fout)

	do i=1,nj
	 do j=1,nj
	  phi(i,j)= phiout(i,j)
	 enddo
	enddo

	 write(6,*)'phi'	
	 do i=1,nj
	  write(6,*)(phi(i,j),j=1,nj)
	 enddo

	 write(6,*)'after phi'
        call TDATE(iday,imonth,iyear)
        call TIMDY(ihour,iminute,isec)
        time04=CPSEC()
        write(6,*)'Current day: ',iday,imonth,iyear
        write(6,*)'Begin time:',ihour,iminute,isec

*****************************************************************
c 6666	continue
c Update sigma
	write(6,*)'before sigma'

	call updatesig(sigma,ntem,ns,sigout,fsout)
	do i=1,nj
	 do j=1,nj
	  sigma(i,j)= sigout(i,j)
	 enddo
	enddo

	 write(6,*)'sigma'	
	 do i=1,nj
	  write(6,*)(sigma(i,j),j=1,nj)
	 enddo

	write(6,*)'after sigma'

        call TDATE(iday,imonth,iyear)
        call TIMDY(ihour,iminute,isec)
        time05=CPSEC()
        write(6,*)'Current day: ',iday,imonth,iyear
        write(6,*)'Begin time:',ihour,iminute,isec
*****************************************************************
c Update beta
	write(6,*)'before beta'
	call fbeta(beta)
	write(6,*)'after beta'
	 write(6,*)'beta'
	 do j=1,nj
	  write(6,*)(beta(j,iv),iv=1,nv)
	 enddo
        call TDATE(iday,imonth,iyear)
        call TIMDY(ihour,iminute,isec)
        time06=CPSEC()
        write(6,*)'Current day: ',iday,imonth,iyear
        write(6,*)'Begin time:',ihour,iminute,isec
*****************************************************************
c Update eta
	write(6,*)'before eta'
c	call  mat2vec(eta,nc,nj,etagvec)	
c	call DUMCGF(feta,ncnj,etagvec, xsc3,1d-3,
c     &   0, 0.01d0,etavec,g3,fetaval)


	call updateta(eta,ntem,ns,etaout,feout)
	

	do ic=1,nc
	 do j=1,nj
	  eta(ic,j)= etaout(ic,j)
	 enddo
	enddo

	write(6,*)'after eta'
	 write(6,*)'eta'
	 do ic=1,nc
	  write(6,*)(eta(ic,j),j=1,nj)
	 enddo
        call TDATE(iday,imonth,iyear)
        call TIMDY(ihour,iminute,isec)
        time07=CPSEC()
        write(6,*)'Current day: ',iday,imonth,iyear
        write(6,*)'Begin time:',ihour,iminute,isec
*****************************************************************
c6666	continue
	call elboLL(phi,sigma,beta,eta,gamma, 
     & dhat,s2h,wdt,y,elbon)
	elbon=elbon-100d0
	write(6,*)'end of iter',iter
	write(6,*)'elbo',elbon
	elbovec(iter)= elbon
	if(dabs(elbon-elbop) .le. 0.1d0) then
	write(6,*)'Converged'

	endif
*****************************************************************		 	  
       enddo

	write(6,*)'final elbo',elbon
*****************************************************************

      call TDATE(iday,imonth,iyear)
      call TIMDY(ihour,iminute,isec)
      time1=CPSEC()
      write(6,*)'Current day: ',iday,imonth,iyear
      write(6,*)'End time:',ihour,iminute,isec



	do id=1,nm
	 do j=1,nj
	  gsum=0d0
	  do it=1,nt
           gsum=gsum+gamma(id,it,j)
          enddo
	  gbar(j,1)= gsum/dfloat(nt)
	 enddo
 	 call MATMUL(eta,gbar,matr,nc,nj,1)
	 write(6,*)'gbar',(gbar(j,1),j=1,nj)
	 do ic=1,nc
	  yp(id,ic)=0d0
	  matr1(ic)= matr(ic,1)
	 enddo
         write(6,*)(matr1(ic),ic=1,nc)
 	 call whichmax(matr1,nc,ind)
	 write(6,*) ind
	 yp(id,ind)=1d0
	enddo

	write(6,*)'y pred'
	do id=1,nm
	 write(6,*)(yp(id,ic), ic=1,nc)
	enddo



	acc1=0d0
	do id=1,nm
	  if(y(id,1).eq. yp(id,1)) then
	   acc1= acc1+1d0
	  endif
	enddo
	write(6,*)'train count', acc1
	write(6,*)'train accuracy', acc1/(dfloat(nm))

c      	write(6,*)'elbo values'
c	write(6,*)(elbovec(i),i=1,nem)




c 9999	continue

	STOP
	END





*****************************************************************
*****************************************************************
*c Subroutines  for the E step					*
*****************************************************************
*****************************************************************


*****************************************************************
        subroutine fdhatd(ntnj,dhatdvec,fdhr)
*****************************************************************

c This is the objective function to maximize. 
c We consider negative to minimize
c inputs:  nt*nj, dhatvec
c outputs: fdhr


      implicit real*8 (a-h,o-z)
      parameter (nj=3,nv=4,nm=86,nt=57,nc=2)
      integer ntnj,d
      real*8 dhatdvec(nt*nj),fdhr,dhatv(nt,nj)
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 gammad(nt,nj)
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fmd(nt,nj,1),fvd(nt,nj,nj)
      real*8 smd(nt,nj,1),svd(nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9


	 common/vpar/gamma,dhat,s2h
	 common/mpar/phi,sigma,beta,eta
	 common/dat/wdt,y
	 common/gd/gammad


*****************************************************************
c dhatdvec to dhatv
*****************************************************************

c	write(6,*)'dhatvec'
c	write(6,*)(dhatvec(i),i=1,nmntnj)

	call vec2mat(dhatdvec,nt,nj,dhatv)


	call FFBSd(sigma,dhatv,s2h,phi,fmd,fvd,smd,svd)

*****************************************************************


*****************************************************************	          	
c term1
*****************************************************************
	term1=0d0
	 do it=1,nt
	  st1=0d0
           do j=1,nj
	    t11=dexp(smd(it,j,1)+(svd(it,j,j)/2))
	    st1=st1+t11
	   enddo
	   do j=1,nj
	    t12=gammad(it,j)*(smd(it,j,1)-log(st1))
	    term1=term1+t12
	   enddo
	 enddo


*****************************************************************
c term 3 
***************************************************************** 
      term3=0d0

	do it=1,1
	 do i=1,nj
	  smmat(i,1)=smd(it,i,1)
	  smmat1(i,1)=0d0
	 enddo
	 call MATMUL(phi,smmat1,mt31,nj,nj,1)
	 do i=1,nj
	  mt32(i,1)=smmat(i,1)-mt31(i,1)
	 enddo
	 call MATTRAN(mt32,tmt32,nj,1)
	 call DLINDS(nj,sigma,nj,siginv,nj)
	 call MATMUL(tmt32,siginv,mt33,1,nj,nj) 
	 call MATMUL(mt33,mt32,mt34,1,nj,1)
	 term3=term3+ mt34(1,1)
	enddo
	
       do it=2,nt
	 do i=1,nj
	  smmat(i,1)=smd(it,i,1)
	  smmat1(i,1)=smd(it-1,i,1)
	 enddo
	 call MATMUL(phi,smmat1,mt31,nj,nj,1)
	 do i=1,nj
	  mt32(i,1)=smmat(i,1)-mt31(i,1)
	 enddo
	 call MATTRAN(mt32,tmt32,nj,1)
	 call DLINDS(nj,sigma,nj,siginv,nj)
	 call MATMUL(tmt32,siginv,mt33,1,nj,nj) 
	 call MATMUL(mt33,mt32,mt34,1,nj,1)
	 term3=term3+ mt34(1,1)
	enddo


     
      

*****************************************************************
c fdhr
*****************************************************************

       fdhr = -(term1-(0.5*term3))
	
       return
       end

*****************************************************************
c output  dhatd
***************************************************************** 


*****************************************************************
	subroutine fdhat(nmntnj,dhatvec,fdhr)
*****************************************************************

c This is the objective function to maximize. 
c We consider negative to minimize in dbcpol
c inputs:  nj, dhatvec
c outputs: fdhr


      implicit real*8 (a-h,o-z)
      parameter (nj=3,nv=4,nm=86,nt=57,nc=2)
      integer nmntnj
      real*8 dhatvec(nm*nt*nj),fdhr,dhatv(nm,nt,nj)
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)


	 common/vpar/gamma,dhat,s2h
	 common/mpar/phi,sigma,beta,eta
	 common/dat/wdt,y

*****************************************************************
c dhatvec to dhat
*****************************************************************

c	write(6,*)'dhatvec'
c	write(6,*)(dhatvec(i),i=1,nmntnj)

	call vec2arr(dhatvec,nm,nt,nj,dhatv)

c	do id=1,nm
c	 do it=1,nt
c	  write(6,*)(dhatv(id, it,j),j=1,nj)
c	enddo
c     enddo

*****************************************************************

	call FFBS(sigma,dhatv,s2h,phi,fm,fv,sm,sv)

*****************************************************************


*****************************************************************	          	
c term1
*****************************************************************
	call term1f(sm,sv,gamma,term1)

*****************************************************************
c term 3 
***************************************************************** 
	call term3f(sm,phi,sigma,term3) 
     
      

*****************************************************************
c fdhr
*****************************************************************

       fdhr = -(term1-(0.5*term3))
	
       return
       end

*****************************************************************
c output  dhat
***************************************************************** 




*****************************************************************
      subroutine fs2h(ns,s2hv,fs2)
*****************************************************************


c This is the objective function to maximize. 
c We consider negative to minimize in dbcpol
c inputs:  1, s2h
c outputs: fs2

      implicit real*8 (a-h,o-z)
      parameter (nj=3,nv=4,nm=86,nt=57,nc=2)
      integer ns
      real*8 s2hv,fs2
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj)
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21, phit(nj,nj)
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43

	 common/vpar/gamma,dhat,s2h
	 common/mpar/phi,sigma,beta,eta
	 common/dat/wdt,y

*****************************************************************

	call FFBS(sigma,dhat,s2hv,phi,fm,fv,sm,sv)

*****************************************************************

	  
*****************************************************************	          	
c term1
*****************************************************************
	call term1f(sm,sv,gamma,term1)

*****************************************************************	 
c term2
*****************************************************************
 	call term2f(sv,phi,sigma,term2)


*****************************************************************
c term 3 
***************************************************************** 
	call term3f(sm,phi,sigma,term3) 
     
      
*****************************************************************
c term 4
*****************************************************************

	call term4f(sv,term4) 

*****************************************************************
c fs2
*****************************************************************

       fs2= -(term1-(0.5*term2)-(0.5*term3)+(0.5*term4))
	
       return
       end

*****************************************************************
c output s2h
*****************************************************************


*****************************************************************
	subroutine fgamma(gamu)
*****************************************************************


      implicit real*8 (a-h,o-z)
      parameter (nj=3,nv=4,nm=86,nt=57,nc=2)
      real*8 gam1(nm,nt,nj), gamu(nm,nt,nj) 
      real*8 gamsum(nm,nt)
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)	
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 prod,sum1,t7,t8
      real*8 g31,g32,g3,g2,g1
      real*8 vec1(nj)

	 common/vpar/gamma,dhat,s2h
	 common/mpar/phi,sigma,beta,eta
	 common/dat/wdt,y


*****************************************************************

	call FFBS(sigma,dhat,s2h,phi,fm,fv,sm,sv)

*****************************************************************


c id loop
      do id=1,nm
c Define hTgamma dt
       g31=0d0
       do ic=1,nc
        prod=1d0
        do it=1,nt
         sum1=0d0
         do j=1,nj
          t7=dexp(eta(ic,j)/dfloat(nt))
          t8 = t7*gamma(id,it,j)
          sum1=sum1+t8
         enddo
         prod=prod*sum1
        enddo
        g31=g31+prod
       enddo
c	write(6,*)'g31',g31

c Define hj
      do j=1,nj
	do k=1,nj
	 vec1(k)=0d0
	enddo
       vec1(j)=1d0
       g32=0d0
       do ic=1,nc
        prod=1d0
        do it=1,nt
         sum1=0d0
         do i=1,nj
          t7=dexp(eta(ic,i)/dfloat(nt))
          t8 = t7*vec1(i)
          sum1=sum1+t8
         enddo
         prod=prod*sum1
        enddo
        g32=g32+prod
       enddo
c	write(6,*)'g32',g32
        g3=g32/g31

c	write(6,*)'g3',g3
        do ic=1,nc
         if(y(id,ic).eq. 1d0) then
          g2=eta(ic,j)/dfloat(nt)
c	  write(6,*)ic,g2
         endif
        enddo
c	write(6,*)'g2',g2
c it loop
        do it=1,nt
         g1=0d0
         do iv=1,nv
          g1= g1+wdt(id,it,iv)*beta(j,iv)
         enddo
	 gam1(id,it,j)=g1* dexp(sm(id,it,j,1)+g2-g3)
        enddo
       enddo
      enddo

c	write(6,*)(gam1(1,1,j),j=1,nj)
	do id=1,nm
	 do it=1,nt
	  gamsum(id,it)=0d0
	  do j=1,nj
           gamsum(id,it)=gamsum(id,it)+gam1(id,it,j)
	  enddo
	 enddo
	enddo

	do id=1,nm
	 do it=1,nt
	  do j=1,nj
              gamu(id,it,j)=gam1(id,it,j)/gamsum(id,it)
	  enddo
	 enddo
	enddo

      return
      end 
      

*****************************************************************
c output gamma
*****************************************************************



*****************************************************************
*****************************************************************
*c Subroutines  for the M step					*
*****************************************************************
*****************************************************************


*****************************************************************
	subroutine fbeta(betau)
c this is to update beta in the M-step
*****************************************************************


      implicit real*8 (a-h,o-z)
      parameter (nj=3,nv=4,nm=86,nt=57,nc=2)
      real*8 beta1(nj,nv),betau(nj,nv)
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)	
      real*8 sum1(nj)

	 common/vpar/gamma,dhat,s2h
	 common/mpar/phi,sigma,beta,eta
	 common/dat/wdt,y


        do j=1,nj
	 sum1(j)=0d0
	 do iv=1,nv
	 beta1(j,iv)=0d0
	  do id=1,nm
	   do it=1,nt
	    beta1(j,iv)=beta1(j,iv)+gamma(id,it,j)*wdt(id,it,iv)
	   enddo
	  enddo
	  sum1(j)=sum1(j)+beta1(j,iv)
	 enddo
	enddo	

        do j=1,nj
	 do iv=1,nv
	  betau(j,iv)=beta1(j,iv)/sum1(j)
	 enddo
	enddo

      return
      end 

*********************************************************************   
c output beta
********************************************************************* 



*****************************************************************
	subroutine fphi(njnj,phivec,fp)
c this is to update phi in the M-step
*****************************************************************

c This is the objective function to maximize. 
c We consider negative to minimize in dbcpol
c inputs:  njnj, phivec
c outputs: fp
      
      implicit real*8 (a-h,o-z)
      parameter (nj=3,nv=4,nm=86,nt=57,nc=2)
      integer njnj
      real*8 phivec(nj*nj),fp,phiv(nj,nj)
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 term1,t1,st1,t11,t12, phitv(nj,nj)
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 fpold,ev1,fp1
      real*8 phi1(nj,nj),tsum
      complex*16 eigenp(nj),sum

	 common/vpar/gamma,dhat,s2h
	 common/mpar/phi,sigma,beta,eta
	 common/dat/wdt,y
 	 common/fphi1/fp1


*****************************************************************
c phivec to phi
*****************************************************************

c	write(6,*)'phivec'
c	write(6,*)(phivec(i),i=1,njnj)

	call vec2mat(phivec,nj,nj,phiv)

c	do i=1,nj
c	 write(6,*)(phiv(i,j),j=1,nj)
c	enddo



*****************************************************************

	call FFBS(sigma,dhat,s2h,phiv,fm,fv,sm,sv)

*****************************************************************
	  
*****************************************************************	          	
c term1
*****************************************************************
	call term1f(sm,sv,gamma,term1)


*****************************************************************	 
c term2
*****************************************************************
 	call term2f(sv,phiv,sigma,term2)

*****************************************************************
c term 3 
***************************************************************** 
	call term3f(sm,phiv,sigma,term3) 
     
      
*****************************************************************
c term 4
*****************************************************************

	call term4f(sv,term4) 
	tsum=0d0
	do i=1,nj
	 do j=1,nj
	 tsum=tsum+dabs(phiv(i,j))
	 enddo
	enddo
*****************************************************************
c fp
*****************************************************************

       fp= -(term1-(0.5*term2)-(0.5*term3)+(0.5*term4))
c       fp=fp+(dfloat(nm)*tsum)
       return
       end


*********************************************************************   
c output phi
*********************************************************************  




*****************************************************************
	subroutine feta(ncnj,etavec,fe)
c this is to update eta in the M-step
*****************************************************************

c This is the objective function to maximize.
c We consider negative to minimize in dbcpol
c inputs:  ncnj,etavec
c outputs: fe
      
      	
      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      integer ncnj
      real*8 etavec(nc*nj),fe,etav(nc,nj)
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2, tsum

	 common/vpar/gamma,dhat,s2h
	 common/mpar/phi,sigma,beta,eta
	 common/dat/wdt,y




*****************************************************************
c etavec to eta
*****************************************************************

c	write(6,*)'etavec'
c	write(6,*)(etavec(i),i=1,ncnj)

        call vec2mat(etavec,nc,nj,etav)

c	do i=1,nc
c	 write(6,*)(etav(i,j),j=1,nj)
c	enddo


*****************************************************************
c term6
*****************************************************************

	call term6f(gamma,etav,y,term6)

*****************************************************************	
c term7
*****************************************************************
	call term7f(etav,gamma,term7) 

c	penalty term
	tsum=0d0
	do ic=1,nc
	 do j=1,nj
	 tsum=tsum+(etav(ic,j)*etav(ic,j))
	 enddo
	enddo


*****************************************************************
c fe
*****************************************************************

      fe= -(term6-term7-tsum)
	
      return
      end



*****************************************************************
c output eta
*****************************************************************


***************************************************************** 
	subroutine fsig(nk,sigvec,fs)
c this is to update sigma in the M-step
*****************************************************************

c This is the objective function to maximize. 
c We consider negative to minimize in dbcpol
c inputs:  nk, sigvec
c outputs: fs

	
      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      integer nk
      real*8 sigvec(nj+(nj*(nj-1)/2)),fs,sigmav(nj,nj)
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj),fs1,fs11
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 fsold,eigens(nj),evmin, tsum

	 common/vpar/gamma,dhat,s2h
	 common/mpar/phi,sigma,beta,eta
	 common/dat/wdt,y
	 common/fsig1/fs11

*****************************************************************
c sigvec to sigma
*****************************************************************

c	write(6,*)'sigvec'
c	write(6,*)(sigvec(i),i=1,njnj)

	call svec2mat(sigvec,nj,sigmav)

c	do i=1,nj
c	 write(6,*)(sigmav(i,j),j=1,nj)
c	enddo

c Check if sigma is pd--- all eigen values >0

	call DEVLSF(nj,sigmav,nj,eigens)
	evmin = eigens(1)
        do j= 2,nj
	 evmin= dmin1(evmin,eigens(j))	
	enddo
	
	if(evmin.le.0d0) then
	 fs=fs11+1d0
	 go to 999	
	endif
	

c	if(not pd) set fs=fold
c	return to end
*****************************************************************

	call FFBS(sigmav,dhat,s2h,phi,fm,fv,sm,sv)

*****************************************************************


*****************************************************************	          	
c term1
*****************************************************************
	call term1f(sm,sv,gamma,term1)

*****************************************************************	 
c term2
*****************************************************************
 	call term2f(sv,phi,sigmav,term2)

*****************************************************************
c term 3 
***************************************************************** 
	call term3f(sm,phi,sigmav,term3) 
     
*****************************************************************
c term 4
*****************************************************************

	call term4f(sv,term4) 

*****************************************************************
c term 5
*****************************************************************
 	call term5f(sigmav,term5) 
	
	tsum=0d0
	do i=1,nj
	 do j=1,nj
	 tsum=tsum+dabs(sigmav(i,j))
	 enddo
	enddo

*****************************************************************
c fs
*****************************************************************

       fs1= -(term1-(0.5d0*term2)-(0.5d0*term3))
	  fs= fs1-((0.5d0*term4)-(0.5d0*nm*nt*term5))
c	fs11=fs+(dfloat(nt)*tsum)
	fs11=fs
	 
999    continue	
	fs11=fs
       return
       end

*****************************************************************
c output sigma
*****************************************************************

*****************************************************************
*****************************************************************
*c Subroutines to calculate ELBO on the L			*
      subroutine elboLL(phi,sigma,beta,eta,gamma, 
     & dhat,s2h,wdt,y,elbo)
*****************************************************************
*****************************************************************

c Input  : phi,sigma,beta,eta,gamma,dhat,s2hwdt,y
c Output : ELBO
	
      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9
      real*8 fs1,fs2,elbo

	
*****************************************************************

	call FFBS(sigma,dhat,s2h,phi,fm,fv,sm,sv)

*****************************************************************


*****************************************************************	          	
c term1
*****************************************************************
	call term1f(sm,sv,gamma,term1)

*****************************************************************	 
c term2
*****************************************************************
 	call term2f(sv,phi,sigma,term2)

*****************************************************************
c term 3 
***************************************************************** 
	call term3f(sm,phi,sigma,term3) 
     
*****************************************************************
c term 4
*****************************************************************

	call term4f(sv,term4) 

*****************************************************************
c term 5
*****************************************************************
 	call term5f(sigma,term5) 
	
*****************************************************************
c term 6
*****************************************************************
	call term6f(gamma,eta,y,term6)

*****************************************************************	
c term7
*****************************************************************
	call term7f(eta,gamma,term7) 

*****************************************************************	
c term8
*****************************************************************

	call term8f(gamma,term8) 


*****************************************************************	
c term9
*****************************************************************
	call term9f(gamma,wdt,beta,term9) 




*****************************************************************
c ELBO
*****************************************************************

       fs1= -(term1-(0.5d0*term2)-(0.5d0*term3))
       fs2= fs1-((0.5d0*term4)-(0.5d0*nm*nt*term5))
       elbo= fs2-(term6-term7-term8+term9)
	 
	 
	
       return
       end

*****************************************************************
c output elbo
*****************************************************************




*****************************************************************   
	subroutine mat2vec(mat,m,n,mvec)
*****************************************************************
c inputs:  mat,m,n
c outputs: mvec

      implicit real*8 (a-h,o-z)
      integer m,n
      real*8 mat(m,n),mvec(m*n)
	
      do im=1,m
       do in=1,n
        imn=((im-1)*n)+in
        mvec(imn)= mat(im,in)
       enddo
      enddo

      return
      end 



***************************************************************** 
	subroutine vec2mat(vec,m,n,vmat)
*****************************************************************
c inputs:  vec,m,n
c outputs: vmat

      implicit real*8 (a-h,o-z)
      integer m,n
      real*8 vmat(m,n),vec(m*n)
	
      do im=1,m
       do in=1,n
        imn=((im-1)*n)+in
        vmat(im,in)=vec(imn)
       enddo
      enddo

      return
      end 

*****************************************************************   
	subroutine smat2vec(mat,m,mvec)
*****************************************************************
c inputs:  mat,m
c outputs: mvec

      implicit real*8 (a-h,o-z)
      integer m
      real*8 mat(m,m),mvec(m+(m*(m-1)/2))


      ik=1	
      do im=1,m
       do j=im,m
        mvec(ik)= mat(im,j)
	ik=ik+1
       enddo
      enddo

      return
      end 


*****************************************************************   
	subroutine svec2mat(mvec,m,mat)
*****************************************************************
c inputs:  mvec,m
c outputs: mat

      implicit real*8 (a-h,o-z)
      integer m
      real*8 mat(m,m),mvec(m+(m*(m-1)/2))


      ik=1	
      do im=1,m
       do j=im,m
        mat(im,j)=mvec(ik)
        mat(j,im)=mvec(ik)
	ik=ik+1
       enddo
      enddo

      return
      end 

*****************************************************************  
	subroutine MAKEI(n, matiden)
*****************************************************************
c inputs:  n
c outputs: matiden

      implicit real*8 (a-h,o-z)
      integer n
	real*8 matiden(n,n)
	
      do i=1,n
       do j=1,n
	   matiden(i,j)=0d0
	   matiden(i,i)=1d0
       enddo
      enddo

      return
      end 
 

***************************************************************** 
	subroutine MATMUL(a,b,c,na,nb,nc)
*****************************************************************
       implicit real*8 (a-h,o-z)
       real*8 a(na,nb),b(nb,nc),c(na,nc)
            do ia=1,na
             do ic=1,nc
               c(ia,ic)=0.0d0
                 do ib=1,nb
                   c(ia,ic)=c(ia,ic)+a(ia,ib)*b(ib,ic)
                 enddo
             enddo
            enddo
          return
          end


*****************************************************************
      subroutine MATTRAN(x,y,n,p)
*****************************************************************
c routine to take the transpose of an n*p matrix x 
c and return p*n matrix y
	implicit real*8 (a-h,o-z)
      integer n,p,i,j
      double precision x(n,p),y(p,n)
 
      do i=1,n
         do j=1,p
            y(j,i)=x(i,j)
         enddo
      enddo
      return
      end


*****************************************************************
      subroutine MATTRACE(a,na,tr)
*****************************************************************
c routine to obtain trace of a n*n matrix a 
c and return trace tr 
	implicit real*8 (a-h,o-z)
      integer n
      real*8 a(na,na), tr
 
	tr=0d0
      do i=1,na
         tr=tr+a(i,i)
      enddo

      return
      end

*****************************************************************  
	subroutine dotprod(a,b,n,result)
*****************************************************************
c inputs:  a,b,n
c outputs: result

	implicit real*8 (a-h,o-z)
	integer n
	real*8 a(n),b(n),result
	
	result=0d0
	do i=1,n
    	 result=result+a(i)*b(i)
	enddo

	return
	end 


*****************************************************************
	subroutine arr2vec(arr,m,n,k,vec)
*****************************************************************  

c inputs:  arr,m,n,k
c outputs: vec

      implicit real*8 (a-h,o-z)
      integer m,n,k
      real*8 arr(m,n,k), vec(m*n*k)
	
      do im=1,m
       do in=1,n
	do ik=1,k
         imnk=((im-1)*n*k)+((in-1)*k)+ik
         vec(imnk)= arr(im,in,ik)
	enddo
       enddo
      enddo

      return
      end 


*****************************************************************
	subroutine vec2arr(vec,m,n,k,arr)
*****************************************************************

c inputs:  vec,m,n,k
c outputs: arr

      implicit real*8 (a-h,o-z)
      integer m,n,k
      real*8 vec(m*n*k),arr(m,n,k)
	
      do im=1,m
       do in=1,n
	do ik=1,k
         imnk=((im-1)*n*k)+((in-1)*k)+ik
         arr(im,in,ik)=vec(imnk)
	enddo
       enddo
      enddo

      return
      end 

*****************************************************************
	subroutine FFBS(sigma,dhat,s2h,phi,fm,fv,sm,sv)
*****************************************************************
c Input  : sigma,dhat,s2h,phi
c Output : fm,fv,sm,sv


      implicit real*8 (a-h,o-z)
      parameter (nj=3,nv=4,nm=86,nt=57,nc=2)
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 sigma(nj,nj),gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 phi(nj,nj)
      real*8 phit(nj,nj),K(nj,nj),IK(nj,nj), dh(nj,1)
      real*8 m1(nj,nj),m2(nj,nj),m3(nj,nj),m4(nj,nj),m5(nj,nj)
      real*8 m6(nj,nj),m7(nj,1),m8(nj,1),m9(nj,nj)
      real*8 m10(nj,nj),m11(nj,nj),m12(nj,nj)
      real*8 m13(nj,nj),m14(nj,nj),m15(nj,nj)
      real*8 m16(nj,nj),m17(nj,nj),m18(nj,1)
      real*8 m19(nj,1),m20(nj,nj),m21(nj,nj)
      real*8 J1(nj,nj), tJ(nj,nj),m22(nj,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 fm1(nj,1),fm2(nj,1)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 m0(nj),V0(nj,nj),mp(nj,1),Vp(nj,nj)
      real*8 mn(nj,1),Vn(nj,nj),V1(nj,nj)


*****************************************************************
c initialize m0 and V0
*****************************************************************
      do i=1,nj
	m0(i)=0d0
        do j=1,nj
	 V0(i,j)=0d0
	enddo
       V0(i,i)=0.01d0
      enddo

c	write(6,*)'V0'
c            do j=1,nj
c             write(6,*)(V0(i,j),I=1,nj)
c	    enddo


*****************************************************************
c find filter mean fm(nm,nt,nj,1) and filtervcov fv(nm,nt,nj,nj)
*****************************************************************


      do id=1,nm
       do i=1,nj
	   mp(i,1)=m0(i)
	   do j=1,nj
	    Vp(i,j)=V0(i,j)
	   enddo
	  enddo

       do it=1,nt

c         write(6,*)'mp'
c          write(6,*)(mp(j,1),j=1,nj)

	   do j=1,nj
	    dh(j,1)=dhat(id,it,j)
	   enddo
	   call MATMUL(phi,Vp,m1,nj,nj,nj)
	   call MATTRAN(phi,phit,nj,nj)
	   call MATMUL(m1,phit,m2,nj,nj,nj)
	 
c	write(6,*)'m2 matrix'
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)m2(i,j)
c	     enddo
c	   enddo

c	write(6,*)'sigma matrix'
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)sigma(i,j)
c	     enddo
c	   enddo


        
	   do i=1,nj
	    do j=1,nj
	     m3(i,j)= m2(i,j)+sigma(i,j)
	     m4(i,j)= m2(i,j)+sigma(i,j)
             m4(i,i)= m2(i,i)+sigma(i,i)+s2h
	    enddo
	   enddo

c	write(6,*)'m3'
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)m3(i,j)
c	     enddo
c	   enddo

        
	 
c	write(6,*)'m4 matrix before dlinds call'
c            do j=1,nj
c             write(6,*)(m4(i,j),i=1,nj)
c	    enddo

c	   call DLINRG(nj,m4,nj,m5,nj)
	   call DLINDS(nj,m4,nj,m5,nj)


c	write(6,*)'Filtering iter id,it',id,it

c	write(6,*)'m5 matrix after dlinds call'
c	     do j=1,nj
c             write(6,*)(m5(i,j),i=1,nj)
c	    enddo
   
           Call MATMUL(m3,m5,K,nj,nj,nj)
	   do i=1,nj
	    do j=1,nj
	     IK(i,j)= -K(i,j)
	     IK(i,i)= 1d0-K(i,i)
	    enddo
	   enddo
	   call MATMUL(IK,phi,m6,nj,nj,nj)
	   call MATMUL(m6,mp,m7,nj,nj,1)
	   call MATMUL(K,dh,m8,nj,nj,1)
	   do j=1,nj
	    fm(id,it,j,1)=m7(j,1)+m8(j,1)
          enddo
	  call MATMUL(IK,m3,m9,nj,nj,nj)
          do i=1,nj
	    do j=1,nj
	     fv(id,it,i,j)=m9(i,j)
	    enddo
	   enddo
c	   write(6,*)'fm'
c           write(6,*)(fm(id,it,j,1),j=1,nj)
       
	   do i=1,nj
	    mp(i,1)=fm(id,it,i,1)
	    do j=1,nj
	     Vp(i,j)= fv(id,it,i,j)
            enddo
	   enddo
	  enddo
      enddo

*****************************************************************	   
c find smooth mean sm(nm,nt,nj,1) and smooth vcov sv(nm,nt,nj,nj)
*****************************************************************

      do id=1,nm
       do i=1,nj
	 sm(id,nt,i,1)=fm(id,nt,i,1)
	 mn(i,1)=fm(id,nt,i,1)
	 do j=1,nj
	  sv(id,nt,i,j)=fv(id,nt,i,j)
	  Vn(i,j)=fv(id,nt,i,j)
	 enddo
	enddo
c	write(6,*)'Vn'
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)Vn(i,j)
c	     enddo
c	   enddo

	do it=nt-1,1,-1
	  do i=1,nj
	   do j=1,nj
	    V1(i,j)=fv(id,it,i,j)
	   enddo
	  enddo
	 call MATMUL(V1,phit,m10,nj,nj,nj)
	 call MATMUL(phi,V1,m11,nj,nj,nj)
	 call MATMUL(m11,phit,m12,nj,nj,nj)

c	write(6,*)'m12 matrix '
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)m12(i,j)
c	     enddo
c	   enddo

	 do i=1,nj
	  do j=1,nj
	   m13(i,j)=m12(i,j)+sigma(i,j)
	  enddo
	 enddo

c	write(6,*)'m13 matrix before dlinds call'
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)m13(i,j)
c	     enddo
c	   enddo

	 call DLINDS(nj,m13,nj,m14,nj)

c	write(6,*)'Smoothing Iter id,it',id,it
c	write(6,*)'m14 matrix after dlinds call'
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)m14(i,j)
c	     enddo
c	   enddo


	 call MATMUL(m10,m14,m15,nj,nj,nj)
	 do i=1,nj
	  do j=1,nj
	   J1(i,j)=m15(i,j)
	  enddo
	 enddo
	 call MATMUL(J1,phi,m16,nj,nj,nj)     
	 do i=1,nj
	  do j=1,nj
	   m17(i,j)= -m16(i,j)
	   m17(i,i)= 1d0-m16(i,i)
	  enddo
	 enddo
	 do i=1,nj
	  fm1(i,1)= fm(id,it,i,1)
	  fm2(i,1)= fm(id,it+1,i,1)
	 enddo
	 call MATMUL(m17,fm1,m18,nj,nj,1)  
	 call MATMUL(J1,fm2,m19,nj,nj,1) 
	 do i=1,nj
	  sm(id,it,i,1)=m18(i,1)+m19(i,1)
         enddo
	 do i=1,nj
	  do j=1,nj
	   m20(i,j)= Vn(i,j)-m12(i,j)-sigma(i,j)
	  enddo
	 enddo
	 call MATMUL(J1,m20,m21,nj,nj,nj)
	 call MATTRAN(J1,tJ,nj,nj)
	 call MATMUL(m21,tJ,m22,nj,nj,nj)
	 do i=1,nj
	  do j=1,nj
	  sv(id,it,i,j)=V1(i,j)+ m22(i,j)
	  enddo
	 enddo
	 do i=1,nj
	  mn(i,1)= sm(id,it,i,1)
	  do j=1,nj
	   Vn(i,j)=sv(id,it,i,j)
	  enddo
	 enddo   
        enddo
      enddo



      return
      end 
*****************************************************************


*****************************************************************	          	
	subroutine term1f(sm,sv,gamma,term1)
*****************************************************************

      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9

	term1=0d0
	do id=1,nm
	 do it=1,nt
	  st1=0d0
           do j=1,nj
	    t11=dexp(sm(id,it,j,1)+(sv(id,it,j,j)/2))
	    st1=st1+t11
	   enddo
	   do j=1,nj
	    t12=gamma(id,it,j)*(sm(id,it,j,1)-log(st1))
	    term1=term1+t12
	   enddo
	 enddo
	enddo
c	write(6,*)'term1',term1

      return
      end 

*****************************************************************	 
	subroutine term2f(sv,phi,sigma,term2)
*****************************************************************

      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9

      term2=0d0
      do id=1,nm
	  do it=1,nt
	   do i=1,nj
	    do j=1,nj
             svmat(i,j)=sv(id,it,i,j)
	    enddo
	   enddo
	   call MATMUL(phi,svmat,mt21,nj,nj,nj)
	   call MATTRAN(phi,phit,nj,nj)
	   call MATMUL(mt21,phit,mt22,nj,nj,nj)
	   do i=1,nj
	    do j=1,nj
             mt23(i,j)=sv(id,it,i,j)+mt22(i,j)
	    enddo
	   enddo
	   call DLINDS(nj,sigma,nj,siginv,nj)
	   call MATMUL(siginv,mt23,mt24,nj,nj,nj)
	   call MATTRACE(mt24,nj,t21)
	   term2=term2+t21
	  enddo
        enddo

      return
      end 



*****************************************************************
	subroutine term3f(sm,phi,sigma,term3) 
*****************************************************************

      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9


      term3=0d0
      do id=1,nm
       do it=1,1
	do i=1,nj
	 smmat(i,1)=sm(id,it,i,1)
	 smmat1(i,1)=0d0
	enddo
	call MATMUL(phi,smmat1,mt31,nj,nj,1)
	do i=1,nj
	 mt32(i,1)=smmat(i,1)-mt31(i,1)
	enddo
	call MATTRAN(mt32,tmt32,nj,1)
	call DLINDS(nj,sigma,nj,siginv,nj)
	call MATMUL(tmt32,siginv,mt33,1,nj,nj) 
	call MATMUL(mt33,mt32,mt34,1,nj,1)
	term3=term3+ mt34(1,1)
       enddo
       do it=2,nt
	do i=1,nj
	 smmat(i,1)=sm(id,it,i,1)
	 smmat1(i,1)=sm(id,it-1,i,1)
	enddo
	call MATMUL(phi,smmat1,mt31,nj,nj,1)
	do i=1,nj
	 mt32(i,1)=smmat(i,1)-mt31(i,1)
	enddo
	call MATTRAN(mt32,tmt32,nj,1)
	call DLINDS(nj,sigma,nj,siginv,nj)
	call MATMUL(tmt32,siginv,mt33,1,nj,nj) 
	call MATMUL(mt33,mt32,mt34,1,nj,1)
	term3=term3+ mt34(1,1)
       enddo
      enddo
      
c      write(6,*)'term3',term3
      return
      end 



*****************************************************************
 	subroutine term4f(sv,term4) 
*****************************************************************

      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9

      term4=0d0
      do id=1,nm
       do it=1,nt
        do i=1,nj
         do j=1,nj
          mat41(i,j)=sv(id,it,i,j)
         enddo
        enddo

c	write(6,*)'id,it',id,it
c	write(6,*)'mat41'
c	   do i=1,nj
c	     do j=1,nj
c             write(6,*) mat41(i,j)
c	     enddo
c	   enddo

        call DLFTDS(nj,mat41,nj,mat42,nj)
        call DLFDDS(nj,mat42,nj,t41,t42)
        t43=t41*(10**t42)
        term4=term4+dlog(t43)
       enddo
       enddo
c      write(6,*)'term4',term4
      return
      end 

*****************************************************************
 	subroutine term5f(sigma,term5) 
*****************************************************************

      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9

      call DLFTDS(nj,sigma,nj,mat51,nj)
      call DLFDDS(nj,mat51,nj,t51,t52)
	t53=t51*(10**t52)
       term5=dlog(t53)
c      write(6,*)'term5',term5

      return
      end 


*****************************************************************
	subroutine term6f(gamma,eta,y,term6) 
*****************************************************************
      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9


c find gamma mean

    	 do i=1,nm
       	  do j=1,nj
           gammasum=0d0
            do it=1,nt
             gammasum=gammasum+gamma(i,it,j)
            enddo
            gammamean(i,j)=gammasum/dfloat(nt)
         enddo
        enddo


      term6=0d0
       do id=1,nm
        do j=1,nj
         gamma1(j)=gammamean(id,j)
         do ic=1,nc
          if(y(id,ic).eq.1d0) then
           eta1(j)=eta(ic,j)
          endif
         enddo
        enddo
        call dotprod(eta1,gamma1,nj,result12)
        term6=term6+result12
       enddo
c	write(6,*)'term6',term6

      return
      end 



*****************************************************************	
	subroutine term7f(eta,gamma,term7) 
*****************************************************************

      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9

      term7=0d0
      do id=1,nm
       result2=0d0
       do ic=1,nc
        prod=1d0
        do it=1,nt
         sum1=0d0
         do j=1,nj
         t7=dexp(eta(ic,j)/dfloat(nt))
          t8 = t7*gamma(id,it,j)
	  sum1=sum1+t8
         enddo
         prod=prod*sum1
        enddo
        result2=result2+prod
       enddo
       term7=term7+dlog(result2)
      enddo
c	write(6,*)'term7',term7

      return
      end 


*****************************************************************	
	subroutine term8f(gamma,term8) 
*****************************************************************

      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9


      term8=0d0
      do id=1,nm
       do it=1,nt
	do j=1,nj
         term8=term8+(gamma(id,it,j)*dlog(gamma(id,it,j)))
        enddo
       enddo
      enddo

c      write(6,*)'term8',term8


      return
      end 



*****************************************************************	
	subroutine term9f(gamma,wdt,beta,term9) 
*****************************************************************

      implicit real*8 (a-h,o-z)
      parameter(nj=3,nv=4,nm=86,nt=57,nc=2)	
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 gamma(nm,nt,nj),dhat(nm,nt,nj),s2h
      real*8 sigma(nj,nj),phi(nj,nj),beta(nj,nv),eta(nc,nj)
      real*8 fm(nm,nt,nj,1),fv(nm,nt,nj,nj)
      real*8 sm(nm,nt,nj,1),sv(nm,nt,nj,nj)
      real*8 phit(nj,nj)
      real*8 term1,t1,st1,t11,t12
      real*8 term2,svmat(nj,nj),mt21(nj,nj)
      real*8 mt22(nj,nj),mt23(nj,nj),siginv(nj,nj)
      real*8 mt24(nj,nj),t21
      real*8 term3,smmat(nj,1),smmat1(nj,1),mt31(nj,1)
      real*8 mt32(nj,1),tmt32(1,nj),mt33(1,nj),mt34(1,1)
      real*8 term4, mat41(nj,nj), mat42(nj,nj),t41,t42,t43
      real*8 term5,mat51(nj,nj),mat52(nj,nj),t51,t52,t53
      real*8 eta1(nj),gamma1(nj)
      real*8 gammasum,gammamean(nm,nj)
      real*8 term6,term7
      real*8 prod,sum1,t7,t8,result2
      real*8 term8,term9


      term9=0d0
      do id=1,nm
       do it=1,nt
	do j=1,nj
         do iv=1,nv
	  term9=term9+gamma(id,it,j)*wdt(id,it,iv)*dlog(beta(j,iv))
	 enddo
        enddo
       enddo
      enddo

c	write(6,*)'term9',term9


      return
      end 


*****************************************************************
      subroutine FFBSd(sigma,dhatd,s2h,phi,fmd,fvd,smd,svd)
*****************************************************************
c Input  : sigma,dhatd,s2h,phi
c Output : fmd,fvd,smd,svd


      implicit real*8 (a-h,o-z)
      parameter (nj=3,nv=4,nm=86,nt=57,nc=2)
      integer d
      real*8 y(nm,nc),wdt(nm,nt,nv)			
      real*8 sigma(nj,nj),gamma(nm,nt,nj),dhatd(nt,nj),s2h
      real*8 phi(nj,nj)
      real*8 phit(nj,nj),K(nj,nj),IK(nj,nj), dh(nj,1)
      real*8 m1(nj,nj),m2(nj,nj),m3(nj,nj),m4(nj,nj),m5(nj,nj)
      real*8 m6(nj,nj),m7(nj,1),m8(nj,1),m9(nj,nj)
      real*8 m10(nj,nj),m11(nj,nj),m12(nj,nj)
      real*8 m13(nj,nj),m14(nj,nj),m15(nj,nj)
      real*8 m16(nj,nj),m17(nj,nj),m18(nj,1)
      real*8 m19(nj,1),m20(nj,nj),m21(nj,nj)
      real*8 J1(nj,nj),tJ(nj,nj),m22(nj,nj)
      real*8 fmd(nt,nj,1),fvd(nt,nj,nj)
      real*8 fm1(nj,1),fm2(nj,1)
      real*8 smd(nt,nj,1),svd(nt,nj,nj)
      real*8 m0(nj),V0(nj,nj),mp(nj,1),Vp(nj,nj)
      real*8 mn(nj,1),Vn(nj,nj),V1(nj,nj)


*****************************************************************
c initialize m0 and V0
*****************************************************************
      do i=1,nj
	m0(i)=0d0
        do j=1,nj
	 V0(i,j)=0d0
	enddo
       V0(i,i)=0.01d0
      enddo

c	write(6,*)'V0'
c            do j=1,nj
c             write(6,*)(V0(i,j),i=1,nj)
c	    enddo


*****************************************************************
c find filter mean fmd(nt,nj,1) and filtervcov fvd(nt,nj,nj)

       do i=1,nj
	   mp(i,1)=m0(i)
	   do j=1,nj
	    Vp(i,j)=V0(i,j)
	   enddo
	  enddo

       do it=1,nt

c         write(6,*)'mp'
c          write(6,*)(mp(j,1),j=1,nj)

	   do j=1,nj
	    dh(j,1)=dhatd(it,j)
	   enddo
	   call MATMUL(phi,Vp,m1,nj,nj,nj)
	   call MATTRAN(phi,phit,nj,nj)
	   call MATMUL(m1,phit,m2,nj,nj,nj)

	   do i=1,nj
	    do j=1,nj
	     m3(i,j)= m2(i,j)+sigma(i,j)
	     m4(i,j)= m2(i,j)+sigma(i,j)
             m4(i,i)= m2(i,i)+sigma(i,i)+s2h
	    enddo
	   enddo

c	write(6,*)'m4 matrix before dlinds call'
c            do j=1,nj
c             write(6,*)(m4(i,j),i=1,nj)
c	    enddo

c	   call DLINRG(nj,m4,nj,m5,nj)
	   call DLINDS(nj,m4,nj,m5,nj)




c	write(6,*)'m5 matrix after dlinds call'
c	     do j=1,nj
c             write(6,*)(m5(i,j),i=1,nj)
c	    enddo
   
           Call MATMUL(m3,m5,K,nj,nj,nj)
	   do i=1,nj
	    do j=1,nj
	     IK(i,j)= -K(i,j)
	     IK(i,i)= 1d0-K(i,i)
	    enddo
	   enddo
	   call MATMUL(IK,phi,m6,nj,nj,nj)
	   call MATMUL(m6,mp,m7,nj,nj,1)
	   call MATMUL(K,dh,m8,nj,nj,1)
	   do j=1,nj
	    fmd(it,j,1)=m7(j,1)+m8(j,1)
          enddo
	  call MATMUL(IK,m3,m9,nj,nj,nj)
          do i=1,nj
	    do j=1,nj
	     fvd(it,i,j)=m9(i,j)
	    enddo
	   enddo
       
	   do i=1,nj
	    mp(i,1)=fmd(it,i,1)
	    do j=1,nj
	     Vp(i,j)= fvd(it,i,j)
            enddo
	   enddo
	  enddo

c	write(6,*)'Filter done'

*****************************************************************	   
c find smooth mean smd(nt,nj,1) and smooth vcov svd(nt,nj,nj)
*****************************************************************

       do i=1,nj
	 smd(nt,i,1)=fmd(nt,i,1)
	 mn(i,1)=fmd(nt,i,1)
	 do j=1,nj
	  svd(nt,i,j)=fvd(nt,i,j)
	  Vn(i,j)=fvd(nt,i,j)
	 enddo
	enddo
c	write(6,*)'Vn'
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)Vn(i,j)
c	     enddo
c	   enddo

	do it=nt-1,1,-1
	  do i=1,nj
	   do j=1,nj
	    V1(i,j)=fvd(it,i,j)
	   enddo
	  enddo
	 call MATMUL(V1,phit,m10,nj,nj,nj)
	 call MATMUL(phi,V1,m11,nj,nj,nj)
	 call MATMUL(m11,phit,m12,nj,nj,nj)

c	write(6,*)'m12 matrix '
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)m12(i,j)
c	     enddo
c	   enddo

	 do i=1,nj
	  do j=1,nj
	   m13(i,j)=m12(i,j)+sigma(i,j)
	  enddo
	 enddo

c	write(6,*)'m13 matrix before dlinds call'
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)m13(i,j)
c	     enddo
c	   enddo
c	write(6,*) it

	 call DLINDS(nj,m13,nj,m14,nj)

c	write(6,*)'m14 matrix after dlinds call'
c	   do i=1,nj
c	     do j=1,nj
c              write(6,*)m14(i,j)
c	     enddo
c	   enddo


	 call MATMUL(m10,m14,m15,nj,nj,nj)
	 do i=1,nj
	  do j=1,nj
	   J1(i,j)=m15(i,j)
	  enddo
	 enddo
	 call MATMUL(J1,phi,m16,nj,nj,nj)     
	 do i=1,nj
	  do j=1,nj
	   m17(i,j)= -m16(i,j)
	   m17(i,i)= 1d0-m16(i,i)
	  enddo
	 enddo
	 do i=1,nj
	  fm1(i,1)= fmd(it,i,1)
	  fm2(i,1)= fmd(it+1,i,1)
	 enddo
	 call MATMUL(m17,fm1,m18,nj,nj,1)  
	 call MATMUL(J1,fm2,m19,nj,nj,1) 
	 do i=1,nj
	  smd(it,i,1)=m18(i,1)+m19(i,1)
         enddo
	 do i=1,nj
	  do j=1,nj
	   m20(i,j)= Vn(i,j)-m12(i,j)-sigma(i,j)
	  enddo
	 enddo
	 call MATMUL(J1,m20,m21,nj,nj,nj)
	 call MATTRAN(J1,tJ,nj,nj)
	 call MATMUL(m21,tJ,m22,nj,nj,nj)
	 do i=1,nj
	  do j=1,nj
	  svd(it,i,j)=V1(i,j)+ m22(i,j)
	  enddo
	 enddo
	 do i=1,nj
	  mn(i,1)= smd(it,i,1)
	  do j=1,nj
	   Vn(i,j)=svd(it,i,j)
	  enddo
	 enddo   
        enddo



      return
      end 
*****************************************************************



*****************************************************************
      subroutine updatephi(phi,ntem,ns,phiout,fout)
*****************************************************************
	implicit real*8 (a-h,o-z)
	integer ntem,ns,ri,ci,neps,maxeval
        parameter (nj=3,nv=4,nm=86,nt=57,nc=2)
	real*8 tempr(1),dif1(1),dif2(1)
	real*8 phi(nj,nj),phi1(nj,nj),phiout(nj,nj)
	real*8 phivec(nj*nj),phi1vec(nj*nj)
	real*8 fp,fp1,fout,fold
	real*8 v(nj,nj)
	real*8 p1,p2,c1,eps,rt
	real*8 temp,count,ratio
	real*8 diff


	c1=1d0
	eps=0.001d0
	rt=0.25d0
	temp=1d0
	maxeval=3
	diff=0.001d0
	
	do i=1,nj
	 do j=1,nj
	  phiout(i,j)=phi(i,j)
	  v(i,j)= 0.001d0 
	 enddo
	enddo
	call  mat2vec(phi,nj,nj,phivec) 
	call  fphi(nj*nj,phivec,fp)
	fout=fp
	
	fold=fout
	iter=1
       do while(iter .le. maxeval)
	do in=1,ntem
	 count=0d0
	 do is=1,ns
	  call perturbp(phi,phi1,v)
	  call  mat2vec(phi1,nj,nj,phi1vec) 
	  call  fphi(nj*nj,phi1vec,fp1)
          if(fp1<fp) then
	   do i=1,nj
	    do j=1,nj
	     phi(i,j)=phi1(i,j)
	    enddo
	   enddo
	   dif1(1)= fp-fp1
	   fp=fp1
	   count=count+1d0
	  endif
          if(fp1<fout) then
	   do i=1,nj
	    do j=1,nj
	     phiout(i,j)=phi1(i,j)
	    enddo
	   enddo
	   dif2(1)= fout-fp1
	   fout=fp1
	  endif
          if(fp1.ge.fp) then
	   p1=DRNUNF()
           p2=dexp((fp1-fp)/temp)
	   if(p2 .gt. p1) then
	    do i=1,nj
	     do j=1,nj
	      phi(i,j)=phi1(i,j)
	     enddo
	    enddo
	    fp=fp1
	   endif
	  endif
	 enddo
	 ratio=count/dfloat(ns)
	 if(ratio .gt. 0.6d0) then
	  do i=1,nj
	   do j=1,nj
	    v(i,j)= v(i,j)*(1d0+c1*(ratio-0.6d0)/0.4d0)
	    if(v(i,j) .gt. diff) then
	     v(i,j)= diff
	    endif
	   enddo
	  enddo
	 else if(ratio .lt. 0.4d0) then
	  do i=1,nj
	   do j=1,nj
	    v(i,j)= v(i,j)/(1d0+c1*((0.4d0-ratio)/0.4d0))
	   enddo
	  enddo
	 endif
        enddo
	iter=iter+1
	 if(dif1(1)<eps .and. dif2(1)<eps) then
	  go to 4444
	 else
	  temp=rt*temp
	 do i=1,nj
	  do j=1,nj
	   phi(i,j)=phiout(i,j)
	  enddo
	 enddo
	endif
       enddo

 4444 continue
      return
      end 
*****************************************************************




*****************************************************************
      subroutine perturbp(phi,phi1,v1)
*****************************************************************
        implicit real*8 (a-h,o-z)
        parameter (nj=3,nv=4,nm=20,nt=57,nc=2)
	real*8 phi(nj,nj),phi1(nj,nj),tempr(1)
	real*8 ev1
	complex*16 eigenp(nj)
	integer ri,ci
	real*8 v1(nj,nj)



	 do i=1,nj
	  do j=1,nj
	   tempr(1)= -1d0+ (DRNUNF()*2d0)
           phi1(i,j)=phi(i,j)+(v1(i,j)*tempr(1))
	  enddo
         enddo

	call DEVLRG(nj,phi1,nj,eigenp)
	ev1 = cdabs(eigenp(1))
        do j= 2,nj
	 ev1= dmax1(ev1,cdabs(eigenp(j)))	
	enddo


	do while(ev1 .ge. 1)
         do i=1,nj
	  do j=1,nj
	   tempr(1)= -1d0+ (DRNUNF()*2d0)
           phi1(i,j)=phi(i,j)+(v1(i,j)*tempr(1))
	  enddo
         enddo
	 call DEVLRG(nj,phi1,nj,eigenp)
	 ev1 = cdabs(eigenp(1))
         do j= 2,nj
	  ev1= dmax1(ev1,cdabs(eigenp(j)))	
	 enddo
	enddo

      return
      end 
*****************************************************************




*****************************************************************
      subroutine updatesig(sigma,ntem,ns,sigout,fsout)
*****************************************************************
	implicit real*8 (a-h,o-z)
	integer ntem,ns,ri,ci,neps,maxeval
        parameter (nj=3,nv=4,nm=20,nt=57,nc=2)
	real*8 tempr(1),dif1(1),dif2(1)
	real*8 sigma(nj,nj),sigma1(nj,nj), sigout(nj,nj)
	real*8 sigvec(nj+(nj*(nj-1)/2)),sig1vec(nj+(nj*(nj-1)/2))
	real*8 fs,fs1,fout,fold
	real*8 v(nj+(nj*(nj-1)/2))
	real*8 p1,p2
	real*8 temp,count(nj+(nj*(nj-1)/2)),ratio
	real*8 ub,lb,diff


	c1=1d0
	eps=0.001d0
	rt=0.25d0
	maxeval=3
	temp=1d0
	diff=0.1d0
	nk=nj+(nj*(nj-1)/2) 
	
	

	do i=1,nj
	 do j=1,nj
	  sigout(i,j)=sigma(i,j)
	 enddo
	enddo
	call  smat2vec(sigma,nj,sigvec) 
	call  fsig(nk,sigvec,fs)
	fout=fs
	
	do ink=1,nk
	 v(ink)= 0.001d0 
	enddo
	fold=fout

	iter=1
        do while(iter .le. maxeval)
	do in=1,ntem
	 count=0d0
	 do is=1,ns
	  call perturbs(sigma,sigma1,v)
	  call  smat2vec(sigma1,nj,sig1vec) 
	  call  fsig(nk,sig1vec,fs1)
          if(fs1<fs) then
	   do i=1,nj
	    do j=1,nj
	     sigma(i,j)=sigma1(i,j)
	    enddo
	   enddo
	   dif1(1)= fs-fs1
	   fs=fs1
	   count=count+1d0
	  endif
          if(fs1<fout) then
	   do i=1,nj
	    do j=1,nj
	     sigout(i,j)=sigma1(i,j)
	    enddo
	   enddo
	   do i=1,nj
	    do j=1,nj
	     sigma(i,j)=sigma1(i,j)
	    enddo
	   enddo
	   dif2(1)= fout-fs1
	   fout=fs1
	   fs=fs1
	  endif
          if(fs1.ge.fs) then
	   p1=DRNUNF()
           p2= dexp((fs-fs1)/temp)
	   if(p2 .gt. p1) then
	    do i=1,nj
	     do j=1,nj
	      sigma(i,j)=sigma1(i,j)
	     enddo
	    enddo
	    fs=fs1
	   endif
	  endif
	 enddo
	 ratio=count(ink)/dfloat(ns)
	 if(ratio .gt. 0.6d0) then
	  do ink=1,nk
	   v(ink)= v(ink)*(1d0+c1*(ratio-0.6d0)/0.4d0)
	    if(v(ink) .gt. diff) then
	     v(ink)= diff
	    endif
	  enddo
	 else if(ratio .lt. 0.4d0) then
	  do ink=1,nk
           v(ink)= v(ink)/(1d0+c1*((0.4d0-ratio)/0.4d0))
	  enddo
	 endif 
	enddo
	iter=iter+1
	 if(dif1(1)<eps .and. dif2(1)<eps) then
	  go to 5555
	 else
	  temp=rt*temp
	  do i=1,nj
	   do j=1,nj
	    sigma(i,j)=sigout(i,j)
	   enddo
	  enddo
	 endif
	enddo
 5555 continue
      return
      end 

*****************************************************************
      subroutine perturbs(sigma,sigma1,v1)
*****************************************************************
        implicit real*8 (a-h,o-z)
        parameter (nj=3,nv=4,nm=20,nt=57,nc=2)
	real*8 sigma(nj,nj),sigma1(nj,nj),tempr(1)
	real*8 eigens(nj),evmin
	integer ri,nk
	real*8 sigvec(nj+(nj*(nj-1)/2)),sigvec1(nj+(nj*(nj-1)/2))
	real*8 v1(nj+(nj*(nj-1)/2) )



	 call  smat2vec(sigma,nj,sigvec)
	 nk=nj+(nj*(nj-1)/2) 


	 do ik=1,nk
	  tempr(1)= -1d0+ (DRNUNF()*2d0)
	  sigvec1(ik)=sigvec(ik)+(v1(ik)*tempr(1))
	 enddo

	 
	 call svec2mat(sigvec1,nj,sigma1)

	do i=1,nj
	 sigma1(i,i)=dabs(sigma1(i,i))
	enddo

	call DEVLSF(nj,sigma1,nj,eigens)
	evmin = eigens(1)
        do j= 2,nj
	 evmin= dmin1(evmin,eigens(j))	
	enddo

	do while(evmin.le.0d0)
	 do ik=1,nk
	  tempr(1)= -1d0+ (DRNUNF()*2d0)
	  sigvec1(ik)=sigvec(ik)+(v1(ik)*tempr(1))
	 enddo
	 call svec2mat(sigvec1,nj,sigma1)
	 do i=1,nj
	  sigma1(i,i)=dabs(sigma1(i,i))
	 enddo
	 call DEVLSF(nj,sigma1,nj,eigens)
	 evmin = eigens(1)
         do j= 2,nj
	  evmin= dmin1(evmin,eigens(j))	
	 enddo
	enddo

      return
      end 

*****************************************************************

*****************************************************************
      subroutine updateta(eta,ntem,ns,etaout,feout)
*****************************************************************
	implicit real*8 (a-h,o-z)
	integer ntem,ns,ri,ci,neps,maxeval
        parameter (nj=3,nv=4,nm=20,nt=57,nc=2)
	real*8 tempr(1),dif1(1),dif2(1)
	real*8 eta(nc,nj),etaout(nc,nj),eta1(nc,nj)
	real*8 etavec(nc*nj),etavec1(nc*nj)
	real*8 fe,feout,fe1,fold
	real*8 v(nc,nj)
	real*8 p1,p2,c1,eps,rt
	real*8 temp,count,ratio,diff


	c1=1d0
	eps=0.001d0
	rt=0.25d0
	maxeval=3
	temp=1d0
	diff=0.1d0

	do inc=1,nc
	 do j=1,nj
	  etaout(inc,j)=eta(inc,j)
	  v(inc,j)=0.001d0 
	 enddo
	enddo

	call  mat2vec(eta,nc,nj,etavec) 
	call  feta(nc*nj,etavec,fe)
	feout=fe
        fold=feout
	iter=1
        do while(iter .le. maxeval)
	 do in=1,ntem
	  do inc=1,nc
	   do j=1,nj
	    count=0d0
	   enddo
	  enddo
	  do is=1,ns
	   call perturbe(eta,eta1,v)
	   call  mat2vec(eta1,nc,nj,etavec1) 
	   call  feta(nc*nj,etavec1,fe1)
           if(fe1 .le. fe) then
	    do inc=1,nc
	     do j=1,nj
              eta(inc,j)=eta1(inc,j)
	     enddo
	    enddo
	    dif1(1)=fe-fe1
	    fe=fe1
	    count=count+1d0
	   endif
           if(fe1 .le. feout) then
	    do inc=1,nc
	     do j=1,nj
              etaout(inc,j)=eta1(inc,j)
	     enddo
	    enddo
	    do inc=1,nc
	     do j=1,nj
              eta(inc,j)=eta1(inc,j)
	     enddo
	    enddo
	    dif2(1)= feout-fe1
	    feout=fe1
	    fe=fe1
	   endif
           if(fe1.ge.fe) then
	    p1=DRNUNF()
            p2= dexp((fe-fe1)/temp)
	    if(p2 .gt. p1) then
	     do inc=1,nc
	      do j=1,nj
               eta(inc,j)=eta1(inc,j)
	      enddo
	     enddo
	     fe=fe1
	    endif
	   endif
	  enddo
	  ratio=count/dfloat(ns)
	  if(ratio .gt. 0.6d0) then
	  do inc=1,nc
	   do j=1,nj
	    v(inc,j)= v(inc,j)*(1d0+c1*(ratio-0.6d0)/0.4d0)
	    if(v(inc,j) .gt. diff) then
	     v(inc,j)= diff
	    endif
	   enddo
	  enddo
	 else if(ratio .lt. 0.4d0) then
	  do inc=1,nc
	   do j=1,nj
	     v(inc,j)= v(inc,j)/(1d0+c1*((0.4d0-ratio)/0.4d0))
	    enddo
	   enddo 
	  endif
	 enddo
	 iter=iter+1
	 if(dif1(1)<eps .and. dif2(1)<eps) then
	  go to 7777
	 else
	  temp=rt*temp
	  do inc=1,nc
	   do j=1,nj
	    eta(inc,j)=etaout(inc,j)
	   enddo
	  enddo
	 endif
	enddo

 7777 continue
      return
      end 
*****************************************************************



*****************************************************************
      subroutine perturbe(eta,eta1,v1)
*****************************************************************
        implicit real*8 (a-h,o-z)
        parameter (nj=3,nv=4,nm=20,nt=57,nc=2)
	real*8 eta(nc,nj),eta1(nc,nj),tempr(1)
	integer ri,ci
	real*8 v1(nc,nj)

	 do ic=1,nc
	  do j=1,nj
	   tempr(1)= -1d0+ (DRNUNF()*2d0)
           eta1(ic,j)= eta(ic,j)+(v1(ic,j)*tempr(1))
	  enddo
         enddo
	  
      return
      end 
*****************************************************************



*****************************************************************
	subroutine updatedhatd(dhatd,ntem,ns,dhatdout,fdout)
*****************************************************************
	implicit real*8 (a-h,o-z)
	integer ntem,ns,ri,ci,neps,maxeval
        parameter (nj=3,nv=4,nm=20,nt=57,nc=2)
	real*8 tempr(1),dif1(1),dif2(1)
	real*8 dhatd(nt,nj),dhatdout(nt,nj),dhatd1(nt,nj)
	real*8 dhatdvec(nt*nj)
	real*8 fdhr,fdout,fdhr1,fold
	real*8 v(nt,nj)
	real*8 p1,p2,c1,eps,rt
	real*8 temp,count,ratio,diff


	c1=1d0
	eps=0.001d0
	rt=0.25d0
	maxeval=3
	temp=1d0
	diff=1d0


	do it=1,nt
	 do j=1,nj
	  dhatdout(it,j)= dhatd(it,j)
	  v(it,j)=0.1d0
	 enddo
	enddo

	call  mat2vec(dhatd,nt,nj, dhatdvec) 
	call  fdhatd(nt*nj,dhatdvec,fdhr)
	fdout=fdhr

	fold= fdout
        iter=1
        do while(iter .le. maxeval)
	 do in=1,ntem
	  count=0d0
	  do is=1,ns
	   call perturbd(dhatd,dhatd1,v)
	   call  mat2vec(dhatd1,nt,nj,dhatdvec) 
           call fdhatd(nt*nj,dhatdvec,fdhr1)
           if(fdhr1<fdhr) then
	    do it=1,nt
	     do j=1,nj
	      dhatd(it,j)= dhatd1(it,j)
	     enddo
	    enddo	    
	    dif1(1)= fdhr-fdhr1
	    fdhr=fdhr1
	    count=count+1d0
	   endif
           if(fdhr1<fdout) then
	    do it=1,nt
	     do j=1,nj
	      dhatdout(it,j)= dhatd1(it,j)
	     enddo
	    enddo
	    do it=1,nt
	     do j=1,nj
	      dhatd(it,j)= dhatd1(it,j)
	     enddo
	    enddo	    
	    dif2(1)= fdout-fdhr1
	    fdout=fdhr1
	    fdhr=fdhr1
	   endif
           if(fdhr1 .ge. fdhr) then
	    p1=DRNUNF()
            p2= dexp((fdhr-fdhr1)/temp)
	    if(p2 .gt. p1) then
	     do it=1,nt
	      do j=1,nj
	       dhatd(it,j)= dhatd1(it,j)
	      enddo
	     enddo
	     fdhr=fdhr1
	    endif
	   endif
	  enddo
	  ratio=count/dfloat(ns)
	  if(ratio .gt. 0.6d0) then
	   do it=1,nt
	    do j=1,nj
	     v(it,j)= v(it,j)*(1d0+c1*(ratio-0.6d0)/0.4d0)
	    enddo
	   enddo
	  else if(ratio .lt. 0.4d0) then
	   do it=1,nt
	    do j=1,nj
	     v(it,j)= v(it,j)/(1d0+c1*((0.4d0-ratio)/0.4d0))
	     if(v(it,j) .gt. diff) then
	      v(it,j)= diff
	     endif
	    enddo
	   enddo
	  endif
	 enddo
	 iter=iter+1
	 if(dif1(1)<eps .and. dif2(1)<eps) then
	  go to 8888
	 else
	  temp=rt*temp
	 do it=1,nt
	  do j=1,nj
	   dhatd(it,j)= dhatdout(it,j)
	  enddo
	 enddo
	endif
       enddo

 8888 continue
      return
      end 
*****************************************************************


*****************************************************************
      subroutine perturbd(dhatd,dhatd1,v1)
*****************************************************************
        implicit real*8 (a-h,o-z)
        parameter (nj=3,nv=4,nm=20,nt=57,nc=2)
	real*8 dhatd1(nt,nj), dhatd(nt,nj),tempr(1)
	integer ri,ci
	real*8 v1(nt,nj)

	do it=1,nt
	 do j=1,nj
	  tempr(1)= -1d0+ (DRNUNF()*2d0)
	  dhatd1(it,j)=dhatd(it,j)+(v1(it, j)*tempr(1))
	 enddo
	enddo

      return
      end 
*****************************************************************




*****************************************************************
	subroutine whichmax(a,n,ind)
*****************************************************************

	integer n,ind
	real*8 a(n),maxval

	ind=1
	maxval=a(1)
	do i=2,n
	 if(a(i) .gt. maxval) then
	  ind=i
	  maxval=a(i)
	 endif 
	enddo

      return
      end 
*****************************************************************











