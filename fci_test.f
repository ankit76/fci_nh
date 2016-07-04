CABHVFCI.  A DETERMINANT BASED FULL CONFIGURATION INTERACTION PROGRAM.
C1   P.J. KNOWLES, N.C. HANDY
CREF. IN COMP. PHYS. COMMUN. 54 (1989) 75
CFile: fci.f
      PROGRAM FCI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      COMMON /CHAMIL/ ISYM,NNCI
      COMMON /TAPES / INP,IOUT
      DIMENSION MULTT(8,8)
      DATA MULTT /1,2,3,4,5,6,7,8, 2,1,4,3,6,5,8,7,
     2            3,4,1,2,7,8,5,6, 4,3,2,1,8,7,6,5,
     4            5,6,7,8,1,2,3,4, 6,5,8,7,2,1,4,3,
     6            7,8,5,6,3,4,1,2, 8,7,6,5,4,3,2,1/
      INP = 5
      IOUT = 6
      WRITE (IOUT,35)
35    FORMAT(' PROGRAM * FCI  (Determinant based Full CI)',5X,
     >'Author: P.J. Knowles, 1984'/)
C...  SET UP DEFAULT OR INITIAL VALUES
      MEMORY = 10000000
      IPRTIM = 2
      IPRINT = -1
      IPRDIA = 2
      IPRHAM = 2
      IOPTIO = 6
      MAXIT = 25
      THR = 1D-5
      THRRES = 0.05
      NROOT = 1
      NACT = 0
      NELEC = 0
      MS2 = 0
      ISYM = 1
      DO 2 I=1,8
      NOCC(I) = 0
C...  D2H MULTIPLICATION TABLE
      DO 2 J=1,8
2     MULT(J,I) = MULTT(J,I)
      DO 3 I=1,256
3     ITYPEA(I)=1
C
C...  READ INPUT
      CALL INPDAT
C
C...  INITIALIZE MEMORY HANDLING
      I = INICOR(MEMORY)
C
C...  NUMBERS OF ALPHA, BETA ELECTRONS
      NA = (NELEC+MS2)/2
      NB = (NELEC-MS2)/2
      DO 36 I=1,NACT
36    NOCC(ITYPEA(I)) = NOCC(ITYPEA(I))+1
      CALL INITC
C
      WRITE (IOUT,25) ISYM
25    FORMAT(/' Space symmetry:',T30,I3)
      WRITE (IOUT,21) MAXIT,THR,THRRES,NROOT
21    FORMAT(
     4' Maximum iterations:',   T30,I3/
     5' Convergence threshold:',T30,F12.7/
     8' Output threshold:',     T30,F12.7/
     A' Number of roots:',      T30,I3)
      NNCI = NCI(ISYM)
      IV1 = ICORR(NNCI)
      IV2 = ICORR(NNCI)
      CALL DAVIDS (Q(IV1),Q(IV2),IOPTIO)
      CALL CORLSR(IV1)
C
      CALL ENDCOR
      END
      SUBROUTINE FCMXM (D,Z,E,NK,NIJ,NOPCNT)
CSUBR MATRIX MULTIPLICATION KERNEL FOR FULL CI PROGRAM
CSUBR  E = D * Z
CSUBR  Z IS INTEGRAL MATRIX AND NOT SPARSE, DIMENSION NIJ=NUMBER OF
CSUBR            ORBITAL PAIRS OF GIVEN SYMMETRY
CSUBR  D,E LEADING DIMENSION NK=NUMBER OF DETERMINANTS OF A SYMMETRY
CSUBR            BLOCK, OR SUB-BATCH OF THIS
CSUBR  D IS VERY SPARSE IN 1ST ITERATION OF CI, GRADUALLY RISES TO
CSUBR             USUALLY 50-80% POPULATED
CSUBR OPTIMUM FORM OF THIS ROUTINE IS DEPENDENT ON MACHINE ARCHITECTURE
CSUBR USUALLY NK>NIJ, SO BETTER TO VECTORISE ALONG NK, BUT SOMETIMES
CSUBR BETTER TO EXPLOIT SPARSITY OF D
CSUBR ROUTINE SHOULD RETURN NOPCNT -- NUMBER OF FLOATING POINT OPS --
CSUBR IF ADDITIONAL TIMING INFORMATION IS REQUIRED
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NHALF=30)
      DIMENSION D(NK,NIJ),Z(NIJ,NIJ),E(NK,NIJ)
      DO 1 J=1,NK*NIJ
1     E(J,1) = 0D0
C...  EVALUATE SPARSITY OF D
      IDS = 0
      DO 7 J=1,NK*NIJ
7     IF (D(J,1).NE.0D0) IDS=IDS+1
      IZS = NIJ**2
      IF (IDS*(NIJ+NHALF).LT.IZS*(NK+NHALF)) THEN
C...    USE SPARSITY OF D
        NOPCNT = NOPCNT + 2*IDS*NIJ
C        CALL MXMA (Z,1,NIJ, D,NK,1, E,NK,1, NIJ,NIJ,NK)
        DO 2 K=1,NK
        DO 2 J=1,NIJ
        IF (D(K,J).EQ.0D0) GOTO 2
        DO 21 I=1,NIJ
21      E(K,I) = E(K,I) + D(K,J)*Z(J,I)
2       CONTINUE
      ELSE
C...    VECTORISE ALONG D
        NOPCNT = NOPCNT + 2*IZS*NK
C        CALL MXMA (D,1,NK, Z,NIJ,1, E,1,NK, NK,NIJ,NIJ)
        DO 3 I=1,NIJ
        DO 3 J=1,NIJ
C        IF (Z(J,I).EQ.0D0) GOTO 3
        DO 31 K=1,NK
31      E(K,I) = E(K,I) + D(K,J)*Z(J,I)
3       CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE INPDAT
CSUBR DATA INPUT FOR PORTABLE FULL CI PROGRAM
CSUBR NAMELIST INPUT MAY REQUIRE ADJUSTMENT ON SOME MACHINES
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      COMMON /CHAMIL/ ISYM,NNCI
      COMMON /TAPES / INP,IOUT
      INTEGER*4 ORBSYM(255)
c
c ONEFIL: Name of the one-electron integral file
c TWOFIL: Name of the two-electron integral file
c
      character*72 onefil,twofil
      common/filnam/onefil,twofil
c
      EQUIVALENCE (ORBSYM(1),ITYPEA(1))
      NAMELIST /FCI/ NORB,NELEC,MS2,ISYM,ORBSYM
     >                 ,IPRINT,IPRDIA,IPRHAM,IPRTIM
     >                         ,INT,MEMORY,CORE,MULT
     >                 ,MAXIT,THR,THRRES,NROOT,
     >                 onefil,twofil 
      INT = -1
      READ (INP,FCI)
      IF (INT.NE.-1) INP=INT
      NACT = NORB
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE INPINT (Z1,Z2)
CSUBR INTEGRAL INPUT FOR PORTABLE FULL CI PROGRAM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      COMMON /TAPES / INP,IOUT
      DIMENSION Z2(*),Z1(*)
c
      character*72 onefil,twofil
      common/filnam/onefil,twofil
      character*72 filen
c
      logical exst
c
      CALL FZERO (Z1,NPAIR(1))
      CALL FZERO (Z2,NINT2)
c
c One-electron integral related
c
       inquire(file=onefil,exist=exst)
       if(.not.exst)then
         write(iout,*)'sr inpint: 1-elec. integral file not found'
         write(iout,*)'file name',onefil
         stop
       end if
c
       iunit=100
       filen=onefil
       open(unit=iunit,file=filen,status='old',form='unformatted',err=1)
       goto 2
c
 1     write(iout,*)'sr inpint: error opening file',filen
       stop
c
 2     continue
       rewind(iunit)
       read(iunit)nxorb,n1x,nzero
       if(nxorb.ne.nact)then
         write(iout,*)'sr inpint: nxorb.ne.nact in 1-elec file'
         write(iout,*)'nxorb,nact',nxorb,nact
         stop
       end if
c
       do 3 ixx=1,nzero
         read(iunit)ijorb,z
         call unpack(ijorb,i,j)
         if(max(i,j).gt.nact)then
           write(iout,*)'sr inpint: orbital indices out of range'
           write(iout,*)'i,j',i,j
           stop
         end if
         IJ=IQ(INTADT+I-1+(J-1)*NACT)
         Z1(IJ)=Z
 3     continue
c
c Read the core energy as well
c
       read(iunit)core
c
       close(unit=iunit,status='keep')
c
c Two-electron integral related
c
       inquire(file=twofil,exist=exst)
       if(.not.exst)then
         write(iout,*)'sr inpint: 2-elec. integral file not found'
         write(iout,*)'file name',twofil
         stop
       end if
c
       filen=twofil
       open(unit=iunit,file=filen,status='old',form='unformatted',
     : err=1)
       rewind(iunit)
 4     continue
       read(iunit,end=6,err=7)ijorb,klorb,z
       call unpack(ijorb,i,j)
       call unpack(klorb,k,l)
c       write(iout,*)'i,j,k,l,z',i,j,k,l,z
c      if(i.ge.25.and.j.ge.25) write(iout,*)'i,j,k,l,z',i,j,k,l,z       
       if(max(i,j,k,l).gt.nact)then
         write(iout,*)'sr inpint: orbital indices out of range'
         write(iout,*)'i,j,k,l',i,j,k,l
         stop
       end if
       isymij = mult(itypea(i),itypea(j))
       ij=iq(intadt+i-1+(j-1)*nact)
       kl=iq(intadt+k-1+(l-1)*nact)
       z2(ij+(kl-1)*npair(isymij)+intoff(isymij))=z
       z2(kl+(ij-1)*npair(isymij)+intoff(isymij))=z
c           
       goto 4
c
 6     write(iout,*)'sr inpint: done reading integrals'
       close(iunit,status='keep')
       return
c
 7     write(iout,*)'sr inpint: error reading 2-el file'
       stop
c
       end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE INPINT2 (Z1,Z2,BUF,IBUF,NBUF)
cSUBR INTEGRAL INPUT FOR PORTABLE FULL CI PROGRAM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      COMMON /TAPES / INP,IOUT
      DIMENSION Z2(*),Z1(*)
      DIMENSION BUF(NBUF),IBUF(NBUF)

      character*72 onefil,twofil
      common/filnam/onefil,twofil
      character*72 filen
c
      logical exst
c
      CALL FZERO (Z1,NPAIR(1))
      CALL FZERO (Z2,NINT2)
c
c One-electron integral related
c
       inquire(file=onefil,exist=exst)
       if(.not.exst)then
         write(iout,*)'sr inpint: 1-elec. integral file not found'
         write(iout,*)'file name',onefil
         stop
       end if
c
       iunit=100
       filen=onefil
       open(unit=iunit,file=filen,status='old',form='unformatted',
     : err=1)
       goto 2
c
 1     write(iout,*)'sr inpint: error opening file',filen
       stop
c
 2     continue
       rewind(iunit)
       read(iunit)nxorb,n1x,nzero
       if(nxorb.ne.nact)then
         write(iout,*)'sr inpint: nxorb.ne.nact in 1-elec file'
         write(iout,*)'nxorb,nact',nxorb,nact
         stop
       end if
c
       read(iunit)(ibuf(i),i=1,nzero)
       read(iunit)(buf(i),i=1,nzero)
       do 3 ixx=1,nzero
         call unpack(ibuf(ixx),i,j)
         ij=iq(intadt+i-1+(j-1)*nact)
         z1(ij)=buf(ixx)
 3     continue
c
       close(unit=iunit,status='keep')
c
c Two-electron integral related
c
       inquire(file=twofil,exist=exst)
       if(.not.exst)then
         write(iout,*)'sr inpint: 2-elec. integral file not found'
         write(iout,*)'file name',twofil
         stop
       end if
c
       filen=twofil
       open(unit=iunit,file=filen,status='old',form='unformatted',
     : err=1)
       rewind(iunit)
c
 4     continue
       istage=0
       read(iunit,end=6,err=7)nzero,ijorb
       call unpack(ijorb,i,j)
       istage=1
       read(iunit,err=7)(ibuf(ixx),ixx=1,nzero)
       istage=2
       read(iunit,err=7)(buf(ixx),ixx=1,nzero)
       do 5 ixx=1,nzero
         call unpack(ibuf(ixx),k,l)
         if(max(i,j,k,l).gt.nact)then
           write(iout,*)'sr inpint: 2-el orbital indices out of bound'
           write(iout,*)'norb,i,j,k,l',nact,i,j,k,l
           isymij= mult(itypea(i),itypea(j))
           ij=iq(intadt+i-1+(j-1)*nact)
           kl=iq(intadt+k-1+(l-1)*nact)
           z2(ij+(kl-1)*npair(isymij)+intoff(isymij))=buf(ixx)
           z2(kl+(ij-1)*npair(isymij)+intoff(isymij))=buf(ixx)
         end if
 5     continue
c           
       goto 4
c
 6     close(iunit,status='keep')
       return
c
 7     write(iout,*)'sr inpint: error reading 2-el file, istage=',
     : istage
       stop
c
       end 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE DAVIDS (V1,V2,IOPTIO)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /TAPES / INP,IOUT
C
C...  DAVIDSON DIAGONALISER, SINGLE ROOT
C...  WARNING, NOT NECESSARILY SAFE WITH DEGENERACIES, ETC.
      PARAMETER (IFILD=92,IFIL1=93,IFIL2=94,MAXITR=1000)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /CHAMIL/ ISYMH,N
      DIMENSION V1(N),V2(N)
      DIMENSION ZM(MAXITR,MAXITR),ZV(MAXITR,MAXITR),R(MAXITR)
C
      TRIAL = 0
      RESULT = 0
      MAXITS = MIN(MAXIT,MAXITR,N)
      CALL FCDIAG (V2)
      IF (IOPTIO.GE.5) CALL OUTVEC (V2,N,'DIAGONAL ELEMENTS')
C
      OPEN (unit=IFILD,file='filed.tmp',STATUS='unknown',
     : FORM='UNFORMATTED')
      OPEN (unit=IFIL1,file='file1.tmp',STATUS='unknown',
     : FORM='UNFORMATTED')
      OPEN (unit=IFIL2,file='file2.tmp',STATUS='unknown',
     : FORM='UNFORMATTED')
      REWIND IFILD
      REWIND IFIL1
      REWIND IFIL2
C
      WRITE (IFILD) V2
      REWIND IFILD
      CALL TRIALV (V1,V2)
C
      WRITE (IOUT,'(/'' It  Tr    CPU'',2X,''Convergence'',3X,
     >''Energy''/)')
      DO 190 IT=1,MAXITS
      WRITE (IFIL1) V1
C...  CALL THE CI PROGRAM V2 = H . V1  + V2
      DO 30 I=1,N
30    V2(I)=0.0D0
      CALL FSIGMA (V1,V2)
      WRITE (IFIL2) V2
C
      REWIND IFIL1
      DO 50 JT=1,IT
      READ (IFIL1) V1
      ZM(JT,IT) = 0.0D0
      DO 40 I=1,N
40    ZM(JT,IT) = ZM(JT,IT) + V1(I)*V2(I)
50    ZM(IT,JT) = ZM(JT,IT)
      REWIND IFIL1
      REWIND IFIL2
C
C...  SMALL MATRIX DIAGONALISER
      DO 60 JT=1,IT
      DO 60 KT=1,IT
60    ZV(KT,JT) = ZM(KT,JT)
      CALL DIAG2 (MAXITR,IT,R,ZV)
C
      TEST = 0
      DO 61 IROOT=1,MIN(IT,NROOT)
      TES = ABS(ZV(IT,IROOT))
      IF (TES.LT.TEST) GOTO 61
      ITRACK = IROOT
      TEST = TES
61    CONTINUE
      IF (IT.EQ.1) EREF=R(1)
      WRITE (IOUT,70) IT,ITRACK,SECOND(),TEST
     >,(R(IROOT)+CORE,IROOT=1,MIN(NROOT,IT))
70    FORMAT(I3,I4,F7.1,F12.8,5F17.9)
      IF (TEST.LT.THR) GOTO 210
C
C...  MAKE RESIDUAL
      DO 80 I=1,N
80    V1(I) = 0.0D0
      DO 100 JT=1,IT
      READ (IFIL1) V2
      DO 90 I=1,N
90    V1(I) = V1(I) + (-R(ITRACK)*ZV(JT,ITRACK))*V2(I)
      READ (IFIL2) V2
      DO 100 I=1,N
100   V1(I) = V1(I) + ZV(JT,ITRACK) * V2(I)
      READ (IFILD) V2
      IF (IT.EQ.1) THEN
      DO 110 JJ = 1,N
      IF (ABS(V2(JJ)-EREF).LT.1D-9) V2(JJ)=1D20
110   CONTINUE
      END IF
      REWIND IFILD
      DO 130 I=1,N
130   V1(I) = V1(I) / (V2(I)-R(1)+1.0D-10)
C
C...  ORTHOGONALISE
140   REWIND IFIL1
      DO 160 JT=1,IT
      READ (IFIL1) V2
      ZZ = 0.0D0
      DO 150 I=1,N
150   ZZ = ZZ + V1(I)*V2(I)
      DO 160 I=1,N
160   V1(I) = V1(I) + (-ZZ) * V2(I)
      ZZ = 0.0D0
      DO 170 I=1,N
170   ZZ = ZZ + V1(I)**2
      IF (ZZ.LT.1D-20) GOTO 210
      ZZ = 1.0 / SQRT(ZZ)
      DO 180 I=1,N
180   V1(I) = V1(I) * ZZ
C...  REORTHOGONALISE IF A LOT HAS BEEN ANNIHILATED
      IF (ZZ.GT.1D3) GOTO 140
C
190   CONTINUE
      IT = MAXITS
      WRITE (IOUT,200)
200   FORMAT(/1X,'*** Convergence not achieved in max iterations')
C
210   DO 280 IROOT=1,NROOT
      REWIND IFIL1
      DO 220 I=1,N
220   V1(I) = 0.0D0
      DO 230 JT=1,IT
      READ (IFIL1) V2
      DO 230 I=1,N
230   V1(I) = V1(I) + ZV(JT,IROOT) * V2(I)
      WRITE (IOUT,'('' State'',I2,5X,''Energy'',F20.12)')
     > IROOT,R(IROOT)+CORE
      IF (IROOT.EQ.1) WRITE (IOUT,'('' Correlation energy'',F20.12)')
     > R(1)-EREF
      WRITE (IOUT,'(/'' Final CI vector''/)')
      CALL FCVECO (V1,ISYMH,THRRES)
280   CONTINUE
C
      CLOSE (IFILD,STATUS='DELETE')
      CLOSE (IFIL1,STATUS='DELETE')
      CLOSE (IFIL2,STATUS='DELETE')
C
      RETURN
      END
      SUBROUTINE TRIALV (V1,V2)
CSUBR SUBROUTINE TO GENERATE TRIAL CI VECTOR
CSUBR COULD BE REPLACED BY E.G. READING VECTOR IF REQUIRED
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /TAPES / INP,IOUT
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      COMMON /CHAMIL/ ISYMH,N
      DIMENSION IMIN(64)
      DIMENSION V1(N),V2(N)
      TRIAL = 0
      IF (NINT(TRIAL).EQ.0) THEN
C...  FIND THE DETERMINANT WITH LOWEST ENERGY
      VMIN = 1D20
      DO 10 I=1,N
      IF (V2(I).GT.VMIN) GOTO 10
      IMIN(1) = I
      VMIN = V2(I)
10    CONTINUE
C...  CONSTRUCT A SPIN EIGENFUNCTION CONTAINING THIS DETERMINANT
      CALL FCSPAD (V1,IMIN(1))
      WRITE (IOUT,20)
20    FORMAT(/1X,'Initial configuration generated:')
      NREF = 0
      DO 40 I=1,N
      IF (V1(I).EQ.0.0D0) GOTO 40
      WRITE (IOUT,30) I,V1(I),V2(I)+CORE
      NREF = NREF+1
      IMIN(NREF) = I
30    FORMAT(I8,2F15.7)
40    CONTINUE
      ELSE
C... READ FROM FILE
      CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE FCVECO (V,ISYM,THRESH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /TAPES / INP,IOUT
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      CHARACTER*130 TEST
      DIMENSION V(*)
      COMMON /CLOCAL/ ICGA(32),ICGB(32)
      N = NCI(ISYM)
      IBASE = ICORR(0)
      INTERA = ICORI(NA*NACT)
      CALL INITI (IQ(INTERA),NA,NACT)
      INTERB = ICORI(NB*NACT)
      CALL INITI (IQ(INTERB),NA,NACT)
      IC = ICORI(NSTRAA)
      ITZ = ICORI(NSTRBB)
      CALL FCMIC (IQ(IC),IQ(ITZ),ISYM)
      DO 60 II=1,NCI(ISYM)
      IF (ABS(V(II)).LT.THRESH) GOTO 60
      CALL FCSTRG (II,ISYM,ICGA,ICGB)
      WRITE (IOUT,'(1X,F20.12,64I3)') V(II)
     1,(ICGA(I)+NFROZ,I=1,NA)
     2,(ICGB(I)+NFROZ,I=1,NB)
60    CONTINUE
      WRITE (IOUT,'(1X,''/EOF'')')
      CALL CORLSR(IBASE)
      RETURN
      END
      SUBROUTINE INITC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /TAPES / INP,IOUT
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      DIMENSION ICG(32)
C...  LENGTHS OF CI EXPANSIONS
C...  NSTRA, NSTRB ARE NUMBERS OF ALPHA,BETA STRINGS OF EACH SYMMETRY
C...  NCI ARE NUMBERS OF DETERMINANTS OF EACH SYMMETRY
C...  NPAIR ARE NUMBERS OF ORBITAL PAIRS OF EACH SYMMETRY
C.... INTOFF OFFSET FOR SYMMETRY BLOCKS OF TWO ELECTRON INTEGRALS
      NSTRAA = 0
      NSTRBB = 0
      MAXAA  = 0
      MAXBB  = 0
      DO 100 ISYM=1,8
      ICG(1) = 0
      NSTRA(ISYM) = -1
80    NSTRA(ISYM) = NSTRA(ISYM)+1
      CALL FCSTRS (NA,ICG,ISYM,ITYPEA)
      IF (ICG(1).NE.0) GOTO 80
      NSTRAA = NSTRAA + NSTRA(ISYM)
      MAXAA = MAX(MAXAA,NSTRA(ISYM))
      ICG(1) = 0
      NSTRB(ISYM) = -1
90    NSTRB(ISYM) = NSTRB(ISYM)+1
      CALL FCSTRS (NB,ICG,ISYM,ITYPEA)
      IF (ICG(1).NE.0) GOTO 90
      NSTRBB = NSTRBB + NSTRB(ISYM)
      MAXBB = MAX(MAXBB,NSTRB(ISYM))
100   CONTINUE
      DO 110 ISYM=1,8
      NCI(ISYM) = 0
      DO 110 ISYMB=1,8
110   NCI(ISYM) = NCI(ISYM) + NSTRA(MULT(ISYMB,ISYM))*NSTRB(ISYMB)
C
C...  ORBITAL PAIR ADDRESSING, ETC..
      DO 120 I=1,256
120   ICF(I) = (I*(I-1))/2
      INTADT = ICORI(NACT*NACT)
      NINT2=0
      MAXPAR = 0
      DO 140 ISYM=1,8
      INTOFF(ISYM) = NINT2
      NPAIR(ISYM)=0
      NPB = 0
      DO 130 I=1,NACT
      DO 130 J=1,I
      IF (ISYM.EQ.MULT(ITYPEA(I),ITYPEA(J))) THEN
      NPAIR(ISYM) = NPAIR(ISYM)+1
      IQ(INTADT+(I-1)*NACT+(J-1)) = NPAIR(ISYM)
      IQ(INTADT+(I-1)+NACT*(J-1)) = NPAIR(ISYM)
      END IF
130   CONTINUE
      MAXPAR = MAX(NPAIR(ISYM),MAXPAR)
140   NINT2 = NINT2 + NPAIR(ISYM)**2
C...  UPPER BOUND ON NUMBER OF SINGLE REPLACEMENTS FROM ANY STRING
      MAXREP = MIN(MAXPAR,MAX(NA*(NACT+1-NA),NB*(NACT+1-NB)))
C
      WRITE (IOUT,141)
     > NACT,(NOCC(I),I=1,8)
     >,NA+NB,(NA-NB)*0.5,NPAIR,NSTRA,NSTRB,NCI
141   FORMAT(
     1' Active orbitals:',              T30,I3,'(',8I2,')'/
     2' Active electrons:',             T30,I3/
     3' Spin quantum number:',          T30,F5.1/
     4' Orbital pairs:',                T30,8I6/
     5' Strings:',                      T30,8I6/T30,8I6/
     6' Determinants:',                 T30,4I12/T30,4I12)
C      WRITE (IOUT,'('' ACTIVE ORBITAL SYMMETRIES:'')')
C      WRITE (IOUT,150) (ITYPEA(I),I=1,NACT)
C150   FORMAT(1X,40I2)
C
C...  LOAD INTEGRALS
      INTAA = ICORR(NINT2)
      INTHAA = ICORR(NACT)
      INTJAA = ICORR(NACT**2)
      INTKAA = ICORR(NACT**2)
      IBASE = ICORR(0)
      IONEA = ICORR(NPAIR(1))
      CALL INPINT (Q(IONEA),Q(INTAA))
C
C...  INTEGRALS FOR DIAGONAL ELEMENTS OF HAMILTONIAN
      DO 10 I=1,NACT
10    Q(INTHAA-1+I) = Q(IONEA-1+IQ(INTADT+(I-1)*(NACT+1)))
      DO 11 I=1,NACT
      DO 11 J=1,NACT
      IJ = IQ(INTADT+(I-1)*NACT+J-1)
      II = IQ(INTADT+(I-1)*NACT+I-1)
      JJ = IQ(INTADT+(J-1)*NACT+J-1)
      ISYMIJ = MULT(ITYPEA(I),ITYPEA(J))
      Q(INTJAA+(I-1)*NACT+J-1) = Q(INTAA+(II-1)*NPAIR(1)+JJ-1)
11    Q(INTKAA+(I-1)*NACT+J-1) =
     >     Q(INTAA+INTOFF(ISYMIJ)+(IJ-1)*(NPAIR(ISYMIJ)+1))
C
C...  MODIFY INTEGRALS SUCH THAT HAMILTONIAN BECOMES
C     SUM(IJKL)  E(IJ) E(KL) (IJ|KL)"
C     (NO 1 ELECTRON TERM REMAINS)
C
C...  (IJ|KL)' = 0.5 * (IJ|KL)
      DO 20 IJKL=1,NINT2
20    Q(INTAA-1+IJKL) = Q(INTAA-1+IJKL)*0.5D0
C
C...  H'(IJ) = H(IJ) - SUM(K) (IK|KJ)'
      IJ = IONEA-1
      DO 40 I=1,NACT
      DO 40 J=1,I
      IF (ITYPEA(I).NE.ITYPEA(J)) GOTO 40
      IJ = IJ+1
      DO 30 K=1,NACT
      ISYMIK = MULT(ITYPEA(I),ITYPEA(K))
      IK =IQ(INTADT+(I-1)*NACT+K-1)-1
      JK =IQ(INTADT+(J-1)*NACT+K-1)-1
30    Q(IJ) = Q(IJ) - Q(INTAA+INTOFF(ISYMIK)+IK*NPAIR(ISYMIK)+JK)
40    CONTINUE
C
C...  (IJ|KL)" = (IJ|KL)'+(1/2*NELEC)*(DEL(IJ)*H'(KL)+DEL(KL)*H'(IJ))
      FAC = 2*(NA+NB)
      FAC = 1D0/FAC
      DO 70 K=1,NACT
      KK = IQ(INTADT+(K-1)*(NACT+1))
      DO 50 IJ=1,NPAIR(1)
50    Q(INTAA+(KK-1)*NPAIR(1)+IJ-1) = Q(INTAA+(KK-1)*NPAIR(1)+IJ-1)
     > + FAC*Q(IONEA-1+IJ)
      DO 60 IJ=1,NPAIR(1)
60    Q(INTAA+(IJ-1)*NPAIR(1)+KK-1) = Q(INTAA+(IJ-1)*NPAIR(1)+KK-1)
     > + FAC*Q(IONEA-1+IJ)
70    CONTINUE
      CALL CORLSR(IBASE)
      WRITE (IOUT,'('' Core energy:'',T30,F20.12)') CORE
      RETURN
      END
      SUBROUTINE FSIGMA (C,S)
CSUBR RETURNS S = H * C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /TAPES / INP,IOUT
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      LOGICAL TIMING
      DIMENSION IOFFCI(8),C(*),S(*)
      COMMON /CHAMIL/ ISYMH,NNCI
      TIMING = IPRTIM.GE.0
       NOPMXM = 0
       NOPMX  = 0
       NOPOTH = 0
       NOPCCA = 0
       NOPCCB = 0
       NLOCCA = 0
       NLOCCB = 0
       NOPONE = 0
       T0 = SECOND()
       TMXM = 0
       TCCA = 0
       TCCB = 0
       TONE = 0
      IBASE = ICORR(0)
      NNCI = 0
      DO 10 ISYMB=1,8
      ISYMA = MULT(ISYMH,ISYMB)
      IOFFCI(ISYMB) = NNCI
      IF (IPRHAM.GT.0) WRITE (IOUT,*)'ISYMB,ISYMA,NNCI,NSTRA,NSTRB '
     1,ISYMB,ISYMA,NNCI,NSTRA(ISYMA),NSTRB(ISYMB)
10    NNCI = NNCI + NSTRA(ISYMA)*NSTRB(ISYMB)
      IF (IPRHAM.GT.0) CALL OUTVEC (C,NNCI,'c passed to sigma')
      IC = ICORI(NSTRAA)
      ITZ = ICORI(NSTRBB)
      CALL FCMIC (IQ(IC),IQ(ITZ),ISYMH)
      DO 170 ISYMA=1,8
      DO 170 ISYMB=1,8
      IBASES = ICORR(0)
      ISYMR = MULT(ISYMH,MULT(ISYMA,ISYMB))
      ISYMBE = MULT(ISYMB,ISYMR)
      NSBE = NSTRB(ISYMBE)
      NSA = NSTRA(ISYMA)
      NSB = NSTRB(ISYMB)
      NPAIRR = NPAIR(ISYMR)
      IF (NSA.EQ.0.OR.NSB.EQ.0.OR.NPAIRR.EQ.0) GOTO 170
      MAXREA = MIN(NA*(NACT-NA+1),NPAIRR)
      MAXREB = MIN(NB*(NACT-NB+1),NPAIRR)
      MAXR2A = MAXREA * 2
      MAXR3A = MAXREA * 3
      MAXR2B = MAXREB * 2
      MAXR3B = MAXREB * 3
      IBB = ICORR(0)
      IWA = ICORI(MAXR3A*NSA)
      NWA = ICORI(NSA)
      IWB = ICORI(MAXR3B*NSB)
      NWB = ICORI(NSB)
C...  CONSTRUCT ALL SINGLE REPLACEMENTS
      T1 = SECOND()
      INTER = ICORI(NA*NACT)
      CALL INITI (IQ(INTER),NA,NACT)
      CALL ONELAL (NA,ITYPEA,IQ(INTER),IQ(INTADT),ISYMA,ISYMR
     1 ,IQ(IWA),IQ(NWA),IQ(IC),MAXREA)
      CALL CORLSI(INTER)
      INTER = ICORI(NB*NACT)
      CALL INITI (IQ(INTER),NB,NACT)
      CALL ONELAL (NB,ITYPEA,IQ(INTER),IQ(INTADT),ISYMB,ISYMR
     1 ,IQ(IWB),IQ(NWB),IQ(ITZ),MAXREB)
      TONE = TONE+SECOND()-T1
      CALL CORLSI(INTER)
C...  DETERMINE BLOCKING
      NMAT = 2
      MWORD = ICORRM()/NMAT
      WORDS = MWORD
      NWORD = SQRT(WORDS/NPAIRR)
      NABLK = (MIN(NWORD,NSA)/2)*2+1
      DO 20 NBBLK=(NSB/2)*2+1,1,-2
      IF (NABLK*NBBLK*NPAIRR.LE.MWORD) GOTO 30
20    CONTINUE
30    IDA = ICORR(NABLK*NBBLK*NPAIRR)
      IEA = ICORR(NABLK*NBBLK*NPAIRR)
      IDB = IDA
      IEB = IEA
C...  LOOP OVER BLOCKS
      DO 160 IOFFA=0,NSA-1,NABLK
      NAA = MIN(NABLK,NSA-IOFFA)
      DO 160 IOFFB=0,NSB-1,NBBLK
      NBB = MIN(NBBLK,NSB-IOFFB)
      NAABB = NAA*NBB
      IF (TIMING) NOPOTH = NOPOTH+NAABB*NPAIRR
      CALL FZERO (Q(IDA),NAABB*NPAIRR)
C...  ALPHA REPLACEMENTS
      IF (TIMING) THEN
       T1 = SECOND()
       ICNT = 0
      END IF
      IATABB = IWA+IOFFA*MAXR3A
      IDAA = IDA-1-NAABB
      DO 61 IA=1,NAA
      IATAB = IATABB
      IF (TIMING) ICNT = ICNT+IQ(NWA+IOFFA-1+IA)
      DO 60 IWW=1,IQ(NWA+IOFFA-1+IA)
      ICCC = IOFFB + IQ(IATAB)
      IAAA = IDAA+IQ(IATAB+MAXR2A)*NAABB
      IF (IQ(IATAB+MAXREA).LT.0) THEN
C...  EQUATION (18)
      DO 40 IB=1,NBB
40    Q(IAAA+IB) = - C(ICCC+IB)
      ELSE
C...  EQUATION (18)
      DO 50 IB=1,NBB
50    Q(IAAA+IB) =   C(ICCC+IB)
      END IF
60    IATAB = IATAB+1
      IDAA = IDAA + NBB
61    IATABB = IATABB + MAXR3A
      IF (TIMING) THEN
       TCCA = TCCA + SECOND()-T1
       NOPCCA = NOPCCA + ICNT*NBB
       NLOCCA = NLOCCA + ICNT
      END IF
C...  BETA REPLACEMENTS
      IF (TIMING) THEN
       ICNT = 0
       T1 = SECOND()
      END IF
      IBTABB = IWB+IOFFB*MAXR3B
      ICCCC = IOFFCI(ISYMBE)+(IOFFA-1)*NSBE
      IDBB = IDB-NBB-NAABB
      DO 91 IB=1,NBB
      IBTAB = IBTABB
      IF (TIMING) ICNT = ICNT+IQ(NWB+IOFFB-1+IB)
      DO 90 IWW=1,IQ(NWB-1+IOFFB+IB)
      ICCC = ICCCC+IQ(IBTAB)
      IBBB = IDBB+IQ(IBTAB+MAXR2B)*NAABB
      IF (IQ(IBTAB+MAXREB).LT.0) THEN
C...  EQUATION (19)
      DO 70 IA=1,NAA
70    Q(IBBB+IA*NBB) = Q(IBBB+IA*NBB) - C(ICCC+IA*NSBE)
      ELSE
C...  EQUATION (19)
      DO 80 IA=1,NAA
80    Q(IBBB+IA*NBB) = Q(IBBB+IA*NBB) + C(ICCC+IA*NSBE)
      END IF
90    IBTAB = IBTAB + 1
      IDBB = IDBB + 1
91    IBTABB = IBTABB + MAXR3B
      IF (TIMING) THEN
       TCCB = TCCB + SECOND()-T1
       NOPCCB = NOPCCB + ICNT*NAA
       NLOCCB = NLOCCB + ICNT
      END IF
C
C...  MATRIX MULTIPLICATION WITH INTEGRALS
C...  EQUATION (20)
      IF (TIMING) T1 = SECOND()
      IF (IPRHAM.GT.3)
     1CALL OUTSQR (Q(IDA),NAABB,NAABB,NPAIRR,'DA')
      IF (IPRHAM.GT.3)
     1CALL OUTSQR (Q(INTAA+INTOFF(ISYMR)),NPAIRR,NPAIRR,NPAIRR,'INTAA')
      CALL FCMXM (Q(IDA),Q(INTAA+INTOFF(ISYMR)),Q(IEA),NAABB,NPAIRR,
     >NOPMXM)
      NOPMX = NOPMX+2*NAABB*NPAIRR**2
      IF (IPRHAM.GT.3)
     1CALL OUTSQR (Q(IEA),NAABB,NAABB,NPAIRR,'EA')
      IF (TIMING) TMXM = TMXM+SECOND()-T1
C
C...  ALPHA REPLACEMENTS
      IF (TIMING) THEN
       T1 = SECOND()
       ICNT = 0
      END IF
      IATABB = IWA+IOFFA*MAXR3A
      IEAA = IEA-1-NAABB
      DO 121 IA=1,NAA
      IATAB = IATABB
      IF (TIMING) ICNT = ICNT+IQ(NWA+IOFFA-1+IA)
      DO 120 IWW=1,IQ(NWA+IOFFA-1+IA)
      ICCC = IOFFB + IQ(IATAB)
      IAAA = IEAA+IQ(IATAB+MAXR2A)*NAABB
      IF (IQ(IATAB+MAXREA).LT.0) THEN
C...  EQUATION (21)
      DO 100 IB=1,NBB
100   S(ICCC+IB) = S(ICCC+IB) - Q(IAAA+IB)
      ELSE
C...  EQUATION (21)
      DO 110 IB=1,NBB
110   S(ICCC+IB) = S(ICCC+IB) + Q(IAAA+IB)
      END IF
120   IATAB = IATAB + 1
      IEAA = IEAA + NBB
121   IATABB = IATABB + MAXR3A
      IF (TIMING) THEN
       TCCA = TCCA + SECOND()-T1
       NOPCCA = NOPCCA + ICNT*NBB
       NLOCCA = NLOCCA + ICNT
      END IF
C...  BETA REPLACEMENTS
      IF (TIMING) THEN
       T1 = SECOND()
       ICNT = 0
      END IF
      IBTABB = IWB+IOFFB*MAXR3B
      ICCCC = IOFFCI(ISYMBE)+(IOFFA-1)*NSBE
      IEBB = IEB-NBB-NAABB
      DO 151 IB=1,NBB
      IBTAB = IBTABB
      IF (TIMING) ICNT = ICNT+IQ(NWB+IOFFB-1+IB)
      DO 150 IWW=1,IQ(NWB-1+IOFFB+IB)
      ICCC = ICCCC+IQ(IBTAB)
      IBBB = IEBB+IQ(IBTAB+MAXR2B)*NAABB
      IF (IQ(IBTAB+MAXREB).LT.0) THEN
C...  EQUATION (22)
      DO 130 IA=1,NAA
130   S(ICCC+IA*NSBE) = S(ICCC+IA*NSBE) - Q(IBBB+IA*NBB)
      ELSE
C...  EQUATION (22)
      DO 140 IA=1,NAA
140   S(ICCC+IA*NSBE) = S(ICCC+IA*NSBE) + Q(IBBB+IA*NBB)
      END IF
150   IBTAB = IBTAB + 1
      IEBB = IEBB + 1
151   IBTABB = IBTABB + MAXR3B
      IF (TIMING) THEN
       TCCB = TCCB + SECOND()-T1
       NOPCCB = NOPCCB + ICNT*NAA
       NLOCCB = NLOCCB + ICNT
      END IF
160   CONTINUE
      CALL CORLSR(IBASES)
170   CONTINUE
      CALL CORLSR(IBASE)
      IF (TIMING) THEN
      TTOT = SECOND()-T0
      TOTH = TTOT-TMXM-TCCA-TCCB-TONE
      NOPTOT = NOPOTH+NOPMXM+NOPCCA+NOPCCB+NOPONE
      CALL TPRINT ('String c. coeff',NOPONE,TONE,' ',0)
      CALL TPRINT ('Coupl coeffs(A)',NOPCCA,TCCA,'Av. v length'
     >,NOPCCA/NLOCCA)
      CALL TPRINT ('Coupl coeffs(B)',NOPCCB,TCCB,'Av. v length'
     >,NOPCCB/NLOCCB)
      IAVSP = NINT(1D2-(1D2*NOPMXM)/NOPMX)
      CALL TPRINT ('Matrix multiply',NOPMXM,TMXM,'Av. % sparsity',IAVSP)
      CALL TPRINT ('Other          ',NOPOTH,TOTH,' ',0)
      CALL TPRINT ('Total          ',NOPTOT,TTOT,' ',0)
      END IF
      IF (IPRHAM.GT.0)
     1CALL OUTVEC (S,NNCI,'S RETURNED BY SIGMA')
      RETURN
      END
      SUBROUTINE ONELAL (N,ITYPE,INTER,INTAD,ISYM,ISYMR,IW,NW,ITZ,MAXRE)
CSUBR CONSTRUCT AND STORE ALL ONE ELECTRON REPLACEMENTS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C...  MAXIMUM NUMBERS OF ORBITALS IN ANY SYMMETRY MUST BE LESS THANTHIS!
      PARAMETER (MXORB=64)
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      DIMENSION ITYPE(*),INTER(N,*),INTAD(*)
     1  ,IW(*),NW(*),ITZ(*)
C     1  ,IW(MAXRE,3,*),NW(*),ITZ(*)
      DIMENSION ICG(32),ICGI(33),JSTOR1(MXORB),JSTOR2(MXORB)
     >,JSTOR3(MXORB)
      MAXRE2 = MAXRE*2
      MAXRE3 = MAXRE*3
C...  SET COUNTERS TO ZERO AND SET UP LOOK-UP FOR BRAS
      NINSYM = IBINOM(NACT,N)
      ITZI = ICORI(NINSYM)
      NI = 0
      ICGI(1) = 0
      DO 23 II=1,2999999
      IISYM = 0
      CALL FCSTRS (N,ICGI,IISYM,ITYPE)
      IF (ICGI(1).EQ.0) GOTO 24
      IF (IISYM.NE.ISYM) GOTO 23
      NI=NI+1
      IQ(ITZI-1+II) = NI
23    CONTINUE
24    CONTINUE
      DO 25 I=1,NI
25    NW(I) = 0
      IF (NI.EQ.0) GOTO 999
C...  LOOP OVER INTERMEDIATE N-1 ELECTRON STRINGS K
      ICG(1) = 0
      DO 50 ISK=1,2999999
      ISYMK = 0
      CALL FCSTRS (N-1,ICG,ISYMK,ITYPE)
      IF (ICG(1).EQ.0) GOTO 999
      ICG(N) = NACT+1
      ISYMIO = MULT(ISYMK,ISYM)
      ISYMJO = MULT(ISYMR,ISYMIO)
      ISYMJ = MULT(ISYM,ISYMR)
C...  CONSTRUCT CREATIONS TO KETS J (SYMMETRY ISYMJ)
C      DO 777 IREP=1,10
      IREPL = 1
      IADR = 1
      IPARIT = 1
      DO 10 J=1,N-1
10    IADR = IADR + INTER(J+1,ICG(J))
      NJ = 0
      DO 30 JO=1,NACT
      IF(JO.EQ.ICG(IREPL))THEN
      IADR = IADR+INTER(IREPL,JO)-INTER(IREPL+1,JO)
      IREPL = IREPL+1
      IPARIT = -IPARIT
      ELSE IF (ITYPE(JO).EQ.ISYMJO) THEN
C...  GET CONNECTED N ELECTRON DETERMINANT J
C...  STORE INFO
      NJ = NJ+1
      JSTOR1(NJ) = ITZ(IADR+INTER(IREPL,JO))
      JSTOR2(NJ) = IPARIT
      JSTOR3(NJ) = JO
      END IF
30    CONTINUE
C777   CONTINUE
      IF (NJ.EQ.0) GOTO 50
C...  CONSTRUCT CREATIONS TO BRAS I (SYMMETRY ISYM)
      IREPL = 1
      IADR = ITZI
      IPARIT = 1
      DO 110 J=1,N-1
110   IADR = IADR + INTER(J+1,ICG(J))
      DO 130 IO=1,NACT
      IF(IO.EQ.ICG(IREPL))THEN
      IADR = IADR+INTER(IREPL,IO)-INTER(IREPL+1,IO)
      IREPL = IREPL+1
      IPARIT = -IPARIT
      ELSE IF (ITYPE(IO).EQ.ISYMIO) THEN
C...  GET CONNECTED N ELECTRON DETERMINANT I
      I = IQ(IADR+INTER(IREPL,IO))
      IOFFS = (I-1)*MAXRE3+NW(I)
      IOO = (IO-1)*NACT
C...  STORE INFO FOR EACH J
      IF (IPARIT.GT.0) THEN
CDIR$ SHORTLOOP
CDIR$ IVDEP
C$DIR NO_RECURRENCE
*VOPTION INDEP VEC
      DO 45 J=1,NJ
      IW(IOFFS       +J) = JSTOR1(J)
      IW(IOFFS+MAXRE2+J) = JSTOR3(J)+IOO
C      IW(IOFFS+MAXRE2+J) = INTAD(JSTOR3(J)+IOO)
45    IW(IOFFS+MAXRE +J) = JSTOR2(J)
      ELSE
CDIR$ SHORTLOOP
CDIR$ IVDEP
C$DIR NO_RECURRENCE
*VOPTION INDEP VEC
      DO 46 J=1,NJ
      IW(IOFFS       +J) = JSTOR1(J)
      IW(IOFFS+MAXRE2+J) = JSTOR3(J)+IOO
C      IW(IOFFS+MAXRE2+J) = INTAD(JSTOR3(J)+IOO)
46    IW(IOFFS+MAXRE +J) = -JSTOR2(J)
      END IF
      NW(I) = NW(I) + NJ
      END IF
130   CONTINUE
50    CONTINUE
999   CALL CORLSI(ITZI)
C...  GATHER LOOP FOR PACKED INTEGRAL ADDRESSES
      IOFFS = MAXRE2
      DO 147 I=1,NI
      DO 148 J=1,NW(I)
148   IW(IOFFS+J) = INTAD(IW(IOFFS+J))
147   IOFFS = IOFFS + MAXRE3
      RETURN
      END
      SUBROUTINE FCDIAG (S)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /TAPES / INP,IOUT
      COMMON /CHAMIL/ ISYMH,NNCI
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      DIMENSION S(*)
      COMMON /CLOCAL/ ICGA(32),ICGB(32)
      IF (IPRTIM.GE.0) T0 = SECOND()
      NOP = 0
      IF (IPRDIA.GT.0) THEN
C      WRITE (IOUT,*) 'COMPUTE DIAGONAL ELEMENTS'
      CALL OUTVEC (Q(INTHAA),NACT,'HAA')
      CALL OUTSQR (Q(INTJAA),NACT,NACT,NACT,'INTJAA')
      CALL OUTSQR (Q(INTKAA),NACT,NACT,NACT,'INTKAA')
      END IF
C      IF (IPRDIA.GT.1) WRITE (IOUT,*) 'MAXAA,MAXBB ',MAXAA,MAXBB
      IZB = ICORR(NACT*MAXBB)
      IZA = ICORR(NACT*MAXAA)
      IZ  = ICORR(NACT*MAXBB)
      IF  = ICORR(MAXBB)
      ICCC=0
      DO 220 ISYMB=1,8
      ISYMA=MULT(ISYMB,ISYMH)
      NSA = NSTRA(ISYMA)
      NSB = NSTRB(ISYMB)
      IF (IPRDIA.GT.1) WRITE (IOUT,*) 'ISYMB,ISYMA,NSA,NSB ',ISYMB
     1,ISYMA,NSA,NSB
      IF (NSA*NSB.LE.0) GOTO 220
      CALL FZERO (Q(IZB),NACT*NSB)
C..   STORE BETA STRING OCCUPANCIES
      ICGB(1) = 0
      DO 20 IB=1,NSB
      CALL FCSTRS (NB,ICGB,ISYMB,ITYPEA)
      DO 10 II=1,NB
10    Q(IZB-1+IB+(ICGB(II)-1)*NSB)=1
20    CONTINUE
      IF (IPRDIA.GT.2) CALL OUTSQR (Q(IZB),NSB,NSB,NACT,'ZB')
C
C..   STORE ALPHA STRING OCCUPANCIES
      ICGA(1) = 0
      CALL FZERO (Q(IZA),NACT*NSA)
      DO 40 IA=1,NSA
      CALL FCSTRS (NA,ICGA,ISYMA,ITYPEA)
      DO 40 II=1,NA
40    Q(IZA-1+IA+(ICGA(II)-1)*NSA) = 1
      IF (IPRDIA.GT.2) CALL OUTSQR (Q(IZA),NSA,NSA,NACT,'ZA')
C
      DO 50 IAB=1,NSA*NSB
50    S(ICCC+IAB) = 0
C
C...   SPECIAL CODE FOR RHF ORBITALS
      DO 210 IA=1,NSA
C...  OCCUPATION NUMBERS OF DETERMINANTS
      DO 100 I=1,NACT
      DO 100 IB=1,NSB
100   Q(IZ-1+IB+(I-1)*NSB) = Q(IZB-1+IB+(I-1)*NSB)
     1   + Q(IZA-1+IA+(I-1)*NSA)
      IF (IPRDIA.GT.3) CALL OUTSQR (Q(IZ),NSB,NSB,NACT,'Z')
C
C...  CONSTANT EXCHANGE CONTRIBUTION & COULOMB
C...  1/2  SUM(IJ) N(I) N(J) ((II|JJ)-0.5(IJ|JI))
C...  EQUATIONS (23), (25)
      IJ = -1
      DO 130 I=1,NACT
      ZZ = Q(INTHAA-1+I)
      IF (IPRDIA.GT.3) WRITE (IOUT,*) '1567;I,HAA(I) ',I,ZZ
      DO 110 IB=1,NSB
110   S(ICCC+IB) = S(ICCC+IB) + ZZ*Q(IZ-1+IB+(I-1)*NSB)
      DO 130 J=1,I
      IJ = (I-1)*NACT-1+J
      ZZ = Q(INTJAA+IJ)-0.5D0*Q(INTKAA+IJ)
      IF (I.EQ.J) ZZ = ZZ*0.5D0
      IF (IPRDIA.GT.3) WRITE (IOUT,*) '1567; I,J,ZZ ',I,J,ZZ
      DO 120 IB=1,NSB
120    S(ICCC+IB) = S(ICCC+IB)
     1 + ZZ*Q(IZ-1+IB+(I-1)*NSB)*Q(IZ-1+IB+(J-1)*NSB)
130   CONTINUE
C
      IF (IPRDIA.GT.3) CALL OUTVEC (S(ICCC+1),NSB,'S AFTER CONSTS')
C...  VV = 4*S**2, SEE EQUATION (26)
      VV = (NA-NB)**2
C..   Z TO HOLD 1 IF SINGLY OCC, 0 OTHERWISE
C...  (CAPITAL N = N(2-N), EQUATION (25))
      CALL FZERO (Q(IF),NSB)
      DO 150 I=1,NACT
      DO 140 IB=1,NSB
140   Q(IZ-1+IB+(I-1)*NSB) = Q(IZ-1+IB+(I-1)*NSB)
     1 * (2.0D0-Q(IZ-1+IB+(I-1)*NSB))
      DO 150 IB=1,NSB
150   Q(IF-1+IB) = Q(IF-1+IB) + Q(IZ-1+IB+(I-1)*NSB)
      DO 160 IB=1,NSB
160   Q(IF-1+IB) = MAX(Q(IF-1+IB),1.1D0)
      DO 170 IB=1,NSB
170   Q(IF-1+IB) = (VV-Q(IF-1+IB))/( Q(IF-1+IB) * (Q(IF-1+IB)-1.0D0) )
C..   F NOW HOLDS F FACTOR, EQUATION (26)
C
      IF (IPRDIA.GT.3) CALL OUTVEC (Q(IF),NSB,'F')
      IJ = -1
      DO 200 I=1,NACT
      IF (I.EQ.1) GOTO 190
      DO 180 J=1,I-1
      IJ = (I-1)*NACT-1+J
      ZZ = -0.5D0*Q(INTKAA+IJ)
      DO 180 IB=1,NSB
180   S(ICCC+IB) = S(ICCC+IB)
     1 + ZZ*Q(IF-1+IB)*Q(IZ-1+IB+(I-1)*NSB)*Q(IZ-1+IB+(J-1)*NSB)
190   IJ = (I-1)*NACT-1+I
      ZZ = -0.25D0*Q(INTKAA+IJ)
      IF (IPRDIA.GT.3) WRITE (IOUT,*) '1875; I,ZZ ',I,ZZ
      DO 200 IB=1,NSB
200   S(ICCC+IB) = S(ICCC+IB) + ZZ*Q(IZ-1+IB+(I-1)*NSB)
210   ICCC = ICCC+NSB
      NOP = NOP+NSA*NSB*(NACT*7+3*NACT*(NACT+1)+5)
220   CONTINUE
      CALL CORLSR(IZB)
      IF (IPRDIA.GE.0) CALL OUTVEC (S,ICCC,'DIAGONAL ELEMENTS')
      IF (IPRTIM.GE.0) CALL TPRINT ('Hamiltonian diag. els.'
     >,NOP,SECOND()-T0,' ',0)
      RETURN
      END
      FUNCTION INICOR(MEMREQ)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER CORMON,CORLSR,CORLSI,ENDCOR
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      COMMON /TAPES / INP,IOUT
      COMMON /CORCTL/ INTREL,LWORD,LTOP,LMAX,LMIN,MREAL,MXBFF,IBASE,MEM
      SAVE ICOUNT
      DATA ICOUNT/0/
C
C...  EXECUTE ENDCOR CODE IF NOT FIRST CALL TO THIS ROUTINE
      IF (ICOUNT.NE.0) CALL FMAIN (Q(IBASE+1),MEM)
      ICOUNT = ICOUNT+1
C
      MEM = MEMREQ
      LENN = MEM+2
C
C
      CALL GMAINV (Q,IBASE,LENN)
      IBASE = IBASE-1
      write(iout,*)'lenn,mem',lenn,mem
      IF (LENN.LT.MEM) THEN
       WRITE (IOUT,'('' Request for'',I9,'' words refused'',
     >  '', '',I9,'' available'')') MEM,LENN
       STOP 'Insufficient memory'
      END IF
      WRITE (IOUT,'('' Variable memory set to '',I10,'' words'')') MEM
      LWORD = MEM+IBASE
      LTOP = IBASE
      LMAX = IBASE
      LMIN = IBASE
      INICOR = LWORD-IBASE
      RETURN
C
      ENTRY ICORIM ()
      ICORIM = (LWORD - LTOP) * INTREL
      RETURN
C
      ENTRY ICORRM ()
      ICORRM = LWORD - LTOP
      RETURN
C
      ENTRY ICORI (NWOR)
      NWORD = IABS(NWOR)
      NW = (NWORD+INTREL-1)/INTREL
      MFR = LWORD - LTOP
      IF (NW .GT. MFR) GOTO 20
      ICORI = LTOP*INTREL+1
      LTOP=LTOP+NW
      IF (NWOR.LT.0) LMAX = MAX0(LMAX,LTOP)
      IF (NWOR.GT.0) LMIN = MAX0(LMIN,LTOP)
      RETURN
C
      ENTRY ICORR (NWOR)
      NWORD = IABS(NWOR)
      NW = NWORD
      MFR = LWORD - LTOP
      IF (NW .GT. MFR) GOTO 20
      ICORR = LTOP+1
      LTOP=LTOP+NW
      IF (NWOR.LT.0) LMAX = MAX0(LTOP,LMAX)
      IF (NWOR.GT.0) LMIN = MAX0(LTOP,LMIN)
      RETURN
C
20    CONTINUE
      WRITE (IOUT,*) 'insufficient memory available - require ',NW,
     >               ' have ',MFR
      IF (NW .NE. NWORD) THEN
         WRITE(IOUT,*) 'the request was for ',NWORD,' integer words'
      ELSE
         WRITE(IOUT,*) 'the request was for real words'
      END IF
       STOP 'Insufficient memory'
      END
      SUBROUTINE CORLSI (IAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      COMMON /TAPES / INP,IOUT
      COMMON /CORCTL/ INTREL,LWORD,LTOP,LMAX,LMIN,MREAL,MXBFF,IBASE,MEM
      LTOP = (IAD-2+INTREL)/INTREL
      RETURN
C
      ENTRY CORLSR (IAD)
      LTOP = IAD - 1
      RETURN
C
      ENTRY CORMON()
      WRITE (IOUT,40) LMIN-IBASE,LWORD-LMIN
40    FORMAT(/1X,'Minimal core high water =',I8,
     1' real numbers (spare core =',I8,')')
      RETURN
C
      ENTRY ENDCOR()
      CALL FMAIN (Q(IBASE+1),MEM)
      WRITE (IOUT,'('' Variable memory released'')')
      RETURN
      END
      FUNCTION FCISTR (IPAR,INTER,ICG,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ICG(N),INTER(N,*)
      DIMENSION ICGG(32)
      DO 10 J=1,N
10    ICGG(J) = ICG(J)
      IF (N.LT.2) GOTO 40
      DO 30 JJ=2,N
      JJM = JJ-1
      I = JJM
      DO 20 II=1,JJM
      IF (ICGG(I).LT.ICGG(I+1)) GOTO 20
      J1=ICGG(I)
      ICGG(I)=ICGG(I+1)
      ICGG(I+1)=J1
      IPAR=-IPAR
20    I=I-1
30    CONTINUE
40    CONTINUE
      FCISTR=1
      DO 50 J=1,N
50    FCISTR = FCISTR+INTER(J,ICGG(J))
      RETURN
      END
      SUBROUTINE FCMIC (IC,ITZ,ISYM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      DIMENSION ICG(32)
      DIMENSION IC(*),ITZ(*)
      DO 1 I=1,8
      NSTRA(I) = 0
1     NSTRB(I) = 0
      NSTRAT = IBINOM(NACT,NA)
      NSTRBT = IBINOM(NACT,NB)
C...  ALPHA STRING INFORMATION
      ICG(1) = 0
      DO 11 MT=1,NSTRAT
      JMT = 0
      CALL FCSTRS (NA,ICG,JMT,ITYPEA)
      NSTRA(JMT) = NSTRA(JMT)+1
C...  STORE SYMMETRY FOR LATER PROCESSING
11    IC(MT) = -JMT
C...  BETA STRING INFORMATION
      ICG(1) = 0
      DO 10 MT=1,NSTRBT
      JMT = 0
      CALL FCSTRS (NB,ICG,JMT,ITYPEA)
      NSTRB(JMT) = NSTRB(JMT)+1
10    ITZ(MT) = NSTRB(JMT)
C
C...  LOOP OVER SYMMETRIES BUILDING TOTAL OFFSETS
      NNCI = 0
      DO 110 ISYMB=1,8
      ISYMA = MULT(ISYM,ISYMB)
      DO 120 MT=1,NSTRAT
      IF (IC(MT).NE.-ISYMA) GOTO 120
      IC(MT) = NNCI
      NNCI = NNCI + NSTRB(ISYMB)
120   CONTINUE
110   CONTINUE
      RETURN
      END
      SUBROUTINE FCSPAD (C,IREF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      COMMON /BIG/ Q(1)
      INTEGER IQ(1)
      EQUIVALENCE (Q(1),IQ(1))
      DIMENSION C(*)
      COMMON /CHAMIL/ ISYMH,NNCI
      DIMENSION ICGA(32),ICGB(32),IDOCA(32)
     1         ,IDOCB(32),IALP(32),IBET(32)
      DO 10 I=1,NNCI
10    C(I)=0
      IBASE = ICORR(0)
      INTERA = ICORI(NA*NACT)
      INTERB = ICORI(NB*NACT)
      IC     = ICORI(NSTRAA)
      ITZ    = ICORI(NSTRBB)
      CALL INITI (IQ(INTERA),NA,NACT)
      CALL INITI (IQ(INTERB),NB,NACT)
      CALL FCMIC (IQ(IC),IQ(ITZ),ISYMH)
      CALL FCSTRG (IREF,ISYMH,ICGA,ICGB)
      NDOC=0
      NALP=0
      NBET=0
      DO 70 I=1,NACT
      ICASE=1
      DO 20 J=1,NA
      IF (ICGA(J).NE.I) GOTO 20
      ICASE=3
      JA=J
20    CONTINUE
      DO 30 J=1,NB
      IF (ICGB(J).NE.I) GOTO 30
      ICASE=ICASE+1
      JB=J
30    CONTINUE
      GOTO (70,40,50,60),ICASE
40    NBET=NBET+1
      IBET(NBET)=JB
      GOTO 70
50    NALP=NALP+1
      IALP(NALP)=JA
      GOTO 70
60    NDOC=NDOC+1
      IDOCA(NDOC)=JA
      IDOCB(NDOC)=JB
70    CONTINUE
      FAC = 1.0D0/SQRT(DBLE(2**NBET))
      IF (NBET.GT.5) STOP 5
      DO 120 I5=1,2
      DO 110 I4=1,2
      DO 100 I3=1,2
      DO 90 I2=1,2
      DO 80 I1=1,2
      IPAR = 1
      MTA = FCISTR (IPAR,IQ(INTERA),ICGA,NA)
      MTB = FCISTR (IPAR,IQ(INTERB),ICGB,NB)
      C(IQ(IC-1+MTA)+IQ(ITZ-1+MTB)) = DBLE(IPAR)*FAC
      IF (NBET.LT.1) GOTO 130
      J1=ICGB(IBET(1))
      ICGB(IBET(1))=ICGA(IALP(1))
      ICGA(IALP(1))=J1
80    CONTINUE
      IF (NBET.LT.2) GOTO 130
      J1=ICGB(IBET(2))
      ICGB(IBET(2))=ICGA(IALP(2))
      ICGA(IALP(2))=J1
90    CONTINUE
      IF(NBET.LT.3) GOTO 130
      J1=ICGB(IBET(3))
      ICGB(IBET(3))=ICGA(IALP(3))
      ICGA(IALP(3))=J1
100   CONTINUE
      IF (NBET.LT.4) GOTO 130
      J1=ICGB(IBET(4))
      ICGB(IBET(4))=ICGA(IALP(4))
      ICGA(IALP(4))=J1
110   CONTINUE
      IF (NBET.LT.5) GOTO 130
      J1=ICGB(IBET(5))
      ICGB(IBET(5))=ICGA(IALP(5))
      ICGA(IALP(5))=J1
120   CONTINUE
130   CONTINUE
      CALL CORLSR (IBASE)
      RETURN
      END
      SUBROUTINE FCSTRG (II,ISYM,ICGA,ICGB)
CSUBR EXTRACT STRINGS ICGA,ICGB FOR CONFIGURATION II IN SYMMETRY ISYM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      DIMENSION ICGA(*),ICGB(*)
      III = II
      DO 10 ISYMB=1,8
      ISYMA = MULT(ISYMB,ISYM)
      NSB = NSTRB(ISYMB)
      NSA = NSTRA(ISYMA)
      IF (III.LE.NSB*NSA) GOTO 20
10    III = III - NSB*NSA
20    ISA = (III-1)/NSB+1
      ISB = III - (ISA-1)*NSB
      ICGA(1) = 0
      DO 30 JSA=1,ISA
30    CALL FCSTRS (NA,ICGA,ISYMA,ITYPEA)
      ICGB(1) = 0
      DO 40 JSB=1,ISB
40    CALL FCSTRS (NB,ICGB,ISYMB,ITYPEA)
      RETURN
      END
      SUBROUTINE FCSTRS (N,ICG,ISYM,ITYPE)
C... GENERATE A NEW ALPHA OR BETA STRING WITH SYMMETRY ISYM
C...  ICG SHOULD BE PRESERVED BETWEEN CALLS
C... IF (ISYM.EQ.0) ON ENTRY, THEN JUST NEXT STRING GENERATED
C...  ITS SYMMETRY RETURNED IN ISYM
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      DIMENSION ICG(*),ITYPE(*)
C
      IF (ICG(1).LE.0) THEN
       DO 10 I=1,N
10     ICG(I)=I
       IF (N.LE.0) THEN
        DO 111 JJ=1,NACT
        IF (ITYPE(JJ).EQ.1) GOTO 112
111     CONTINUE
112     ICG(1) = JJ
       END IF
       GOTO 40
      END IF
C
21    DO 20 IU=N,1,-1
      ICG(IU) = ICG(IU) + 1
      IF (ICG(IU).LE.NACT+IU-N) GOTO 30
20    CONTINUE
      GOTO 60
CDIR$ NEXTSCALAR
*VOPTION NOVEC
30    DO 31 JU=IU+1,N
31    ICG(JU) = ICG(IU)+JU-IU
C
40    JMT = ITYPE(ICG(1))
      DO 50 J=2,N
50    JMT = MULT(JMT,ITYPE(ICG(J)))
      IF (JMT.NE.ISYM.AND.ISYM.NE.0) GOTO 21
      IF (ISYM.EQ.0) ISYM = JMT
      RETURN
60    ICG(1) = 0
      RETURN
      END
      FUNCTION IBINOM (M,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      TOP = M
      NN = MIN0(M-N,N)
      IF (NN.LT.0) GOTO 20
      BOT = NN
      BINOM = 1.0D0
      DO 10 I=1,NN
      BINOM = BINOM * TOP/BOT
      TOP = TOP - 1.0D0
10    BOT = BOT - 1.0D0
      IBINOM = BINOM + 0.1D0
      IF (DABS(BINOM-DBLE(IBINOM)).GT.0.5D0)
     1 STOP 'ARITHMETIC PROBLEM IN COMPUTING BINOMIAL COEFFICIENT'
      RETURN
20    IBINOM = 0
      RETURN
      END
      SUBROUTINE INITI (INTER,NITEM,NBOX)
CSUBR SETS UP PARTIAL WEIGHT ARRAY FOR ADDRESSING BINOMIAL DISTRIBUTIONS
      DIMENSION INTER(NITEM,*)
      N1=NITEM+1
      DO 30 K=1,NITEM
      DO 10 L=1,NBOX
10    INTER(K,L)=0
      DO 20 L=K,NBOX-NITEM+K-1
      INTER(K,L+1) = IBINOM(NBOX-L,NITEM-K)+INTER(K,L)
20    CONTINUE
30    CONTINUE
      DO 40 K=1,NITEM-1
      DO 40 L=K,NBOX-NITEM+K
40    INTER(K,L) = INTER(K,L) - INTER(K+1,L+1)
      DO 50 L=NITEM,NBOX
50    INTER(NITEM,L) = L-NITEM
      RETURN
      END
      SUBROUTINE FZERO (A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(*)
      DO 1 I=1,N
1     A(I) = 0D0
      RETURN
      END
      SUBROUTINE OUTSQR (Q,IDIM,IA,IB,TITLE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*(*) TITLE
      DIMENSION Q(IDIM,1)
      CHARACTER*10 FORMAT
      DATA FORMAT /'(/1X,A000)'/
      WRITE (FORMAT(7:9),'(I3)') 55+LEN(TITLE)/2
      WRITE (6,FORMAT) TITLE
      M=1
      N=8
10    IF (IB.LT.M) RETURN
      N=MIN0(N,IB)
      WRITE (6,30) (I,I=M,N)
      WRITE (6,40)
      DO 20 J=1,IA
      WRITE (6,50) J,(Q    (J,I) ,I=M,N)
20    CONTINUE
      M=M+8
      N=N+8
      GOTO 10
30    FORMAT(//3X,8I14)
40    FORMAT(/)
50    FORMAT(7X,I3,8F14.7)
      END
      SUBROUTINE OUTVEC (P,IB,TITLE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /TAPES / INP,IOUT
      CHARACTER*(*) TITLE
      DIMENSION P(1)
      WRITE (IOUT,'(/1X,(A))') TITLE
      M=1
      N=10
10    IF (IB.LT.M) RETURN
      N=MIN0(N,IB)
      WRITE(IOUT,20) M,N,(P    (I) ,I=M,N)
      M=M+10
      N=N+10
      GOTO 10
20    FORMAT(I4,'-',I3,10F12.6)
      END
      SUBROUTINE DIAG2(M,N,D,X)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXDIM=1000)
      DIMENSION   D(M), X(M,M)
      DIMENSION   E(MAXDIM)
      IF(M.GT.MAXDIM) THEN
        WRITE(6,10) M,MAXDIM
10      FORMAT(' DIMENSION TOO LARGE IN DIAG2:',2I8)
        STOP
      END IF
CSTART VAX UNIX-HP UNIX-CONVEX CRAY IBM UNIVAC
      EPS=2.5E-13
      TOL=2.5D-18
CEND
C
      IF(N.EQ.1) GO TO 400
C
C     HOUSEHOLDER'S REDUCTION
C     SIMULATION OF LOOP DO 150 I=N,2,(-1)
C
      DO 150 NI=2,N
      I=N+2-NI
      L=I-2
      H=0.0
      G=X(I,I-1)
      IF(L) 140,140,20
   20 DO 30 K=1,L
   30 H=H+X(I,K)**2
      S=H+G*G
      IF(S.GE.TOL) GO TO 50
   40 H=0.0
      GO TO 140
   50 IF(H) 140,140,60
   60 L=L+1
      F=G
      G=DSQRT(S)
      IF(F) 75,75,70
   70 G=-G
   75 H=S-F*G
      X(I,I-1)=F-G
      F=0.0
C
      DO 110 J=1,L
      X(J,I)=X(I,J)/H
      S=0.0
      DO 80 K=1,J
   80 S=S+X(J,K)*X(I,K)
      J1=J+1
      IF(J1.GT.L) GO TO 100
      DO 90 K=J1,L
   90 S=S+X(K,J)*X(I,K)
  100 E(J)=S/H
  110 F=F+S*X(J,I)
C
      F=F/(H+H)
C
      DO 120 J=1,L
  120 E(J)=E(J)-F*X(I,J)
C
      DO 130 J=1,L
      F=X(I,J)
      S=E(J)
      DO 130 K=1,J
  130 X(J,K)=X(J,K)-F*E(K)-X(I,K)*S
C
  140 D(I)=H
  150 E(I-1)=G
C
C     ACCUMULATION OF TRANSFORMATION MATRICES
C
  160 D(1)=X(1,1)
      X(1,1)=1.0
      DO 220 I=2,N
      L=I-1
      IF(D(I)) 200,200,170
  170 DO 190 J=1,L
      S=0.0
      DO 180 K=1,L
  180 S=S+X(I,K)*X(K,J)
      DO 190 K=1,L
  190 X(K,J)=X(K,J)-S*X(K,I)
  200 D(I)=X(I,I)
      X(I,I)=1.0
  210 DO 220 J=1,L
      X(I,J)=0.0
  220 X(J,I)=0.0
C
C     DIAGONALIZATION OF THE TRIDIAGONAL MATRIX
C
      B=0.0
      F=0.0
      E(N)=0.0
C
      DO 340 L=1,N
      H=EPS*(DABS(D(L))+DABS(E(L)))
      IF (H.GT.B) B=H
C
C     TEST FOR SPLITTING
C
      DO 240 J=L,N
      IF (DABS(E(J)).LE.B) GOTO 250
  240 CONTINUE
C
C     TEST FOR CONVERGENCE
C
  250 IF(J.EQ.L) GO TO 340
C
C     SHIFT FROM UPPER 2*2 MINOR
C
  260 P=(D(L+1)-D(L))*0.5/E(L)
      R=DSQRT(P*P+1.0)
      IF(P) 270,280,280
  270 P=P-R
      GO TO 290
  280 P=P+R
  290 H=D(L)-E(L)/P
      DO 300 I=L,N
  300 D(I)=D(I)-H
      F=F+H
C
C     QR TRANSFORMATION
C
      P=D(J)
      C=1.0
      S=0.0
C
C     SIMULATION OF LOOP DO 330 I=J-1,L,(-1)
C
      J1=J-1
      DO 330 NI=L,J1
      I=L+J1-NI
      G=C*E(I)
      H=C*P
C
C     PROTECTION AGAINST UNDERFLOW OF EXPONENTS
C
      IF (DABS(P).LT.DABS(E(I))) GOTO 310
      C=E(I)/P
      R=DSQRT(C*C+1.0)
      E(I+1)=S*P*R
      S=C/R
      C=1.0/R
      GO TO 320
  310 C=P/E(I)
      R=DSQRT(C*C+1.0)
      E(I+1)=S*E(I)*R
      S=1.0/R
      C=C/R
  320 P=C*D(I)-S*G
      D(I+1)=H+S*(C*G+S*D(I))
      DO 330 K=1,N
      H=X(K,I+1)
      X(K,I+1)=X(K,I)*S+H*C
  330 X(K,I)=X(K,I)*C-H*S
C
      E(L)=S*P
      D(L)=C*P
      IF (DABS(E(L)).GT.B) GO TO 260
C
C     CONVERGENCE
C
  340 D(L)=D(L)+F
C
C     ORDERING OF EIGENVALUES
C
      NI=N-1
  350 DO 380I=1,NI
      K=I
      P=D(I)
      J1=I+1
      DO 360J=J1,N
      IF(D(J).GE.P) GOTO 360
      K=J
      P=D(J)
  360 CONTINUE
      IF (K.EQ.I) GOTO 380
      D(K) =D(I)
      D(I)=P
      DO 370 J=1,N
      P=X(J,I)
      X(J,I)=X(J,K)
  370 X(J,K)=P
  380 CONTINUE
C
C     FIXING OF SIGN
C
      DO 385 I=1,N
      PM=0
      DO 386 J=1,N
      IF(PM.GT.DABS(X(J,I))) GOTO 386
      PM =DABS(X(J,I))
      K=J
  386 CONTINUE
      IF(X(K,I).GE.0) GOTO 385
      DO 387 J=1,N
  387 X(J,I)=-X(J,I)
  385 CONTINUE
  390 GO TO 410
C
C     SPECIAL TREATMENT OF CASE N = 1
C
  400 D(1)=X(1,1)
      X(1,1)=1.0
  410 RETURN
      END
      SUBROUTINE TPRINT (TITLE,NOP,T,STITLE,NS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /TAPES / INP,IOUT
      COMMON /CFULLR/ CORE,THR,THRRES
      COMMON /CFULL / MULT(8,8),ICF(256),NPAIR(8),NSTRA(8),NSTRB(8)
     >               ,INTOFF(8),ITYPEA(256)
     >               ,NOCC(8),NCI(8),NELEC,MS2,NA,NB,NACT
     >               ,INTAA,INTHAA,INTJAA,INTKAA,NINT2
     >               ,NSTRAA,NSTRBB,MAXAA,MAXBB,MAXREP,MAXPAR,INTADT
     >               ,IPRINT,IPRDIA,IPRHAM,IPRTIM,MAXIT,NROOT,MEMORY
      CHARACTER*(*) TITLE,STITLE
      FLOP = 0
      IF (T.NE.0D0) FLOP = 1D-6 * NOP/T
      IF (STITLE.EQ.' ') THEN
      WRITE (IOUT,666) TITLE,NOP,T,FLOP
      ELSE
      WRITE (IOUT,666) TITLE,NOP,T,FLOP,STITLE,NS
      END IF
666   FORMAT(' Timing: ',A,5X,'Operations',I12,5X,'Time',F10.2,
     >5X,'Mflops',F10.3,:,5X,A,I6)
      RETURN
      END
CFile: standard.f
      SUBROUTINE GMAINV (QQ,IBASE,LENN)
C.....SUBSTITUTE FOR DYNAMIC MEMORY ALLOCATION
C.... THIS ROUTINE MUST BE CALLED WITH QQ FIRST WORD IN COMMON/BIG/
C....  THIS ROUTINE IS NOT FOOLPROOF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXM=8000000)
      COMMON/BIG/Q(MAXM)
      COMMON /CORCTL/ INTREL,ICORCT(8)
C...  INTREL MUST BE NUMBER OF INTEGERS PER REAL
      INTREL = 2
C
      LENN = MIN(MAXM,LENN)
      IBASE = 1
      RETURN
      ENTRY FMAIN (QQ,LENN)
C.....SUBSITUTE FOR MEMORY RELEASE
      RETURN
      END
      FUNCTION SECOND()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C.....SHOUld RETURN CPU TIME IN SECONDS
      SAVE T
      DATA T/0D0/
      SECOND=T
      T=T+.01D0
      RETURN
      END
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE UNPACK (IJ,I,J)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     INPUT    SYMMETRISIERTER INDEX  IJ
C     OUTPUT    I  UND  J
C
      I=(DSQRT(DFLOAT(2*IJ)+0.25D0)+0.4999D0)
      J=IJ-(I*(I-1))/2
      RETURN
c     END !UNPACK
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CFile: cray.f
C      SUBROUTINE GMAINV (Q,IBASE,LENN)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      COMMON /CORCTL/ INTREL,ICORCT(8)
C      DIMENSION Q(*)
C      INTREL = 1
C      CALL HPSHRINK
C      IF (LENN.LE.0) THEN
C        LENN = 0
C        IBASE = 1
C        RETURN
C      END IF
C      CALL HPALLOC(IBAES,LENN,IERR,1)
C      IBASE = IBAES-LOC(Q)+1
C      RETURN
C      ENTRY FMAIN (Q,LENN)
CC     SHOULD CALL WITH LENN SAME VALUE AS GMAINV
C      IF (LENN.LE.0) RETURN
C      IADR = LOC(Q)
C      CALL HPDEALLC(IADR,IERR,1)
C      CALL HPSHRINK
C      RETURN
C      END
CFile: ibm.asm
C         MACRO
C&NAME    ROUTINE
C         GBLA    &ROUTINC
C         GBLC    &ROUTINN
C         AIF   ('&ROUTINC' NE '0').NODEF
CR0       EQU   0
CR1       EQU   1
CR2       EQU   2
CR3       EQU   3
CR4       EQU   4
CR5       EQU   5
CR6       EQU   6
CR7       EQU   7
CR8       EQU   8
CR9       EQU   9
CR10      EQU   10
CR11      EQU   11
CR12      EQU   12
CR13      EQU   13
CR14      EQU   14
CR15      EQU   15
C.NODEF   ANOP
C&ROUTINC SETA &ROUTINC+1
C         LTORG
C&NAME    CSECT
C&ROUTINN SETC '&NAME'
C         B     88(0,R15)
C         DC    X'09',CL9'&NAME'
C         DS    18F
C         STM   R14,R12,12(R13)
C         LR    R12,R13
C         LA    R13,16(0,R15)
C         ST    R12,4(0,R13)
C         ST    R13,8(0,R12)
C         USING &NAME.+16,R13
C         MEND
C         MACRO
C&NAME    EXIT  &ARG
C         AIF   ('&NAME' EQ '').NONAME
C&NAME    EQU   *
C.NONAME  AIF   ('&ARG' EQ '').R15DONE
C         AIF   ('ARG'(1,1) EQ '(').R15REG
C         LA    R15,&ARG
C         AGO   .R15DONE
C.R15REG  ANOP
C.R15DONE ANOP
C         L     R13,4(R13)
C         L     R14,12(R13)
C         MVI   12(R13),255
C         LM    R2,R12,28(R13)
C         BR    R14
C         MEND
C         MACRO
C&NAME    ENTRYPT
C         GBLC  &ROUTINN
C         AIF   ('&ROUTINN' NE '&SYSECT').ERROR1
C         ENTRY &NAME
C         CNOP  0,8
C&NAME    B     14(0,R15)
C         DC    X'09',CL9'&NAME'
C         STM   R14,R12,12(R13)
C         LR    R12,R13
C         LA    R11,&NAME-&SYSECT-16
C         LR    R13,R15
C         SR    R13,R11
C         ST    R12,4(0,R13)
C         ST    R13,8(0,R12)
C         USING &SYSECT.+16,R13
C         MEXIT
C.ERROR1  MNOTE 12,'ENTRYPT must be used in CSECT started with ROUTINE'
C         MEND
C*
CSECOND   ROUTINE
C         CLI   SECONDFL,C'N'
C         BE    SECOND1
C         TTIMER
C         L     R1,SECONDT
C         SR    R1,R0
C         ST    R1,SECONDB
C         LD    0,SECONDF
C         DD    0,=D'384.0E2'
C         EXIT
CSECOND1  MVI   SECONDFL,C'Y'
C         STIMER TASK,TUINTVL=SECONDT
C         SDR   0,0
C         EXIT
CSECONDF  DS    0D
C         DC    X'4E000000'         unnormalised zero
CSECONDB  DS    F
CSECONDT  DC    X'7FFFFFFF'
CSECONDFL DC    C'N'
C         LTORG
C****************************
CCORCTL   COM
CINTREL   DS    F
C         DS    8F
CGMAINV   ROUTINE
C*
C* CALL GMAINV (Q,IAD,LEN)
C* do V-type GETMAIN; length in 8-byte words returned in LEN;
C* IAD returned such that first address is Q(IAD), Q r*8;
C* On entry, LEN should hold the maximum no. of words wanted.
C*
C         L     R12,=A(CORCTL)
C         USING CORCTL,R12
C         LA    2,2
C         ST    2,INTREL
C         LM    R2,R4,0(R1)
C         L     R1,0(R4)
C         SLL   R1,3
C         ST    R1,GMAINV3+4
C         GETMAIN VU,LA=GMAINV3,A=GMAINV4
C         L     R1,GMAINV4+4
C         SRL   R1,3
C         ST    R1,0(R4)              length in 8 byte words
C         L     R1,GMAINV4
C         SR    R1,2
C         SRL   R1,3
C         LA    R1,1(R1)
C         ST    R1,0(R3)              position in real*8 array Q
C         EXIT  0
CGMAINV3  DC    F'8'
C         DS    F
CGMAINV4  DC    F'0,0'
CFMAIN    ENTRYPT
C*
C* CALL FMAIN  (Q,LEN)
C* do FREEMAIN; 1st address at Q(1) (Q is r*8); length=LEN*8 bytes
C*
C         L     R2,0(R1)
C         L     R3,4(R1)
C         L     R3,0(R3)
C         LTR   R3,R3
C         BH    FMAINOK
C         BL    FMAINERR
C         EXIT  0
CFMAINOK  SLL   R3,3
C         FREEMAIN R,LV=(R3),A=(R2)
C         EXIT  0
CFMAINERR WTO   'FMAIN called with length < 0',ROUTCDE=11
C         ABEND X'30A000'
C         LTORG
C         END
CFile: convex.c
C/* c representation of fortran data types*/
C#define INTEGER          int
C#define DOUBLE_PRECISION double
C/* c representation of fortran routines */
C#define SECOND           second_
C#define GMAINV           gmainv_
C#define FMAIN            fmain_
C/* c representation of fortran common block name */
C#define CORCTL           _corctl_
C/* include files */
C#include <time.h>
C#include <sys/resource.h>
C                                                                        
C/* DOUBLE PRECISION FUNCTION SECOND() -- return cpu time in seconds */
CDOUBLE_PRECISION SECOND()
C{ struct rusage ruse;
C  getrusage(RUSAGE_SELF,&ruse);
C  return ruse.ru_utime.tv_sec + 0.000001 * ruse.ru_utime.tv_usec;
C}
C                                                                        
C                                                                        
C/* SUBROUTINE GMAINV(Q,IBASE,N)
C   COMMON /CORCTL/ INTREL,ICORCT(8)
C   set INTREL=number of integers per real
C   allocate n double precision words; returns ibase such that q(ibase)
C   is the first word allocated
C   SUBROUTINE FMAIN(Q,N)
C   release allocated area at q */
C                                                                        
Cstatic DOUBLE_PRECISION *pointer;
C                                                                        
CGMAINV(q,ibase,n) DOUBLE_PRECISION *q; INTEGER *ibase, *n;
C{ extern struct { INTEGER intrel; INTEGER icorct[8]} CORCTL ;
C  CORCTL.intrel = sizeof(DOUBLE_PRECISION)/sizeof(INTEGER);
C  pointer=(DOUBLE_PRECISION *) calloc(*n+2,sizeof(DOUBLE_PRECISION));
C  *ibase=pointer-q+1;
C}
C                                                                        
CFMAIN(q,n) DOUBLE_PRECISION *q; INTEGER *n; {free(pointer);}
CFile: vax.f
C      SUBROUTINE GMAINV (Q,IBASE,LENN)
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      LOGICAL ALIGND
C      COMMON /TAPES/ INP,IOUT
C      COMMON /CORCTL/ INTREL,ICORCT(8)
C      DIMENSION Q(*)
CC...  INTREL MUST BE NUMBER OF INTEGERS PER REAL
C      INTREL = 2
CC
C      IF (LENN.LE.0) THEN
C        LENN = 0
C        IBASE = 1
C        RETURN
C      END IF
CC...  NB MUST PROTECT AGAINST POSSIBILITY OF Q NOT BEING ON 8BYTE BOUNDA
C      ALIGND = (%LOC(Q)/8)*8.EQ.%LOC(Q)
C      LENNN = LENN*8
C      IF (.NOT.ALIGND) LENNN = LENNN+8
C      IBASE = LIB$GET_VM(LENNN,IBAES)
C      IF (.NOT.IBASE) THEN
C        WRITE (IOUT,10)LENN
C10      FORMAT(//' *** REQUEST FOR',I12,
C     1      ' WORDS OF STORE HAS FAILED')
C        CALL LIB$SIGNAL(%VAL(IBASE))
C      END IF
C      IBASE = (IBAES-%LOC(Q)+7)/8+1
C      RETURN
C      ENTRY FMAIN (Q,LENN)
CC...  SHOULD CALL LENN WITH SAME VALUE AS CALLED GMAINV
C      IF (LENN.LE.0) RETURN
C      ALIGND = (%LOC(Q)/8)*8.EQ.%LOC(Q)
C      LENNN = LENN
C      IF (.NOT.ALIGND) LENNN = LENN+1
C      IBAES = ((%LOC(Q))/8)*8
C      IRET = LIB$FREE_VM(LENNN,IBAES)
C      IF (.NOT.IRET) THEN
C       WRITE (IOUT,'('' ERROR IN FREEING STORE'')')
C       CALL LIB$SIGNAL(%VAL(IRET))
C      END IF
C      RETURN
C      END
C      FUNCTION SECOND()
C      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      SAVE INIT
C      DATA INIT/0/
C      IF (INIT.NE.0) THEN
C        CALL LIB$STAT_TIMER(2,ICPU)
C      ELSE
C        CALL LIB$INIT_TIMER
C        INIT = 1
C        ICPU = 0
C      END IF
C      SECOND=DFLOAT(ICPU)*0.01
C      RETURN
C      END
