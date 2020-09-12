!  EXPLICITFVMMESH448.f90 
!
!  FUNCTIONS:
!  EXPLICITFVMMESH448 - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: EXPLICITFVMMESH448
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

! THIS PROGRAM IS SOLVING HEAT TRANSFER IN 2D BY FVM
! EXPLCIT METHOD
! BY MEHRDAD BABA HOSEIN POUR
! BEYOND SKY!

program Explcit_FVM
    
implicit none
    
integer :: NNODE,i,NELEM,ELEM,j,M,G,N
    
DOUBLE PRECISION :: DELT,ALPHA,H,RO,CP,TINF
    
double precision,allocatable,dimension(:):: X,Y,T,BOUNDARY,K,CONS
    
double precision,allocatable,dimension(:,:)::XT,YT,DELL
    
integer,allocatable,dimension(:):: ELEM_1,ELEM_2,ELEM_3,NEIGHBOUR_1,NEIGHBOUR_2,NEIGHBOUR_3,NODE
    
INTEGER,allocatable,dimension(:,:)::NEIGHBOUR
    
NNODE=488
    
NELEM=894

H=28

RO=7800

CP=480

TINF=0
    
allocate (X(NNODE),Y(NNODE),NODE(NNODE),T(NNODE),BOUNDARY(NELEM),ELEM_1(NELEM),ELEM_2(NELEM),ELEM_3(NELEM),XT(3,NELEM),YT(3,NELEM))
    
allocate (NEIGHBOUR(NNODE,7*2),DELL(NNODE,7*2),K(NNODE),CONS(NNODE))
    
Open(1,file='Mesh 488.txt')
    
Do i=1,NNODE
  
	Read(1,*)NODE(i),X(i),Y(i)
    
End do
    
do i=NNODE+1,NELEM+NNODE
  
  	Read(1,*)ELEM,ELEM_1(ELEM),ELEM_2(ELEM),ELEM_3(ELEM)
    
    end do
    
    Do i=1,NNODE
      
	T(i)=50
    
End do

ALPHA=3.5

DELT=0.00000001

Do i=1,NNODE
  
	BOUNDARY(i)=0
    
End do

    do i= 1 ,NELEM

      do j=1,3

        if (j==1)then

      XT(j,i)=X(ELEM_1(i))
      
      YT(j,i)=Y(ELEM_1(i))

      elseif (j==2)then
      
      XT(j,i)=X(ELEM_2(i))
      
      YT(j,i)=Y(ELEM_2(i))
      
      elseif (j==3)then
      
      XT(j,i)=X(ELEM_3(i))
      
      YT(j,i)=Y(ELEM_3(i))

      end if

      end do

      end do

DO N=1,100
        
DO I= 1,NNODE

  M=1

      DO J =1,NELEM

        IF (ELEM_1(J)==NODE(I)) THEN

          NEIGHBOUR(I,M)=ELEM_2(J)

          M=M+1

          NEIGHBOUR(I,M)=ELEM_3(J)

          M=M+1

          ELSEIF (ELEM_2(J)==NODE(I)) THEN

          NEIGHBOUR(I,M)=ELEM_1(J)

          M=M+1

          NEIGHBOUR(I,M)=ELEM_3(J)

          M=M+1

          ELSEIF (ELEM_3(J)==NODE(I)) THEN


          NEIGHBOUR(I,M)=ELEM_1(J)

          M=M+1

          NEIGHBOUR(I,M)=ELEM_2(J)
          
          M=M+1

           END IF

            if ((XT(1,J)==0.AND.XT(2,J)==0))THEN

          BOUNDARY(J)=1
      
      T(ELEM_1(J))=T(ELEM_1(J))+(((H*DELT)/(RO*CP*ABS(YT(1,J)-YT(2,J))))*(T(ELEM_1(J))-TINF))
      
      T(ELEM_2(J))=T(ELEM_2(J))+(((H*DELT)/(RO*CP*ABS(YT(1,J)-YT(2,J))))*(T(ELEM_2(J))-TINF))
          
        ELSEIF((XT(2,J)==0.AND.XT(3,J)==0))THEN

        BOUNDARY(J)=1
      
      T(ELEM_2(J))=T(ELEM_2(J))+(((H*DELT)/(RO*CP*ABS(YT(3,J)-YT(2,J))))*(T(ELEM_2(J))-TINF))
      
      T(ELEM_3(J))=T(ELEM_3(J))+(((H*DELT)/(RO*CP*ABS(YT(3,J)-YT(2,J))))*(T(ELEM_3(J))-TINF))
        
       ELSEIF((XT(1,J)==0.AND.XT(3,J)==0))then

      BOUNDARY(J)=1
      
      T(ELEM_1(J))=T(ELEM_1(J))+(((H*DELT)/(RO*CP*ABS(YT(1,J)-YT(3,J))))*(T(ELEM_1(J))-TINF))
      
      T(ELEM_3(J))=T(ELEM_3(J))+(((H*DELT)/(RO*CP*ABS(YT(3,J)-YT(1,J))))*(T(ELEM_3(J))-TINF))
      
      elseif ((YT(1,J)==0.AND.YT(2,J)==0))THEN

      BOUNDARY(J)=2

     T(ELEM_1(J))=100
     
     T(ELEM_2(J))=100
      
      ELSEIF((YT(2,J)==0.AND.YT(3,J)==0))THEN

      BOUNDARY(J)=2
     
     T(ELEM_2(J))=100
     
     T(ELEM_3(J))=100
      
      ELSEIF((YT(1,J)==0.AND.YT(3,J)==0))then

     BOUNDARY(J)=2

     T(ELEM_1(J))=100
     
     T(ELEM_3(J))=100
      
      elseif ((XT(1,J)==1.AND.XT(2,J)==1))THEN

      BOUNDARY(J)=3

     T(ELEM_1(J))=50
     
     T(ELEM_2(J))=50
      
      ELSEIF((XT(2,J)==1.AND.XT(3,J)==1))THEN

      BOUNDARY(J)=3
     
     T(ELEM_2(J))=50
     
     T(ELEM_3(J))=50
      
      ELSEIF((XT(1,J)==1.AND.XT(3,J)==1))then
      
     BOUNDARY(J)=3

     T(ELEM_1(J))=50
     
     T(ELEM_3(J))=50
      
      elseif ((YT(1,J)==1.AND.YT(2,J)==1))THEN

      BOUNDARY(J)=4

      T(ELEM_1(J))=((1-2*(ALPHA*DELT/(((XT(1,J)-XT(2,J))**2))))*T(ELEM_1(J)))+&
      
      &((ALPHA*DELT/(((XT(1,J)-XT(2,J))**2)))*(T(ELEM_1(J))+T(ELEM_2(J))))
      
      T(ELEM_2(J))=((1-2*(ALPHA*DELT/(((XT(1,J)-XT(2,J))**2))))*T(ELEM_2(J)))+&
      
      &((ALPHA*DELT/(((XT(1,J)-XT(2,J))**2)))*(T(ELEM_1(J))+T(ELEM_2(J))))

      ELSEIF((YT(2,J)==1.AND.YT(3,J)==1))THEN

      BOUNDARY(J)=4
      
      T(ELEM_2(J))=((1-2*(ALPHA*DELT/(((XT(1,J)-XT(2,J))**2))))*T(ELEM_2(J)))+&
      
      &((ALPHA*DELT/(((XT(3,J)-XT(2,J))**2)))*(T(ELEM_3(J))+T(ELEM_2(J))))
      
      T(ELEM_3(J))=((1-2*(ALPHA*DELT/(((XT(1,J)-XT(2,J))**2))))*T(ELEM_3(J)))+&
      
      &((ALPHA*DELT/(((XT(3,J)-XT(2,J))**2)))*(T(ELEM_3(J))+T(ELEM_2(J))))
      
      ELSEIF((YT(1,J)==1.AND.YT(3,J)==1))then
      
      BOUNDARY(J)=4

      T(ELEM_1(J))=((1-2*(ALPHA*DELT/(((XT(1,J)-XT(2,J))**2))))*T(ELEM_1(J)))+&
      
      &((ALPHA*DELT/(((XT(3,J)-XT(1,J))**2)))*(T(ELEM_3(J))+T(ELEM_1(J))))
      
      T(ELEM_3(J))=((1-2*(ALPHA*DELT/(((XT(1,J)-XT(2,J))**2))))*T(ELEM_3(J)))+&
      
      &((ALPHA*DELT/(((XT(3,J)-XT(1,J))**2)))*(T(ELEM_3(J))+T(ELEM_1(J))))

      !ELSE
      
      end if


          END DO

          K(I)=0

          CONS(I)=(M-1)

          DO G=1,M-1

          DELL(I,G)=SQRT(((X(I)-X(NEIGHBOUR(I,G)))**2)+((Y(I)-Y(NEIGHBOUR(I,G)))**2))

          K(I)=K(I)+((ALPHA*DELT/(DELL(I,G)**2))*T(NEIGHBOUR(I,G)))

          END DO

          K(I)=K(I)/2

          T(I)=((1-(CONS(I)*K(I)))*T(I))+K(I)

          END DO
          
END DO

OPEN(UNIT=15,FILE='MESH488_EXPLICIT.PLT')

        WRITE(15,*) 'VARIABLES="X","Y","T"'
        
        WRITE(15,*)"Zone T='Coordinate.'"
        
write(15,*)"I=",NNODE, "J=",NELEM,",F=FEPOINT ET=Triangle"
        
do i=1,NNODE
  
	write(15,*)X(i),Y(i),T(i)
    
end do
    
do J=1,NELEM
  
	write(15,*)ELEM_1(J),ELEM_2(J),ELEM_3(J)
    
end do

close(15)

deallocate (X,Y,T,BOUNDARY,ELEM_1,ELEM_2,ELEM_3,XT,YT,NEIGHBOUR)

end program  Explcit_FVM

