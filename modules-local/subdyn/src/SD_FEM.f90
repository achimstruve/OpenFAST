!..................................................................................................................................
! LICENSING
! Copyright (C) 2013-2016  National Renewable Energy Laboratory
!
!    This file is part of SubDyn.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
MODULE SD_FEM
  USE NWTC_Library
  USE NWTC_LAPACK
  USE SubDyn_Types
  
  IMPLICIT NONE
  
   INTEGER,         PARAMETER  :: LAKi            = R8Ki                  ! Define the kind to be used for LAPACK routines for getting eigenvalues/vectors. Apparently there is a problem with SGGEV's eigenvectors
  
  INTEGER(IntKi),   PARAMETER  :: MaxMemjnt       = 10                    ! Maximum number of members at one joint
  INTEGER(IntKi),   PARAMETER  :: MaxOutChs       = 2000                  ! Max number of Output Channels to be read in
  INTEGER(IntKi),   PARAMETER  :: TPdofL          = 6                     ! 6 degrees of freedom (length of u subarray [UTP])
   
  ! values of these parameters are ordered by their place in SubDyn input file:
  INTEGER(IntKi),   PARAMETER  :: JointsCol       = 4                     ! Number of columns in Joints (JointID, JointXss, JointYss, JointZss)
  INTEGER(IntKi),   PARAMETER  :: ReactCol        = 7                     ! Number of columns in reaction dof array (JointID,RctTDxss,RctTDYss,RctTDZss,RctRDXss,RctRDYss,RctRDZss)
  INTEGER(IntKi),   PARAMETER  :: InterfCol       = 7                     ! Number of columns in interf matrix (JointID,ItfTDxss,ItfTDYss,ItfTDZss,ItfRDXss,ItfRDYss,ItfRDZss)
  INTEGER(IntKi),   PARAMETER  :: MaxNodesPerElem = 2                     ! Maximum number of nodes per element (currently 2)
  INTEGER(IntKi),   PARAMETER  :: MembersCol      = MaxNodesPerElem + 8   ! Number of columns in Members (Member-ID,MJointID1,MJointID2,MPropSetID1,MPropSetID2,Orientation,RotAngle,PointAXss,PointAYss,PointAZss)
  INTEGER(IntKi),   PARAMETER  :: PropSetsCol     = 6                     ! Number of columns in PropSets  (PropSetID,YoungE,ShearG,MatDens,XsecD,XsecT)  !bjj: this really doesn't need to store k, does it? or is this supposed to be an ID, in which case we shouldn't be storing k (except new property sets), we should be storing IDs
  INTEGER(IntKi),   PARAMETER  :: XPropSetsCol    = 14                    ! Number of columns in XPropSets (PropSetID,YoungE,ShearG,MatDens,XsecA,Xsecalphaxx,Xsecalphayy,Xsecalphaxy,Xsecalphazx,Xsecalphazy,XsecJxx,XsecJyy,XsecJxy,XsecJ0)
  INTEGER(IntKi),   PARAMETER  :: COSMsCol        = 10                    ! Number of columns in (cosine matrices) COSMs (COSMID,COSM11,COSM12,COSM13,COSM21,COSM22,COSM23,COSM31,COSM32,COSM33)
  INTEGER(IntKi),   PARAMETER  :: CMassCol        = 5                     ! Number of columns in Concentrated Mass (CMJointID,JMass,JMXX,JMYY,JMZZ)
  
  INTEGER(IntKi),   PARAMETER  :: SDMaxInpCols    = MAX(JointsCol,ReactCol,InterfCol,MembersCol,PropSetsCol,XPropSetsCol,COSMsCol,CMassCol)
    CONTAINS
    
SUBROUTINE NodeCon(Init,p, ErrStat, ErrMsg)
  
!This Subroutine maps nodes to elements
! allocate for NodesConnE and NodesConnN                                                                               
  USE qsort_c_module
  IMPLICIT NONE

  TYPE(SD_InitType),              INTENT( INOUT )  ::Init   
  TYPE(SD_ParameterType),         INTENT( IN    )  ::p  
  INTEGER(IntKi),                 INTENT(   OUT )  :: ErrStat     ! Error status of the operation
  CHARACTER(*),                   INTENT(   OUT )  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
  !local variable
  INTEGER(IntKi) :: SortA(MaxMemJnt,1)  !To sort nodes and elements
  INTEGER(IntKi) :: I,J,K  !counter
  
  ! bjj: shouldn't there be a check that there AREN'T actually more members at one joint than MaxMemJnt?
  ALLOCATE(Init%NodesConnE(Init%NNode, MaxMemJnt+1), STAT=ErrStat)                !the row index is the number of the real node, i.e. ID, 1st col has number of elements attached to node, and 2nd col has element numbers (up to 10)                                    
  IF ( ErrStat /= 0 )  THEN                                                                                                
      ErrMsg = ' Error allocating NodesConnE matrix'                                                    
      ErrStat = ErrID_Fatal                                                                                                         
      RETURN                                                                                                              
  ENDIF                                                                                                                  
  Init%NodesConnE = 0                                                                                                    
                                                                                                                          
  ALLOCATE(Init%NodesConnN(Init%NNode, MaxMemJnt+2), STAT=ErrStat)                                                    
  IF ( ErrStat /= 0 )  THEN                                                                                                
      ErrMsg = ' Error allocating NodesConnN matrix'                                                    
      ErrStat = ErrID_Fatal                                                                                                         
      RETURN                                                                                                              
  ENDIF                                                                                                                  
  Init%NodesConnN = 0                                                                                                    
                                                                                                                          
!   ! find the node connectivity, nodes/elements that connect to a common node                                             
   
   DO I = 1, Init%NNode                                                                                                   
      Init%NodesConnN(I, 1) = NINT( Init%Nodes(I, 1) )      !This should not be needed, could remove the extra 1st column like for the other array                                                                      
                                                                                                                          
      k = 0                                                                                                               
      DO J = 1, Init%NElem                          !This should be vectorized                                                                      
         IF ( ( NINT(Init%Nodes(I, 1))==p%Elems(J, 2)) .OR. (NINT(Init%Nodes(I, 1))==p%Elems(J, 3) ) ) THEN   !If i-th nodeID matches 1st node or 2nd of j-th element                                                                   
            k = k + 1                                                                                                     
            Init%NodesConnE(I, k + 1) = p%Elems(J, 1)                                                                  
            Init%NodesConnN(I, k + 1) = p%Elems(J, 3)                                                                  
            IF ( NINT(Init%Nodes(I, 1))==p%Elems(J, 3) ) Init%NodesConnN(I, k + 1) = p%Elems(J, 2)     !If nodeID matches 2nd node of element                                                                
         ENDIF                                                                                                            
      ENDDO                                                                                                               
                                                                                                                          
      IF( k>1 )THEN ! sort the nodes ascendingly                                                                          
         SortA(1:k, 1) = Init%NodesConnN(I, 3:(k+2))  
         CALL QsortC( SortA(1:k, 1:1) )                                                                                   
         Init%NodesConnN(I, 3:(k+2)) = SortA(1:k, 1)                                                                      
      ENDIF                                                                                                               
                                                                                                                          
      Init%NodesConnE(I, 1) = k    !Store how many elements connect i-th node in 2nd column                                                                                       
      Init%NodesConnN(I, 2) = k                                                                                           
   ENDDO                            

END SUBROUTINE NodeCon

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE SD_Discrt(Init,p, ErrStat, ErrMsg)

   TYPE(SD_InitType),            INTENT(INOUT)  ::Init
   TYPE(SD_ParameterType),       INTENT(INOUT)  ::p
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
      ! local variable
   INTEGER                       :: I, J, n, Node, Node1, Node2, Prop, Prop1, Prop2   
   INTEGER                       :: OldJointIndex(Init%NJoints)
   INTEGER                       :: NNE      ! number of nodes per element
   INTEGER                       :: MaxNXProp
   REAL(ReKi), ALLOCATABLE       :: TempProps(:, :)
   INTEGER, ALLOCATABLE          :: TempMembers(:, :) ,TempReacts(:,:)         
   INTEGER                       :: knode, kelem, kprop, nprop, MID
   REAL(ReKi)                    :: x1, y1, z1, x2, y2, z2, dx, dy, dz, A1, A2, dA, axx1, axx2, daxx, ayy1, ayy2, dayy, axy1, axy2, daxy, azx1, azx2, dazx, azy1, azy2, dazy, Ixx1, Ixx2, dIxx, Iyy1, Iyy2, dIyy, Ixy1, Ixy2, dIxy, Jzz1, Jzz2, dJzz
   LOGICAL                       :: found, CreateNewProp
   INTEGER(IntKi)                :: ErrStat2
   CHARACTER(1024)               :: ErrMsg2
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   ! number of nodes per element
   IF( ( Init%FEMMod .GE. 0 ) .and. (Init%FEMMod .LE. 3) ) THEN
      NNE = 2 
   ELSE
      CALL SetErrStat(ErrID_Fatal, 'FEMMod '//TRIM(Num2LStr(Init%FEMMod))//' not implemented.',ErrStat,ErrMsg,'SD_Discrt')
      RETURN
   ENDIF
   
   
   Init%NNode = Init%NJoints + ( Init%NDiv - 1 )*p%NMembers    ! Calculate total number of nodes according to divisions 
   Init%NElem = p%NMembers*Init%NDiv                           ! Total number of element   
   MaxNXProp   = Init%NPropSets + Init%NXPropSets + Init%NElem*NNE                ! Maximum possible number of property sets (temp): This is property set per element node, for all elements (bjj, added Init%NPropSets to account for possibility of entering many unused prop sets)(bas, added Init%NXPropSets to account for possibility of entering many unused xprop sets)
   
   ! Calculate total number of nodes and elements according to element types
   ! for 3-node or 4-node beam elements
   Init%NNode = Init%NNode + (NNE - 2)*Init%NElem
   !bjj: replaced with max value instead of NNE: Init%MembersCol = Init%MembersCol + (NNE - 2) 
   
   ! check the number of interior modes
   IF ( p%Nmodes .GT. 6*(Init%NNode - Init%NInterf - p%NReact) ) THEN
      CALL SetErrStat(ErrID_Fatal, ' NModes must be less than or equal to '//TRIM(Num2LStr( 6*(Init%NNode - Init%NInterf - p%NReact) )),ErrStat,ErrMsg,'SD_Discrt')
      RETURN
   ENDIF
   
   CALL AllocAry(p%Elems,         Init%NElem,    MembersCol, 'p%Elems',         ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_Discrt')
   
   CALL AllocAry(Init%Nodes,      Init%NNode,    JointsCol,  'Init%Nodes',      ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_Discrt')
   CALL AllocAry(Init%MemberNodes,p%NMembers,    Init%NDiv+1,'Init%MemberNodes',ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_Discrt')  ! for two-node element only, otherwise the number of nodes in one element is different
   CALL AllocAry(Init%MemberElements,p%NMembers,    Init%NDiv+1,'Init%MemberElements',ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_Discrt') 
   CALL AllocAry(Init%BCs,        6*p%NReact,    2,          'Init%BCs',        ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_Discrt') !!!! RRD: THIS MAY NEED TO CHANGE IF NOT ALL NODES ARE RESTRAINED
   CALL AllocAry(Init%IntFc,      6*Init%NInterf,2,          'Init%IntFc',      ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_Discrt')
   
   CALL AllocAry(TempMembers,     p%NMembers,    MembersCol, 'TempMembers',     ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_Discrt') 
   CALL AllocAry(TempProps,       MaxNXProp,      XPropSetsCol,'TempProps',       ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_Discrt') 
   CALL AllocAry(TempReacts,      p%NReact,      ReactCol,   'TempReacts',      ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_Discrt')
   

   IF ( ErrStat >= AbortErrLev ) THEN
      CALL CleanUp_Discrt()
      RETURN
   ENDIF

   ! Initialize Nodes
   Init%Nodes = 0   
   DO I = 1,Init%NJoints
      OldJointIndex(I) = Init%Joints(I, 1)
      Init%Nodes(I, 1) = I
      Init%Nodes(I, 2) = Init%Joints(I, 2)
      Init%Nodes(I, 3) = Init%Joints(I, 3)
      Init%Nodes(I, 4) = Init%Joints(I, 4)
   ENDDO
   
   ! Initialize Elems, starting with each member as an element (we'll take NDiv into account later)
   p%Elems = 0
   DO I = 1, p%NMembers
      p%Elems(I,     1) = I                     ! element/member number (not MemberID)
      p%Elems(I,     6) = Init%Members(I, 1)    !bas, adds the MemberID to each element in order to get the orientation information for each element later
!bjj: TODO: JMJ wants check that YoungE, ShearG, and MatDens are equal in the two properties because we aren't going to interpolate them. This should be less confusing for users.                                                
      
      
         ! loop through the JointIDs for this member and find the corresponding indices into the Joints array
      DO n = 2,3  ! Members column for JointIDs for nodes 1 and 2
         Node = Init%Members(I, n)  ! n=2 or 3
      
            ! ...... search for index of joint whose JointID matches Node ......
         J = 1
         found = .false.      
         DO WHILE ( .NOT. found .AND. J <= Init%NJoints )
            IF ( Node == NINT(Init%Joints(J, 1)) ) THEN
               p%Elems(I, n) = J                ! index of the joint/node n-1 (i.e., nodes 1 and 2)
               found = .TRUE.
            END IF
            J = J + 1
         END DO 
         
         IF ( .NOT. found) THEN
            CALL SetErrStat(ErrID_Fatal,' Member '//TRIM(Num2LStr(I))//' has JointID'//TRIM(Num2LStr(n-1))//' = '//& 
                                   TRIM(Num2LStr(Node))//' which is not in the node list !', ErrStat,ErrMsg,'SD_Discrt');
            CALL CleanUp_Discrt()
            RETURN
         END IF
         
      END DO ! loop through nodes/joints
      
      
         ! loop through the PropSetIDs for this member and find the corresponding indices into the Joints array
      
      ! we're setting these two values:   
      !p%Elems(I, 4) = property set for node 1 (note this sets the YoungE, ShearG, and MatDens columns for the ENTIRE element)   
      !p%Elems(I, 5) = property set for node 2 (note this should be used only for the XsecD and XsecT properties in the element [for a linear distribution from node 1 to node 2 of D and T])
         
      DO n=4,5 ! Member column for MPropSetID1 and MPropSetID2
      
                                                            
         Prop = Init%Members(I, n)  ! n=4 or 5
      
            ! ...... search for index of property set whose PropSetID matches Prop ......
         J = 1
         found = .false.      
         DO WHILE ( .NOT. found .AND. J <= Init%NPropSets )
            IF ( Prop == NINT(Init%PropSets(J, 1)) ) THEN
               p%Elems(I, n) = Init%NXPropSets + J  ! index of the property set n-3 (i.e., property sets 1 and 2)  ! note that previously, this used Prop instead of J, which assumed the list of MemberIDs was sequential, starting at 1.
               found = .TRUE.
            END IF
            J = J + 1
         END DO

         !bas, search also in XPropSets for matching property sets. Note that double occurrence of PropSetIDs is impossible, because it has already been checked within the SubDyn:SD_Input subroutine.
         J = 1
         DO WHILE ( .NOT. found .AND. J <= Init%NXPropSets )
            IF ( Prop == NINT(Init%XPropSets(J, 1)) ) THEN
               p%Elems(I, n) = J                    ! index of the property set n-3 (i.e., property sets 1 and 2)
               found = .TRUE.
            END IF
            J = J + 1
         END DO
         
         IF ( .NOT. found) THEN
            CALL SetErrStat(ErrID_Fatal,' Member '//TRIM(Num2LStr(I))//' has PropSetID'//TRIM(Num2LStr(n-3))//' = '//& 
                                   TRIM(Num2LStr(Prop))//' which is not in the Member X-Section Property data!', ErrStat,ErrMsg,'SD_Discrt');
            CALL CleanUp_Discrt()
            RETURN
         END IF

      END DO ! loop through property ids         
   
   END DO ! loop through members
   
   ! Conversion from circular PropSets to arbitrary XPropSets
   CALL ConvertPropSets(Init)
   
   ! Initialize TempMembers
   TempMembers = p%Elems(1:p%NMembers,:)
   
   ! Initialize Temp property set, first user defined sets
   TempProps = 0
   TempProps(1:Init%NXPropSets, :) = Init%XPropSets   
   
   ! Initialize boundary constraint vector
   ! Change the node number
   
   
    !Allocate array that will be p%Reacts renumbered and ordered so that ID does not play a role, just ordinal position number will count -RRD
   Init%BCs = 0
   TempReacts=0 !INitialize -RRD
   DO I = 1, p%NReact
      Node1 = p%Reacts(I, 1);  !NODE ID
      TempReacts(I,2:ReactCol)=p%Reacts(I, 2:ReactCol)  !Assign all the appropriate fixity to the new Reacts array -RRD
      
      found = .false.
      DO J = 1, Init%NJoints
         IF ( Node1 == NINT(Init%Joints(J, 1)) ) THEN
            Node2 = J
            found = .true.
            TempReacts(I,1)=Node2      !New node ID for p!React  -RRD
            EXIT  !Exit J loop if node found -RRD
         ENDIF
      ENDDO
      
      IF (.not. found) THEN
         CALL SetErrStat(ErrID_Fatal,' React has node not in the node list !', ErrStat,ErrMsg,'SD_Discrt');
         CALL CleanUp_Discrt()
         RETURN
      ENDIF
      
      
      DO J = 1, 6
         Init%BCs( (I-1)*6+J, 1) = (Node2-1)*6+J;
         Init%BCs( (I-1)*6+J, 2) = p%Reacts(I, J+1);
      ENDDO
      
   ENDDO
   p%Reacts=TempReacts   !UPDATED REACTS
      
   ! Initialize interface constraint vector
   ! Change the node number

   Init%IntFc = 0
   DO I = 1, Init%NInterf
      Node1 = Init%Interf(I, 1);
      found = .false.
      DO J = 1, Init%NJoints
         IF ( Node1 == NINT(Init%Joints(J, 1)) ) THEN
            Node2 = J
            found = .true.
         ENDIF
      ENDDO
      
      IF (.not. found) THEN
         CALL SetErrStat(ErrID_Fatal,' Interf has node not in the node list !', ErrStat,ErrMsg,'SD_Discrt');
         CALL CleanUp_Discrt()
         RETURN
      ENDIF
      
      DO J = 1, 6
         Init%IntFc( (I-1)*6+J, 1) = (Node2-1)*6+J;
         Init%IntFc( (I-1)*6+J, 2) = Init%Interf(I, J+1);
      ENDDO
   ENDDO
  
   ! Change numbering in concentrated mass matrix
   DO I = 1, Init%NCMass
      Node1 = NINT( Init%CMass(I, 1) )
      DO J = 1, Init%NJoints
         IF ( Node1 == NINT(Init%Joints(J, 1)) ) THEN
            Init%CMass(I, 1) = J  !bjj: todo: check this. if there is no return after this is found, are we overwritting the value if Node1 == NINT(Init%Joints(J, 1)) is true for multiple Js?
         ENDIF
      ENDDO
   ENDDO

! discretize structure according to NDiv 

knode = Init%NJoints
kelem = 0
kprop = Init%NXPropSets
Init%MemberNodes = 0
Init%MemberElements = 0


IF (Init%NDiv .GT. 1) THEN
   DO I = 1, p%NMembers !the first p%NMembers rows of p%Elems contain the element information
      ! create new node
      Node1 = TempMembers(I, 2)
      Node2 = TempMembers(I, 3)
      
      IF ( Node1==Node2 ) THEN
         CALL SetErrStat(ErrID_Fatal,' Same starting and ending node in the member.', ErrStat,ErrMsg,'SD_Discrt');
         CALL CleanUp_Discrt()
         RETURN
      ENDIF
    
      
      
      Prop1 = TempMembers(I, 4)
      Prop2 = TempMembers(I, 5)
      
      Init%MemberNodes(I,           1) = Node1
      Init%MemberNodes(I, Init%NDiv+1) = Node2
      
      MID = Init%Members(I,1)
      Init%MemberElements(I, 1) = MID ! write current MemberID to Init%MemberElements
      
      IF  ( ( .not. EqualRealNos(TempProps(Prop1, 2),TempProps(Prop2, 2) ) ) &
       .OR. ( .not. EqualRealNos(TempProps(Prop1, 3),TempProps(Prop2, 3) ) ) &
       .OR. ( .not. EqualRealNos(TempProps(Prop1, 4),TempProps(Prop2, 4) ) ) )  THEN
      
         CALL SetErrStat(ErrID_Fatal,' Material E,G and rho in a member must be the same', ErrStat,ErrMsg,'SD_Discrt');
         CALL CleanUp_Discrt()
         RETURN
      ENDIF

      x1 = Init%Nodes(Node1, 2)
      y1 = Init%Nodes(Node1, 3)
      z1 = Init%Nodes(Node1, 4)

      x2 = Init%Nodes(Node2, 2)
      y2 = Init%Nodes(Node2, 3)
      z2 = Init%Nodes(Node2, 4)
      
      dx = ( x2 - x1 )/Init%NDiv
      dy = ( y2 - y1 )/Init%NDiv
      dz = ( z2 - z1 )/Init%NDiv
      
      A1 = TempProps(Prop1, 5)
      axx1 = TempProps(Prop1, 6)
      ayy1 = TempProps(Prop1, 7)
      axy1 = TempProps(Prop1, 8)
      azx1 = TempProps(Prop1, 9)
      azy1 = TempProps(Prop1, 10)
      Ixx1 = TempProps(Prop1, 11)
      Iyy1 = TempProps(Prop1, 12)
      Ixy1 = TempProps(Prop1, 13)
      Jzz1 = TempProps(Prop1, 14)

      A2 = TempProps(Prop2, 5)
      axx2 = TempProps(Prop2, 6)
      ayy2 = TempProps(Prop2, 7)
      axy2 = TempProps(Prop2, 8)
      azx2 = TempProps(Prop2, 9)
      azy2 = TempProps(Prop2, 10)
      Ixx2 = TempProps(Prop2, 11)
      Iyy2 = TempProps(Prop2, 12)
      Ixy2 = TempProps(Prop2, 13)
      Jzz2 = TempProps(Prop2, 14)
      
      dA = ( A2 - A1 )/Init%NDiv
      daxx = ( axx2 - axx1 )/Init%NDiv
      dayy = ( ayy2 - ayy1 )/Init%NDiv
      daxy = ( axy2 - axy1 )/Init%NDiv
      dazx = ( azx2 - azx1 )/Init%NDiv
      dazy = ( azy2 - azy1 )/Init%NDiv
      dIxx = ( Ixx2 - Ixx1 )/Init%NDiv
      dIyy = ( Iyy2 - Iyy1 )/Init%NDiv
      dIxy = ( Ixy2 - Ixy1 )/Init%NDiv
      dJzz = ( Jzz2 - Jzz1 )/Init%NDiv
      
         ! If dA and dAx and dAy and dIxx and dIyy and dJzz are 0, no interpolation is needed, and we can use the same property set for new nodes/elements. otherwise we'll have to create new properties for each new node
      CreateNewProp = .NOT. ( EqualRealNos( dA , 0.0_ReKi ) .AND. &
                              EqualRealNos( daxx , 0.0_ReKi ) .AND. &
                              EqualRealNos( dayy , 0.0_ReKi ) .AND. &
                              EqualRealNos( daxy , 0.0_ReKi ) .AND. &
                              EqualRealNos( dazx , 0.0_ReKi ) .AND. &
                              EqualRealNos( dazy , 0.0_ReKi ) .AND. &
                              EqualRealNos( dIxx , 0.0_ReKi ) .AND. & 
                              EqualRealNos( dIyy , 0.0_ReKi ) .AND. &
                              EqualRealNos( dIxy , 0.0_ReKi ) .AND. &
                              EqualRealNos( dJzz , 0.0_ReKi ) )  
      
      ! node connect to Node1
      knode = knode + 1
      Init%MemberNodes(I, 2) = knode
      CALL GetNewNode(knode, x1+dx, y1+dy, z1+dz, Init)
      
      
      IF ( CreateNewProp ) THEN   
           ! create a new property set 
           ! k, E, G, rho, A, axx, ayy, axy, azx, azy, Ixx, Iyy, Ixy, Ixy, Jzz, Init
           
           kprop = kprop + 1
           CALL GetNewXProp(kprop, TempProps(Prop1, 2), TempProps(Prop1, 3),&
                           TempProps(Prop1, 4), A1+dA, axx1+daxx, ayy1+dayy, axy1+daxy, azx1+dazx, azy1+dazy, Ixx1+dIxx, Iyy1+dIyy, Ixy1+dIxy, Jzz1+dJzz, TempProps)           
           kelem = kelem + 1
           CALL GetNewElem(kelem, Node1, knode, Prop1, kprop, MID, p)  
           nprop = kprop              
      ELSE
           kelem = kelem + 1
           CALL GetNewElem(kelem, Node1, knode, Prop1, Prop1, MID, p)                
           nprop = Prop1 
      ENDIF
      
      Init%MemberElements(I, 2) = kelem ! write first element number to Init%MemberElements for this member
      
      ! interior nodes
      
      DO J = 2, (Init%NDiv-1)
         knode = knode + 1
         Init%MemberNodes(I, J+1) = knode

         CALL GetNewNode(knode, x1 + J*dx, y1 + J*dy, z1 + J*dz, Init)
         
         IF ( CreateNewProp ) THEN   
              ! create a new property set 
              ! k, E, G, rho, A, axx, ayy, axy, azx, azy, Ixx, Iyy, Ixy, Jzz, Init
              
              kprop = kprop + 1
           CALL GetNewXProp(kprop, TempProps(Prop1, 2), TempProps(Prop1, 3),&
                           TempProps(Prop1, 4), A1 + J*dA, axx1 + J*daxx, ayy1 + J*dayy,&
                           axy1 + J*daxy, azx1 + J*dazx, azy1 + J*dazy, Ixx1 + J*dIxx, Iyy1 + J*dIyy, Ixy1 + J*dIxy, Jzz1 + J*dJzz, TempProps)          
              kelem = kelem + 1
              CALL GetNewElem(kelem, knode-1, knode, nprop, kprop, MID, p)
              nprop = kprop                
         ELSE
              kelem = kelem + 1
              CALL GetNewElem(kelem, knode-1, knode, nprop, nprop, MID, p)         
               
         ENDIF
         
         Init%MemberElements(I, J+1) = kelem ! write next element number to Init%MemberElements for this member
      ENDDO
      
      ! the element connect to Node2
      kelem = kelem + 1
      CALL GetNewElem(kelem, knode, Node2, nprop, Prop2, MID, p) 
      Init%MemberElements(I, Init%NDiv + 1) = kelem ! write last element number to Init%MemberElements for this member

   ENDDO ! loop over all members

ELSE ! NDiv = 1

   Init%MemberNodes(1:p%NMembers, 1:2) = p%Elems(1:Init%NElem, 2:3)
   Init%MemberElements(1:p%NMembers, 1) = p%Elems(1:Init%NElem, 6)
   Init%MemberElements(1:p%NMembers, 2) = p%Elems(1:Init%NElem, 1)

ENDIF ! if NDiv is greater than 1

! set the props in Init
Init%NProp = kprop
CALL AllocAry(Init%Props, Init%NProp, XPropSetsCol,  'Init%Props', ErrStat2, ErrMsg2); CALL SetErrStat(ErrStat2,ErrMsg2,ErrStat,ErrMsg,'SD_Discrt')
   IF (ErrStat >= AbortErrLev ) THEN
      CALL SetErrStat(ErrStat2,ErrMsg2, ErrStat,ErrMsg,'SD_Discrt');
      CALL CleanUp_Discrt()
   RETURN
ENDIF

Init%Props = TempProps(1:Init%NProp, :)

CALL CleanUp_Discrt()

RETURN
CONTAINS
!................
   SUBROUTINE CleanUp_Discrt()
   
! deallocate temp matrices
IF (ALLOCATED(TempProps)) DEALLOCATE(TempProps)
IF (ALLOCATED(TempMembers)) DEALLOCATE(TempMembers)
IF (ALLOCATED(TempReacts)) DEALLOCATE(TempReacts)

   END SUBROUTINE CleanUp_Discrt

END SUBROUTINE SD_Discrt
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE GetNewNode(k, x, y, z, Init)

   TYPE(SD_InitType),      INTENT(INOUT) :: Init
   
   INTEGER,                INTENT(IN)    :: k
   REAL(ReKi),             INTENT(IN)    :: x, y, z
   
   Init%Nodes(k, 1) = k
   Init%Nodes(k, 2) = x
   Init%Nodes(k, 3) = y
   Init%Nodes(k, 4) = z


END SUBROUTINE GetNewNode
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE GetNewElem(k, n1, n2, p1, p2, MID, p)

   INTEGER,                INTENT(IN   )   :: k
   INTEGER,                INTENT(IN   )   :: n1
   INTEGER,                INTENT(IN   )   :: n2
   INTEGER,                INTENT(IN   )   :: p1
   INTEGER,                INTENT(IN   )   :: p2
   INTEGER,                INTENT(IN   )   :: MID
   TYPE(SD_ParameterType), INTENT(INOUT)   :: p
  
   
   p%Elems(k, 1) = k
   p%Elems(k, 2) = n1
   p%Elems(k, 3) = n2
   p%Elems(k, 4) = p1
   p%Elems(k, 5) = p2
   p%Elems(k, 6) = MID

END SUBROUTINE GetNewElem
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE GetNewProp(k, E, G, rho, d, t, TempProps)
!RRD-modifying this routine: This routine intends to calculate new member properties in case NDIV>1 ; 1/23/14
   
   INTEGER   , INTENT(IN)   :: k
   REAL(ReKi), INTENT(IN)   :: E, G, rho, d, t
   REAL(ReKi), INTENT(INOUT):: TempProps(:, :)
   
   TempProps(k, 1) = k
   TempProps(k, 2) = E
   TempProps(k, 3) = G
   TempProps(k, 4) = rho
   TempProps(k, 5) = d
   TempProps(k, 6) = t

END SUBROUTINE GetNewProp
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE GetNewXProp(k, E, G, rho, A, axx, ayy, axy, azx, azy, Ixx, Iyy, Ixy, Jzz, TempProps)

   
   INTEGER   , INTENT(IN)   :: k
   REAL(ReKi), INTENT(IN)   :: E, G, rho, A, axx, ayy, axy, azx, azy, Ixx, Iyy, Ixy, Jzz
   REAL(ReKi), INTENT(INOUT):: TempProps(:, :)
   
   TempProps(k, 1) = k
   TempProps(k, 2) = E
   TempProps(k, 3) = G
   TempProps(k, 4) = rho
   TempProps(k, 5) = A
   TempProps(k, 6) = axx
   TempProps(k, 7) = ayy
   TempProps(k, 8) = axy
   TempProps(k, 9) = azx
   TempProps(k, 10) = azy
   TempProps(k, 11) = Ixx
   TempProps(k, 12) = Iyy
   TempProps(k, 13) = Ixy
   TempProps(k, 14) = Jzz

END SUBROUTINE GetNewXProp
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE ConvertPropSets(Init)
! This subroutine converts the circular properties to general cross-section properties and adds them to the XPropSet
   TYPE(SD_InitType),            INTENT(INOUT)  ::Init

   INTEGER                  :: I, J

   REAL(ReKi)               :: Da, t, E, G, rho ! properties of a circular section
   REAL(ReKi)               :: Di, Ra, Ri, Ixx, Iyy, Jzz, A, axx, ayy, kappa, nu, ratioSq ! conversion entities                       

   
   J = Init%NXPropSets + 1 ! start index for circular properties within XPropSet
   ! loop over all circular property sets
   DO I = 1, Init%NPropSets
             
      ! get the current circular properties
      E   = Init%PropSets(I, 2)
      G   = Init%PropSets(I, 3)
      rho = Init%PropSets(I, 4)
      Da  = Init%PropSets(I, 5)
      t  = Init%PropSets(I, 6)
               
      ! conversion from circular to general x properties
      Ra = Da / 2.0_ReKi
      Ri = Ra - t
      Di = Ri * 2
            
      A = Pi_D * (Ra*Ra - Ri*Ri)
      Ixx = 0.25 * Pi_D * (Ra**4-Ri**4)
      Iyy = Ixx
      Jzz = 2.0 * Ixx
            
      ! calculate kappa, which is for a circular cross-section in each direction the same
      IF( Init%FEMMod == 1 ) THEN ! uniform Euler-Bernoulli
          axx = 0
                     
      ELSEIF( Init%FEMMod == 3 ) THEN ! uniform Timoshenko
          ! kappa = 0.53            
               
          ! equation 13 (Steinboeck et al) in SubDyn Theory Manual 
          nu = E / (2.0_ReKi*G) - 1.0_ReKi
          ratioSq = ( Di / Da )**2
          kappa =   ( 6.0 * (1.0 + nu) **2 * (1.0 + ratioSq)**2 ) &
                  / ( ( 1.0 + ratioSq )**2 * ( 7.0 + 14.0*nu + 8.0*nu**2 ) + 4.0 * ratioSq * ( 5.0 + 10.0*nu + 4.0 *nu**2 ) )
          axx = 1 / kappa
      ENDIF
      
      ayy = axx
                  
      ! add circular properties to XPropSet
      Init%XPropSets(J,1) = Init%PropSets(I,1)    ! add circular cross section ID to XPropSet
      Init%XPropSets(J,2) = E                     ! add E-modulus to XPropSet
      Init%XPropSets(J,3) = G                     ! add shear modulus to XPropSet
      Init%XPropSets(J,4) = rho                   ! add density to XPropSet
      Init%XPropSets(J,5) = A                     ! add area to XPropSet
      Init%XPropSets(J,6) = axx                   ! add shear deformation coefficient along x to XPropSet
      Init%XPropSets(J,7) = ayy                   ! add shear deformation coefficient along y to XPropSet
      Init%XPropSets(J,8) = 0.0_ReKi              ! add shear deformation coefficient for coupling between x and y to XPropSet
      Init%XPropSets(J,9) = 0.0_ReKi              ! add torsion-shear coefficient for torsional coupling with respect to x to XPropSet
      Init%XPropSets(J,10) = 0.0_ReKi             ! add torsion-shear coefficient for torsional coupling with respect to y to XPropSet
      Init%XPropSets(J,11) = Ixx                  ! add second area moment of inertia around x axis to XPropSet
      Init%XPropSets(J,12) = Iyy                  ! add second area moment of inertia around y axis to XPropSet
      Init%XPropSets(J,13) = 0.0_ReKi             ! add second area moment of inertia for coupling between x and y axis to XPropSet
      Init%XPropSets(J,14) = Jzz                  ! add torsional moment of inertia to XPropSet
      
      J = J + 1 ! add one to circular cross section index within XPropSet
   ENDDO

   Init%NXPropSets = Init%NXPropSets + Init%NPropSets ! update the number of general circular cross sections to its new value
   
END SUBROUTINE ConvertPropSets
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE AssembleKM(Init,p, ErrStat, ErrMsg)

   TYPE(SD_InitType),            INTENT(INOUT)  ::Init
   TYPE(SD_ParameterType),       INTENT(INOUT)  ::p
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   
   INTEGER                  :: I, J, K, Jn, Kn
   
   INTEGER                  :: NNE        ! number of nodes in one element
   INTEGER                  :: N1, N2     ! starting node and ending node in the element
   INTEGER                  :: P1, P2     ! property set numbers for starting and ending nodes

   REAL(ReKi)               :: E, G, rho, A, axx, ayy, axy, azx, azy, Ixx, Iyy, Ixy, Jzz ! properties of a section
   REAL(ReKi)               :: A1, A2, axx1, axx2, ayy1, ayy2, axy1, axy2, azx1, azx2, azy1, azy2, Ixx1, Ixx2, Iyy1, Iyy2, Ixy1, Ixy2, Jzz1, Jzz2 ! properties of each node from one element
   REAL(ReKi)               :: X1, Y1, Z1, X2, Y2, Z2    ! coordinates of the nodes
   REAL(ReKi)               :: MID                       ! Current MemberID
   REAL(ReKi)               :: psi                       ! Orientation angle of the current cross-section
   REAL(ReKi)               :: DirCos(3, 3)              ! direction cosine matrices
   REAL(ReKi)               :: L                         ! length of the element
   LOGICAL                  :: shear
   REAL(ReKi), ALLOCATABLE  :: Ke(:,:), Me(:, :), FGe(:) ! element stiffness and mass matrices gravity force vector
   INTEGER, ALLOCATABLE     :: nn(:)                     ! node number in element 
   INTEGER                  :: r
   
   
   INTEGER(IntKi)           :: ErrStat2
   CHARACTER(1024)          :: ErrMsg2


   !bas, FEMMod query is located out of the DO loop for efficiency (according to old comment from jj)
   IF  (Init%FEMMod == 2) THEN ! tapered Euler-Bernoulli
      CALL SetErrStat ( ErrID_Fatal, 'FEMMod = 2 is not implemented.', ErrStat, ErrMsg, 'AssembleKM' )
      CALL CleanUp_AssembleKM()
      RETURN
         
   ELSEIF  (Init%FEMMod == 4) THEN ! tapered Timoshenko
      CALL SetErrStat ( ErrID_Fatal, 'FEMMod = 4 is not implemented.', ErrStat, ErrMsg, 'AssembleKM' )
      CALL CleanUp_AssembleKM()
      RETURN

   ELSEIF  (Init%FEMMod == 1) THEN ! uniform Euler-Bernoulli
      Shear = .false.

   ELSEIF  (Init%FEMMod == 3) THEN ! uniform Timoshenko
      Shear = .true.
         
   ELSE
      CALL SetErrStat ( ErrID_Fatal, 'FEMMod is not valid. Please choose from 1, 2, 3, and 4. ', ErrStat, ErrMsg, 'AssembleKM' )
      CALL CleanUp_AssembleKM()
      RETURN
         
   ENDIF
   
      ! for current application
   IF ( (Init%FEMMod .LE. 3) .and. (Init%FEMMod .GE. 0)) THEN
      NNE = 2
   ELSE
      ErrStat = ErrID_Fatal
      ErrMsg = 'Invalid FEMMod in AssembleKM'
      RETURN
   ENDIF                          
   
   ! total degrees of freedom of the system 
   Init%TDOF = 6*Init%NNode
   
         ! Assemble system stiffness and mass matrices with gravity force vector
   
   ALLOCATE( p%ElemProps(Init%NElem), STAT=ErrStat2)
      IF (ErrStat2 /= 0) THEN
         CALL SetErrStat ( ErrID_Fatal, 'Error allocating p%ElemProps', ErrStat, ErrMsg, 'AssembleKM' )
         CALL CleanUp_AssembleKM()
      RETURN
   ENDIF

   CALL AllocAry( Ke,     NNE*6,         NNE*6 , 'Ke',      ErrStat2, ErrMsg2); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AssembleKM' )    ! element stiffness matrix
   CALL AllocAry( Me,     NNE*6,         NNE*6 , 'Me',      ErrStat2, ErrMsg2); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AssembleKM' )    ! element mass matrix 
   CALL AllocAry( FGe,    NNE*6,                 'FGe',     ErrStat2, ErrMsg2); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AssembleKM' )    ! element gravity force vector 
   CALL AllocAry( nn,     NNE,                   'nn',      ErrStat2, ErrMsg2); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AssembleKM' )    ! node number in element array 
   
   CALL AllocAry( Init%K, Init%TDOF, Init%TDOF , 'Init%K',  ErrStat2, ErrMsg2); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AssembleKM' )    ! system stiffness matrix 
   CALL AllocAry( Init%m, Init%TDOF, Init%TDOF , 'Init%M',  ErrStat2, ErrMsg2); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AssembleKM' )    ! system mass matrix 
   CALL AllocAry( Init%FG,Init%TDOF,             'Init%FG', ErrStat2, ErrMsg2); CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AssembleKM' )    ! system gravity force vector 
    
   IF (ErrStat >= AbortErrLev) THEN
      CALL CleanUp_AssembleKM()
      RETURN
   ENDIF
   
   Init%K = 0.0_ReKi
   Init%M = 0.0_ReKi
   Init%FG = 0.0_ReKi 
   
      ! loop over all elements
   DO I = 1, Init%NElem
   
      DO J = 1, NNE
         NN(J) = p%Elems(I, J + 1)
      ENDDO
   
      N1 = p%Elems(I,       2)
      N2 = p%Elems(I, NNE + 1)
      
      P1 = p%Elems(I, NNE + 2)
      P2 = p%Elems(I, NNE + 3)
      
      ! search for current MemberID
      DO J = 1, p%NMembers
         DO K = 1, Init%NDiv
            IF ( Init%MemberElements(J, K) == I ) THEN
               MID = Init%MemberElements(J, 1)
            ENDIF
         
         ENDDO
      
      ENDDO
      
      
      E   = Init%Props(P1, 2)
      G   = Init%Props(P1, 3)
      rho = Init%Props(P1, 4)
      A1  = Init%Props(P1, 5)
      axx1  = Init%Props(P1, 6)
      ayy1  = Init%Props(P1, 7)
      axy1  = Init%Props(P1, 8)
      azx1  = Init%Props(P1, 9)
      azy1  = Init%Props(P1, 10)
      Ixx1  = Init%Props(P1, 11)
      Iyy1  = Init%Props(P1, 12)
      Ixy1  = Init%Props(P1, 13)
      Jzz1  = Init%Props(P1, 14)
      A2  = Init%Props(P2, 5)
      axx2  = Init%Props(P2, 6)
      ayy2  = Init%Props(P2, 7)
      axy2  = Init%Props(P2, 8)
      azx2  = Init%Props(P2, 9)
      azy2  = Init%Props(P2, 10)
      Ixx2  = Init%Props(P2, 11)
      Iyy2  = Init%Props(P2, 12)
      Ixy2  = Init%Props(P2, 13)
      Jzz2  = Init%Props(P2, 14)
      
      X1  = Init%Nodes(N1, 2)
      Y1  = Init%Nodes(N1, 3)
      Z1  = Init%Nodes(N1, 4)
      
      X2  = Init%Nodes(N2, 2)
      Y2  = Init%Nodes(N2, 3)
      Z2  = Init%Nodes(N2, 4)

      CALL GetDirCos(Init, p, X1, Y1, Z1, X2, Y2, Z2, MID, DirCos, L, psi, ErrStat2, ErrMsg2)
         CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'AssembleKM' )
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp_AssembleKM()
            RETURN
         END IF
         
      A = (A1 + A2) / 2.0_ReKi
      axx = (axx1 + axx2) / 2.0_ReKi
      ayy = (ayy1 + ayy2) / 2.0_ReKi
      axy = (axy1 + axy2) / 2.0_ReKi
      azx = (azx1 + azx2) / 2.0_ReKi
      azy = (azy1 + azy2) / 2.0_ReKi
      Ixx = (Ixx1 + Ixx2) / 2.0_ReKi
      Iyy = (Iyy1 + Iyy2) / 2.0_ReKi
      Ixy = (Ixy1 + Ixy2) / 2.0_ReKi
      Jzz = (Jzz1 + Jzz2) / 2.0_ReKi
         
         
      p%ElemProps(i)%Area = A
      p%ElemProps(i)%Length = L
      p%ElemProps(i)%Ixx = Ixx
      p%ElemProps(i)%Iyy = Iyy
      p%ElemProps(i)%Ixy = Ixy
      p%ElemProps(i)%Jzz = Jzz
      p%ElemProps(i)%Shear = Shear
      p%ElemProps(i)%axx = axx
      p%ElemProps(i)%ayy = ayy
      p%ElemProps(i)%axy = axy
      p%ElemProps(i)%azx = azx
      p%ElemProps(i)%azy = azy
      p%ElemProps(i)%YoungE = E
      p%ElemProps(i)%ShearG = G
      p%ElemProps(i)%Rho = rho
      p%ElemProps(i)%psi = psi
      p%ElemProps(i)%DirCos = DirCos
         
         
      CALL ElemK(A, L, Ixx, Iyy, Ixy, Jzz, Shear, axx, ayy, axy, azx, azy, E, G, DirCos, Ke, ErrStat2, ErrMsg2)
         CALL SetErrStat ( ErrStat2, ErrMsg2, ErrStat, ErrMsg, 'ElemK' )
         IF (ErrStat >= AbortErrLev) THEN
            CALL CleanUp_AssembleKM()
            RETURN
         END IF
      CALL ElemM(A, L, Ixx, Iyy, Jzz, rho, DirCos, Me)
      CALL ElemG(A, L, rho, DirCos, FGe, Init%g)                                                                                                                                                               

      
      ! assemble element matrices to global matrices
         
      DO J = 1, NNE
         jn = nn(j)
         
         Init%FG( (jn*6-5):(jn*6) ) = Init%FG( (jn*6-5):(jn*6) ) &
                                    + FGe( (J*6-5):(J*6) )
         
         DO K = 1, NNE
            kn = nn(k)
            
            Init%K( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) = Init%K( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) &
                                                  + Ke( (J*6-5):(J*6), (K*6-5):(K*6) )
                  
            Init%M( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) = Init%M( (jn*6-5):(jn*6), (kn*6-5):(kn*6) ) &
                                                  + Me( (J*6-5):(J*6), (K*6-5):(K*6) )
                     
                     
         ENDDO !K
                     
      ENDDO !J
                     
                     
   ENDDO ! I end loop over elements
               
      
      ! add concentrated mass 
   DO I = 1, Init%NCMass
      DO J = 1, 3
          r = ( NINT(Init%CMass(I, 1)) - 1 )*6 + J
          Init%M(r, r) = Init%M(r, r) + Init%CMass(I, 2)
          
      ENDDO
      DO J = 4, 6
          r = ( NINT(Init%CMass(I, 1)) - 1 )*6 + J
          Init%M(r, r) = Init%M(r, r) + Init%CMass(I, J-1)
      ENDDO

   ENDDO ! I concentrated mass
 
      ! add concentrated mass induced gravity force
   DO I = 1, Init%NCMass
      
      r = ( NINT(Init%CMass(I, 1)) - 1 )*6 + 3
      Init%FG(r) = Init%FG(r) - Init%CMass(I, 2)*Init%g 

   ENDDO ! I concentrated mass induced gravity
   
   CALL CleanUp_AssembleKM()
   RETURN
   
CONTAINS
!..............
   SUBROUTINE CleanUp_AssembleKM()
! deallocate temp matrices
      IF (ALLOCATED(Ke)) DEALLOCATE(Ke)
      IF (ALLOCATED(Me)) DEALLOCATE(Me)
      IF (ALLOCATED(FGe)) DEALLOCATE(FGe)
      IF (ALLOCATED(nn)) DEALLOCATE(nn)
   END SUBROUTINE CleanUp_AssembleKM
   
END SUBROUTINE AssembleKM
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE GetDirCos(Init, p, X1, Y1, Z1, X2, Y2, Z2, MID, DirCos, Le, psi, ErrStat, ErrMsg)
   !This should be from local to global -RRD
   ! bjj: note that this is the transpose of what is normally considered the Direction Cosine Matrix  
   !      in the FAST framework. It seems to be used consistantly in the code (i.e., the transpose 
   !      of this matrix is used later).
   
   TYPE(SD_InitType),            INTENT(IN   )  :: Init
   TYPE(SD_ParameterType),       INTENT(IN   )  :: p
   REAL(ReKi) ,                  INTENT(IN   )  :: MID                     ! Current MemberID
   REAL(ReKi) ,                  INTENT(IN   )  :: X1, Y1, Z1, X2, Y2, Z2  ! (x,y,z) positions of two nodes making up an element
   REAL(ReKi) ,                  INTENT(  OUT)  :: DirCos(3, 3)            ! calculated direction cosine matrix
   REAL(ReKi) ,                  INTENT(  OUT)  :: Le                      ! length of element
   REAL(ReKi) ,                  INTENT(  OUT)  :: psi                     ! Orientation angle according to information from member MID
   
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat                 ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg                  ! Error message if ErrStat /= ErrID_None
   
   REAL(ReKi)                       ::  Dx,Dy,Dz           ! distances between nodes
   INTEGER(IntKi)       :: I                               ! Counter index
   REAL(ReKi)           :: PS(3), PE(3), SES(3)            ! Start point PS, end point PE and its subtraction as vectors. Contains its (x,y,z) coordinates
   REAL(ReKi)           :: PA(3)                           ! 3rd orientation point PA of this member. Contains its (x,y,z) coordinates
   REAL(ReKi)           :: PAp(3), PA_ll(3)                ! Projected orientation point PAp and its local vector from start point PS for this member. Contains its (x,y,z) coordinates
   REAL(ReKi)           :: ie_hat(3), je_hat(3), ke_hat(3) ! Unit vector along x_e, y_e and z_e axis of the current member. Contains its (x,y,z) coordinates
   REAL(ReKi)           :: lambda                          ! Scalar for point projection
   REAL(ReKi)           :: Lexy                            ! Length from PS to PE (projected on global SS-XY-plane)
   REAL(ReKi)           :: O_type                          ! Orientation type
   INTEGER(IntKi)       :: MI                              ! Member index
        
   ErrMsg  = ""
   ErrStat = ErrID_None
   
   ! calculate point distances
   Dy=Y2-Y1
   Dx=X2-X1
   Dz=Z2-Z1
   Lexy = sqrt( Dx**2 + Dy**2 )
   Le   = sqrt( Dx**2 + Dy**2 + Dz**2)
   
   IF ( EqualRealNos(Le, 0.0_ReKi) ) THEN
      ErrMsg = ' Same starting and ending location in the element.'
      ErrStat = ErrID_Fatal
      RETURN
   ENDIF
   
   ! get the index of the current member regardless of their sequence in Init%Members array
   DO I = 1, p%NMembers
      IF ( EqualRealNos(Init%Members(I, 1), MID) ) THEN
         MI = I
      ENDIF
   ENDDO

   DirCos=0.0_ReKi              ! whole matrix set to 0
   O_type = Init%Members(MI, 6) ! get current orientation type (1 = direct psi specification; 2 = 3rd point A specification)

   IF ( EqualRealNos(O_type, 1.0_ReKi) ) THEN ! In the case of user specified angle psi
      psi = Init%Members(MI, 7) * Pi_D / 180.0_ReKi ! convert user defined angle from deg into rad
   
      IF ( EqualRealNos(Lexy, 0.0_ReKi) ) THEN
         IF ( Dz < 0.0_ReKi) THEN  !z is kept along global z
            DirCos(1, 1) =  COS(psi)
            DirCos(1, 2) =  -SIN(psi)
            DirCos(2, 1) =  -SIN(psi)
            DirCos(2, 2) =  -COS(psi)
            DirCos(3, 3) = -1.0_ReKi
         ELSE
            DirCos(1, 1) =  COS(psi)
            DirCos(1, 2) =  -SIN(psi)
            DirCos(2, 1) =  SIN(psi)
            DirCos(2, 2) =  COS(psi)
            DirCos(3, 3) =  1.0_ReKi
         ENDIF 
      ELSE
         DirCos(1, 1) = -Dy / Lexy * COS(psi) - ( ( Dx * Dz ) / ( Lexy * Le ) * SIN(psi) )
         DirCos(1, 2) = -( -Dy / Lexy * SIN(psi) ) - ( ( Dx * Dz ) / ( Lexy * Le ) * COS(psi) )
         DirCos(1, 3) =  Dx/Le
      
         DirCos(2, 1) =  Dx / Lexy * COS(psi) + ( -Dy * Dz ) / ( Lexy * Le ) * SIN(psi)
         DirCos(2, 2) = -( Dx / Lexy * SIN(psi) ) + ( -Dy * Dz ) / ( Lexy * Le ) * COS(psi)
         DirCos(2, 3) =  - (-Dy) / Le
     
         DirCos(3, 1) = Lexy / Le * SIN(psi)
         DirCos(3, 2) = Lexy / Le * COS(psi)
         DirCos(3, 3) = Dz / Le
      ENDIF
    
   ELSEIF ( EqualRealNos(O_type, 2.0_ReKi) ) THEN ! In the case of user defined 3rd point A
      
      !! point projection      
      ! get start, end and orientation point position vector
      PS(1) = X1
      PS(2) = Y1
      PS(3) = Z1
      PE(1) = X2
      PE(2) = Y2
      PE(3) = Z2
      PA(1) = Init%Members(MI, 8)
      PA(2) = Init%Members(MI, 9)
      PA(3) = Init%Members(MI, 10)
      
      ! calculate unit vector along z_e axis
      SES = PE - PS
      ke_hat = SES / SQRT( SES(1)**2.0_ReKi + SES(2)**2.0_ReKi + SES(3)**2.0_ReKi )
      
      
      !! check if point PA lies on the z_e axis. The cases make sure that we get no devision by 0.
      ! case member z_e axis is parallel to global SS-X-axis
      IF ( EqualRealNos(ke_hat(3), 0.0_ReKi) .AND. EqualRealNos(ke_hat(2), 0.0_ReKi) .AND. EqualRealNos(( PA(3) - PS(3) ), 0.0_ReKi) .AND. EqualRealNos(( PA(2) - PS(2) ), 0.0_ReKi)) THEN
         ErrMsg = ' Orientation point A is not allowed to lie on member axis z_e!'
         ErrStat = ErrID_Fatal
         RETURN
         
      ! case member z_e axis is parallel to global SS-Y-axis
      ELSEIF ( EqualRealNos(ke_hat(3), 0.0_ReKi) .AND. EqualRealNos(ke_hat(1), 0.0_ReKi) .AND. EqualRealNos(( PA(3) - PS(3) ), 0.0_ReKi) .AND. EqualRealNos(( PA(1) - PS(1) ), 0.0_ReKi)) THEN
         ErrMsg = ' Orientation point A is not allowed to lie on member axis z_e!'
         ErrStat = ErrID_Fatal
         RETURN

      ! case member z_e axis is parallel to global SS-Z-axis
      ELSEIF ( EqualRealNos(ke_hat(1), 0.0_ReKi) .AND. EqualRealNos(ke_hat(2), 0.0_ReKi) .AND. EqualRealNos(( PA(1) - PS(1) ), 0.0_ReKi) .AND. EqualRealNos(( PA(2) - PS(2) ), 0.0_ReKi)) THEN
         ErrMsg = ' Orientation point A is not allowed to lie on member axis z_e!'
         ErrStat = ErrID_Fatal
         RETURN
      
      ! case member z_e axis is not parallel to global SS-X-Y-Z axes (implies that neither ke_hat(1), ke_hat(2) and ke_hat(3) is 0)
      ELSEIF ( .NOT. EqualRealNos(ke_hat(1), 0.0_ReKi) .AND. .NOT. EqualRealNos(ke_hat(2), 0.0_ReKi) .AND. .NOT. EqualRealNos(ke_hat(3), 0.0_ReKi) ) THEN
         IF ( ( PA(1) - PS(1) ) / ke_hat(1) ==  ( PA(2) - PS(2) ) / ke_hat(2) .AND. &
              ( PA(1) - PS(1) ) / ke_hat(1) ==  ( PA(3) - PS(3) ) / ke_hat(3) ) THEN
            ErrMsg = ' Orientation point A is not allowed to lie on member axis z_e!'
            ErrStat = ErrID_Fatal
            RETURN
            
         ENDIF     

      ENDIF
      
      !! calculate projected point PAp for different cases
      ! case member z_e axis is parallel to global SS-X-axis
      IF ( EqualRealNos(ke_hat(2), 0.0_ReKi) .AND. EqualRealNos(ke_hat(3), 0.0_ReKi) ) THEN
          PAp(2) = PA(2)
          PAp(3) = PA(3)
          PAp(1) = PS(1)
          
      ! case member z_e axis is parallel to global SS-Y-axis
      ELSEIF ( EqualRealNos(ke_hat(1), 0.0_ReKi) .AND. EqualRealNos(ke_hat(3), 0.0_ReKi) ) THEN
          PAp(1) = PA(1)
          PAp(3) = PA(3)
          PAp(2) = PS(2)
          
      ! case member z_e axis is parallel to global SS-Z-axis
      ELSEIF ( EqualRealNos(ke_hat(1), 0.0_ReKi) .AND. EqualRealNos(ke_hat(2), 0.0_ReKi) ) THEN
          PAp(1) = PA(1)
          PAp(2) = PA(2)
          PAp(3) = PS(3)
          
      ! case member z_e axis is parallel to global SS-XY-plane
      ELSEIF ( EqualRealNos(ke_hat(3), 0.0_ReKi) .AND. .NOT. EqualRealNos(ke_hat(1), 0.0_ReKi) .AND. .NOT. EqualRealNos(ke_hat(2), 0.0_ReKi)) THEN
         ! calculate lambda
         lambda = ( ke_hat(1)**2.0_ReKi + ke_hat(2)**2.0_ReKi) / ( ke_hat(1) * PA(1) - ke_hat(1) * PS(1) &
                    + ke_hat(2) * PA(2) - ke_hat(2) * PS(2) )
         ! calculate projected point PAp
         PAp = PA - 1 / lambda * ke_hat
         
      ! case member z_e axis is not parallel to global SS-X-Y-Z axes (implies that neither ke_hat(1), ke_hat(2) and ke_hat(3) is 0)
      ELSE
         ! calculate lambda
         lambda = ( ke_hat(3) +  ke_hat(1)**2.0_ReKi / ke_hat(3) + ke_hat(2)**2.0_ReKi / ke_hat(3)) &
                 / ( PA(3) + ( PA(1) * ke_hat(1) ) / ke_hat(3) + ( PA(2) * ke_hat(2) ) / ke_hat(3) &
                 - PS(3) - ( PS(1) * ke_hat(1) ) / ke_hat(3) - ( PS(2) * ke_hat(2) ) / ke_hat(3) )
         ! calculate projected point PAp
         PAp = PA - 1 / lambda * ke_hat
      ENDIF
      
      !! calculate orientation angle psi according to PAp 
      PA_ll = PAp - PS                                                                             ! calculate the local direction vector PA_ll, which is parallel the the local x_e, y_e plane
      ie_hat = PA_ll / SQRT( PA_ll(1)**2.0_ReKi + PA_ll(2)**2.0_ReKi + PA_ll(3)**2.0_ReKi )            ! calculate unit vector along the local x_e axis
      je_hat = CROSS_PRODUCT(ke_hat, ie_hat)                                                       ! calculate unit vector along the local y_e axis

      ! write unit vectors to direction cosine matrix
      DirCos(1:3,1) = ie_hat(1:3)
      DirCos(1:3,2) = je_hat(1:3)
      DirCos(1:3,3) = ke_hat(1:3)
      
      ! It follows the calculation of psi for the output of each member orientation. This calculation is not relevant for the following program parts.
      IF (.NOT. EqualRealNos(je_hat(3), 0.0_ReKi) ) THEN
         psi = ATAN2( DirCos(3,1), DirCos(3,2) )
         
      ELSE
         IF (.NOT. EqualRealNos(ie_hat(2), 0.0_ReKi) ) THEN ! implies that ie_hat lies not parallel to the global X-Z plane
            IF ( ie_hat(2) > 0.0_ReKi ) THEN
               psi = Pi_D / 2.0_ReKi - ATAN( ie_hat(1) / ie_hat(2) )
               
            ELSE 
               psi = - Pi_D / 2.0_ReKi - ATAN( ie_hat(1) / ie_hat(2) )
               
            ENDIF
            IF ( ke_hat(3) < 0.0_ReKi ) THEN
               psi = psi * (-1)
               
            ENDIF
            
         ELSE ! implies that ie_hat lies parallel to the global X-Z plane, which leads to four possible angles psi configurations
            IF ( je_hat(2) > 0 ) THEN
               IF ( ke_hat(3) > 0 ) THEN
                  psi = 0
                  
               ELSE
                  psi = Pi_D
               
               ENDIF
               
            ELSE
               IF ( ke_hat(3) > 0 ) THEN
                  psi = Pi_D
                  
               ELSE
                  psi = 0
               
               ENDIF 
            
            ENDIF
            
         ENDIF
         
      ENDIF
      
   ELSE
      ErrMsg = ' Member '//TRIM(Num2LStr(MID))//' has the undifined Orientation Type '//TRIM(Num2LStr(O_type))//'. Choose Orientation Type 1 or 2 instead.'
      ErrStat = ErrID_Fatal
      RETURN
   
   ENDIF

END SUBROUTINE GetDirCos
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE ElemK(A, L, Ixx, Iyy, Ixy, Jzz, Shear, axx, ayy, axy, azx, azy, E, G, DirCos, K, ErrStat, ErrMsg)
   ! element stiffness matrix according to Schramm, U.; Rubenchik, V. and Pilkey W. D.: "Beam stiffness matrix based on the elasticity equations", International Journal for Numerical Methods in Engineering, Vol. 40, p. 211-232, 1997.
   ! shear is true  -- non-tapered Timoshenko beam 
   ! shear is false -- non-tapered Euler-Bernoulli beam 

   REAL(ReKi), INTENT( IN)               :: A, L, Ixx, Iyy, Ixy, Jzz, axx, ayy, axy, azx, azy, E, G
   REAL(ReKi), INTENT( IN)               :: DirCos(3,3)
   LOGICAL, INTENT( IN)                  :: Shear
   
   REAL(ReKi), INTENT(OUT)               :: K(12, 12)  !RRD:  Ke and Me  need to be modified if convention of dircos is not followed?
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat                 ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg                  ! Error message if ErrStat /= ErrID_None
   
   REAL(ReKi)                            :: A_bar(3, 3) ! contains the shear, torsional and torsional-shear coupling terms
   REAL(ReKi)                            :: S_bar(3, 3) ! contains the second area moments of inertia
   REAL(ReKi)                            :: I_bar(3, 3) ! modified identety matrix
   REAL(ReKi)                            :: D(3, 3)     ! auxiliary matrix D
   INTEGER(IntKi)                        :: IPIV(3)     ! The pivot indices; for 1 <= i <= min(M,N), row i of the matrix was interchanged with row IPIV(i)
   REAL(R8Ki)                            :: WORK(3)     ! WORK(1) returns the optimal LWORK
   REAL(R8Ki)                            :: D_inv(3, 3) ! inverted auxiliary matrix D
   REAL(ReKi)                            :: H(3, 3)     ! auxiliary matrix H
   REAL(ReKi)                            :: B(3, 3)     ! auxiliary matrix B (used instead of K, proposed in the paper)
   REAL(ReKi)                            :: AM1(3, 3)   ! auxiliary matrix
   REAL(ReKi)                            :: AM2(3, 3)   ! auxiliary matrix
   REAL(ReKi)                            :: AM3(3, 3)   ! auxiliary matrix
   REAL(ReKi)                            :: AM4(3, 3)   ! auxiliary matrix
   REAL(ReKi)                            :: DC(12, 12)  ! direction cosine matrix
   
   ErrMsg  = ""
   ErrStat = ErrID_None
   
   K = 0
   
   ! define I_bar
   I_bar = 0
   I_bar(2, 2) = 1.0_ReKi
   I_bar(3, 3) = 1.0_ReKi
   
   ! define S_bar
   S_bar = 0
   S_bar(1, 1) = 1
   S_bar(2, 2) = E * Ixx
   S_bar(2, 3) = E * Ixy
   S_bar(3, 2) = E * Ixy
   S_bar(3, 3) = E * Iyy
   
   ! define A_bar
   A_bar = 0
   A_bar(1,1) = 1 / (G * Jzz)
   A_bar(1,2) = azx / (G * Jzz)
   A_bar(1,3) = azy / (G * Jzz)
   A_bar(2,1) = azx / (G * Jzz)
   A_bar(3,1) = azy / (G * Jzz)
   A_bar(2,2) = axx / (G * A)
   A_bar(2,3) = axy / (G * A)
   A_bar(3,2) = axy / (G * A)
   A_bar(3,3) = ayy / (G * A)
   
   D = I_bar + 12 / L**2 * MATMUL(S_bar, A_bar)
   D_inv = D
   H = 4 * I_bar + 12 / L**2 * MATMUL(S_bar, A_bar)
   B = 2 * I_bar - 12 / L**2 * MATMUL(S_bar, A_bar)
   
   CALL LAPACK_dgetrf( 3, 3, D_inv, IPIV, ErrStat, ErrMsg ) ! On entry, the M-by-N matrix to be factored. On exit, the factors L and U from the factorization D_inv = P*L*U; the unit diagonal elements of L are not stored.
      IF ( ErrStat/= 0 ) THEN
         ErrStat = ErrID_Fatal
         RETURN
      END IF
    CALL LAPACK_dgetri( 3, D_inv, IPIV, WORK, 3, ErrStat, ErrMsg ) ! !< On entry, the factors L and U from the factorization D_inv = P*L*U as computed by DGETRF. On exit, if INFO = 0, the inverse of the original matrix D_inv.
      IF ( ErrStat/= 0 ) THEN
         ErrStat = ErrID_Fatal
         RETURN
      END IF
   
   AM1 = 12 / L**3 * MATMUL(D_inv, S_bar)
   AM2 = 6 / L**2 * MATMUL(D_inv, S_bar)
   AM3 = 1 / L * MATMUL(MATMUL(H, D_inv), S_bar)
   AM4 = 1 / L * MATMUL(MATMUL(B, D_inv), S_bar)
   
   K(1:2, 1:2) = AM1(2:3, 2:3)
   K(3, 3) = E * A / L
   K(1, 4) = - AM2(2, 3)
   K(1, 5) = AM2(2, 2)
   K(2, 4) = - AM2(3, 3)
   K(2, 5) = AM2(3, 2)
   K(1:2, 6) = AM1(2:3, 1)
   K(1:2, 7:8) = - AM1(2:3, 2:3)
   K(3, 9) = - E * A / L
   K(1, 10) = - AM2(2, 3)
   K(1, 11) = AM2(2, 2)
   K(2, 10) = - AM2(3, 3)
   K(2, 11) = AM2(3, 2)
   K(1:2, 12) = - AM1(2:3, 1)
   K(4, 1) = K(1, 4)
   K(4, 2) = K(2, 4)
   K(5, 1) = K(1, 5)
   K(5, 2) = K(2, 5)
   K(6, 1) = K(1, 6)
   K(6, 2) = K(2, 6)
   K(4, 4) = AM3(3, 3)
   K(4, 5) = - AM3(3, 2)
   K(4, 6) = - AM2(3, 1)
   K(4, 7) = AM2(3, 2)
   K(4, 8) = AM2(3, 3)
   K(4, 10) = AM4(3, 3)
   K(4, 11) = - AM4(3, 2)
   K(4, 12) = AM2(3, 1)
   K(5, 4) = K(4, 5)
   K(5, 5) = AM3(2, 2)
   K(5, 6) = AM2(2, 1)
   K(5, 7) = - AM2(2, 2)
   K(5, 8) = - AM2(2, 3)
   K(5, 10) = - AM4(2, 3)
   K(5, 11) = AM4(2, 2)
   K(5, 12) = - AM2(2, 1)
   K(6, 4) = K(4, 6)
   K(6, 5) = K(5, 6)
   K(6, 6) = AM1(1, 1)
   K(6, 7) = - AM1(1, 2)
   K(6, 8) = - AM1(1, 3)
   K(6, 10) = - AM2(1, 3)
   K(6, 11) = AM2(1, 2)
   K(6, 12) = - AM1(1, 1)
   K(7, 1) = K(1, 7)
   K(7, 2) = K(2, 7)
   K(7, 3) = K(3, 7)
   K(7, 4) = K(4, 7)
   K(7, 5) = K(5, 7)
   K(7, 6) = K(6, 7)
   K(7, 7) = AM1(2, 2)
   K(7, 8) = AM1(2, 3)
   K(7, 10) = AM2(2, 3)
   K(7, 11) = - AM2(2, 2)
   K(7, 12) = AM1(2, 1)
   K(8, 1) = K(1, 8)
   K(8, 2) = K(2, 8)
   K(8, 3) = K(3, 8)
   K(8, 4) = K(4, 8)
   K(8, 5) = K(5, 8)
   K(8, 6) = K(6, 8)
   K(8, 7) = K(7, 8)
   K(8, 8) = AM1(3, 3)
   K(8, 10) = AM2(3, 3)
   K(8, 11) = - AM2(3, 2)
   K(8, 12) = AM1(3, 1)
   K(9, 3) = - E * A / L
   K(9, 9) = E * A / L
   K(10, 1) = K(1, 10)
   K(10, 2) = K(2, 10)
   K(10, 3) = K(3, 10)
   K(10, 4) = K(4, 10)
   K(10, 5) = K(5, 10)
   K(10, 6) = K(6, 10)
   K(10, 7) = K(7, 10)
   K(10, 8) = K(8, 10)
   K(10, 9) = K(9, 10)
   K(10, 10) = AM3(3, 3)
   K(10, 11) = - AM3(3, 2)
   K(10, 12) = AM2(3, 1)
   K(11, 1) = K(1, 11)
   K(11, 2) = K(2, 11)
   K(11, 3) = K(3, 11)
   K(11, 4) = K(4, 11)
   K(11, 5) = K(5, 11)
   K(11, 6) = K(6, 11)
   K(11, 7) = K(7, 11)
   K(11, 8) = K(8, 11)
   K(11, 9) = K(9, 11)
   K(11, 10) = K(10, 11)
   K(11, 11) = AM3(2, 2)
   K(11, 12) = - AM2(2, 1)
   K(12, 1) = K(1, 12)
   K(12, 2) = K(2, 12)
   K(12, 3) = K(3, 12)
   K(12, 4) = K(4, 12)
   K(12, 5) = K(5, 12)
   K(12, 6) = K(6, 12)
   K(12, 7) = K(7, 12)
   K(12, 8) = K(8, 12)
   K(12, 9) = K(9, 12)
   K(12, 10) = K(10, 12)
   K(12, 11) = K(11, 12)
   K(12, 12) = AM1(1, 1)
   
   DC = 0
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   K = MATMUL( MATMUL(DC, K), TRANSPOSE(DC) )
   
   !write(*, *) K - TRANSPOSE(K)

END SUBROUTINE ElemK
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE ElemM(A, L, Ixx, Iyy, Jzz, rho, DirCos, M)
   ! element mass matrix for classical beam elements


   REAL(ReKi), INTENT( IN)               :: A, L, Ixx, Iyy, Jzz, rho
   REAL(ReKi), INTENT( IN)               :: DirCos(3,3)
   
   REAL(ReKi)             :: M(12, 12)
         
   REAL(ReKi)                            :: t, rx, ry, po
   REAL(ReKi)                            :: DC(12, 12)
   
   t = rho*A*L;
   rx = rho*Ixx;
   ry = rho*Iyy;
   po = rho*Jzz*L;   

   M = 0
   
      
   M( 9,  9) = t/3.0
   M( 7,  7) = 13.0*t/35.0 + 6.0*ry/(5.0*L)
   M( 8,  8) = 13.0*t/35.0 + 6.0*rx/(5.0*L)
   M(12, 12) = po/3.0
   M(10, 10) = t*L*L/105.0 + 2.0*L*rx/15.0
   M(11, 11) = t*L*L/105.0 + 2.0*L*ry/15.0
   M( 2,  4) = -11.0*t*L/210.0 - rx/10.0
   M( 1,  5) =  11.0*t*L/210.0 + ry/10.0
   M( 3,  9) = t/6.0
   M( 5,  7) =  13.*t*L/420. - ry/10.
   M( 4,  8) = -13.*t*L/420. + rx/10. 
   M( 6, 12) = po/6.
   M( 2, 10) =  13.*t*L/420. - rx/10. 
   M( 1, 11) = -13.*t*L/420. + ry/10.
   M( 8, 10) =  11.*t*L/210. + rx/10.
   M( 7, 11) = -11.*t*L/210. - ry/10. 
   M( 1,  7) =  9.*t/70. - 6.*ry/(5.*L)
   M( 2,  8) =  9.*t/70. - 6.*rx/(5.*L)
   M( 4, 10) = -L*L*t/140. - rx*L/30. 
   M( 5, 11) = -L*L*t/140. - ry*L/30.
   
   M( 3,  3) = M( 9,  9)
   M( 1,  1) = M( 7,  7)
   M( 2,  2) = M( 8,  8)
   M( 6,  6) = M(12, 12)
   M( 4,  4) = M(10, 10)
   M( 5,  5) = M(11, 11)
   M( 4,  2) = M( 2,  4)
   M( 5,  1) = M( 1,  5)
   M( 9,  3) = M( 3,  9)
   M( 7,  5) = M( 5,  7)
   M( 8,  4) = M( 4,  8)
   M(12,  6) = M( 6, 12)
   M(10,  2) = M( 2, 10)
   M(11,  1) = M( 1, 11)
   M(10,  8) = M( 8, 10)
   M(11,  7) = M( 7, 11)
   M( 7,  1) = M( 1,  7)
   M( 8,  2) = M( 2,  8)
   M(10,  4) = M( 4, 10)
   M(11,  5) = M( 5, 11)
   
   DC = 0
   DC( 1: 3,  1: 3) = DirCos
   DC( 4: 6,  4: 6) = DirCos
   DC( 7: 9,  7: 9) = DirCos
   DC(10:12, 10:12) = DirCos
   
   M = MATMUL( MATMUL(DC, M), TRANSPOSE(DC) )

END SUBROUTINE ElemM


!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------

SUBROUTINE ApplyConstr(Init,p)


   TYPE(SD_InitType), INTENT(INOUT)   :: Init
   TYPE(SD_ParameterType), INTENT(IN)   :: p
   
   INTEGER                  :: I !, J, k
   INTEGER                  :: row_n !bgn_j, end_j, 
   
   DO I = 1, p%NReact*6
      row_n = Init%BCs(I, 1)
      
      IF (Init%BCs(I, 2) == 1) THEN
         Init%K(row_n, :) = 0
         Init%K(:, row_n) = 0
         Init%K(row_n, row_n) = 1
      
         Init%M(row_n, :) = 0
         Init%M(:, row_n) = 0
         Init%M(row_n, row_n) = 0 !0.00001          !what is this???? I changed this to 0.  RRD 7/31
      ENDIF
      
   ENDDO ! I

      

END SUBROUTINE ApplyConstr
!------------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------------
SUBROUTINE ElemG(A, L, rho, DirCos, F, g)
! this routine calculates the lumped forces and moments due to gravity on a given element:
! the element has two nodes, with the loads for both elements stored in array F. Indexing of F is:
!  Fx_n1=1,Fy_n1=2,Fz_n1=3,Mx_n1= 4,My_n1= 5,Mz_n1= 6,
!  Fx_n2=7,Fy_n2=8,Fz_n2=9,Mx_n2=10,My_n2=11,Mz_n2=12
!------------------------------------------------------------------------------------------------------
   REAL(ReKi), INTENT( OUT)           :: F(12)        ! returned loads. positions 1-6 are the loads for node 1; 7-12 are loads for node 2.
   REAL(ReKi), INTENT( IN )           :: A            ! area
   REAL(ReKi), INTENT( IN )           :: g            ! gravity
   REAL(ReKi), INTENT( IN )           :: L            ! element length
   REAL(ReKi), INTENT( IN )           :: rho           
   REAL(ReKi), INTENT( IN )           :: DirCos(3, 3) ! direction cosine matrix (for determining distance between nodes 1 and 2)

   REAL(ReKi)                         :: TempCoeff
   REAL(ReKi)                         :: w            ! weight per unit length

   
   F = 0             ! initialize whole array to zero, then set the non-zero portions
   w = rho*A*g       ! weight per unit length
   
      ! lumped forces on both nodes (z component only):
   F(3) = -0.5*L*w 
   F(9) = F(3)
          
      ! lumped moments on node 1 (x and y components only):
   ! bjj: note that RRD wants factor of 1/12 because of boundary conditions. Our MeshMapping routines use factor of 1/6 (assuming generic/different boundary  
   !      conditions), so we may have some inconsistent behavior. JMJ suggests using line2 elements for SubDyn's input/output meshes to improve the situation.
   TempCoeff = L*L*w/12.0_ReKi ! let's not calculate this twice  
   F(4) = -TempCoeff * DirCos(2,3) ! = -L*w*Dy/12.   !bjj: DirCos(2,3) = Dy/L
   F(5) =  TempCoeff * DirCos(1,3) ! =  L*w*Dx/12.   !bjj: DirCos(1,3) = Dx/L

      ! lumped moments on node 2: (note the opposite sign of node 1 moment)
   F(10) = -F(4)
   F(11) = -F(5)
   !F(12) is 0 for g along z alone
   
   
END SUBROUTINE ElemG
!------------------------------------------------------------------------------------------------------
SUBROUTINE LumpForces(Area1,Area2,crat,L,rho, g, DirCos, F)
!bjj: note this routine is a work in progress, intended for future version of SubDyn. Compare with ElemG.

         !This rountine calculates the lumped gravity forces at the nodes given the element geometry
         !It assumes a linear variation of the dimensions from node 1 to node 2, thus the area may be quadratically varying if crat<>1
   REAL(ReKi), INTENT( OUT)           :: F(12)
   REAL(ReKi), INTENT( IN )           :: Area1,Area2,crat !X-sectional areas at node 1 and node 2, t2/t1 thickness ratio
   REAL(ReKi), INTENT( IN )           :: g !gravity
   REAL(ReKi), INTENT( IN )           :: L !Length of element
   REAL(ReKi), INTENT( IN )           :: rho !density
   REAL(ReKi), INTENT( IN )           :: DirCos(3, 3)

    !LOCALS
   REAL(ReKi)                         :: TempCoeff,a0,a1,a2  !coefficients of the gravity quadratically distributed force

   A1 = -99999  !bjj initialized this to avoid getting warning by Intel Analyzers; needs to be set differently
   A2 = -99999  !bjj initialized this to avoid getting warning by Intel Analyzers; needs to be set differently
   
   !Calculate quadratic polynomial coefficients
   a0=A1
   a2=( (Area1+A2) - (Area1*crat+Area2/crat) )/L**2.  !*x**2
   a1= (Area2-Area1)/L -a2*L                       !*x                   
   
   !Now calculate the Lumped Forces
   F = 0
   F(3) = -(a0*L/2. +a1*L**2/6. +a2*L**3/12. )*rho*g  !Forces along z (must be negative on earth)
   F(9) = -(a0*L/2. +a1*L**2/3. +a2*L**3/4.  )*rho*g  !Forces along z (must be negative on earth)

   !Now calculate the Lumped Moments
   !HERE TO BE COMPLETED FOR THE BELOW
   TempCoeff = 1.0/12.0*g*L*L*rho*Area2  !RRD : I am changing this to >0 sign 6/10/13
      
   !F(4) = TempCoeff*( DirCos(1, 3)*DirCos(2, 1) - DirCos(1, 1)*DirCos(2, 3) ) !These do not work if convnetion on z2>z1, x2>x1, y2>y1 are not followed as I have discovered 7/23
   !F(5) = TempCoeff*( DirCos(1, 3)*DirCos(2, 2) - DirCos(1, 2)*DirCos(2, 3) ) 
   
   !RRD attempt at new dircos which keeps x in the X-Y plane
         F(4) = -TempCoeff * SQRT(1-DirCos(3,3)**2) * DirCos(1,1) !bjj: compare with ElemG() and verify this lumping is consistent
         F(5) = -TempCoeff * SQRT(1-DirCos(3,3)**2) * DirCos(2,1) !bjj: compare with ElemG() and verify this lumping is consistent
    !RRD ends
   F(10) = -F(4)
   F(11) = -F(5)
   !F(12) is 0 for g along z alone

   
END SUBROUTINE LumpForces

END MODULE SD_FEM
