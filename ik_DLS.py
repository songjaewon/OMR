## Programmer: So-Yeong Jeon
## Author(Programmer): So-Yeong Jeon
## Creation Date: 12/07/2011
##
## Description:
##
## This is an implementation of Damped Least Square,
## The method is first used for inverse kinematics by
## Wampler[41] and Nakamura and Hanafusa[34].
## (See the unpublished survey paper,
## 'Introduction to Inverse Kinematics with
## Jacobian Transpose, Pseudoinverse and Damped Least Squares methods'
## written by Samuel R. Buss.)
##
## See "IK_Method_HowToUse.docx" which is accompanied by this code
## 
##
####################################################
##//Double Y-Shape
##joint -p -12 0 0 -n joint0;
##    joint -p -8 0 0 -n joint1;
##      select joint1;        
##      joint -p -4 0 0 -n joint2;
##        joint -p 0 0 0 -n joint3;
##          select joint3;
##          joint -p 4 0 0 -n joint4; 
##            joint -p 4 0 4 -n joint5;
##              joint -p 8 0 4 -n joint6;
##                joint -p 12 0 4 -n joint7;
##                joint -e -sao yup -zso -oj xyz joint7;
##              joint -e -sao yup -zso -oj xyz joint6;
##            joint -e -sao yup -zso -oj xyz joint5;      
##          joint -e -sao yup -zso -oj xyz joint4;
##          select joint3;           
##          joint -p 4 0 0 -n joint8; 
##            joint -p 4 0 -4 -n joint9;
##              joint -p 8 0 -4 -n joint10;
##               joint -p 12 0 -4 -n joint11;
##               joint -e -sao yup -zso -oj xyz joint11;
##              joint -e -sao yup -zso -oj xyz joint10;
##            joint -e -sao yup -zso -oj xyz joint9;   
##          joint -e -sao yup -zso -oj xyz joint8;    
##        joint -e -sao yup -zso -oj xyz joint3; 
##      joint -e -sao yup -zso -oj xyz joint2; 
##      select joint1;        
##      joint -p -4 0 0 -n joint12;
##        joint -p -4 0 -4 -n joint13;         
##          select joint13;
##          joint -p -4 0 -8 -n joint14;
##            joint -p 0 0 -8 -n joint15; 
##              joint -p 0 0 -12 -n joint16; 
##                joint -p 0 0 -16 -n joint17;          
##                joint -e -sao yup -zso -oj xyz joint17; 
##              joint -e -sao yup -zso -oj xyz joint16;
##            joint -e -sao yup -zso -oj xyz joint15;
##          joint -e -sao yup -zso -oj xyz joint14;                    
##          select joint13;           
##          joint -p -4 0 -8 -n joint18;
##            joint -p -8 0 -8 -n joint19;
##              joint -p -8 0 -12 -n joint20;
##               joint -p -8 0 -16 -n joint21;
##                joint -e -sao yup -zso -oj xyz joint21; 
##              joint -e -sao yup -zso -oj xyz joint20;
##            joint -e -sao yup -zso -oj xyz joint19;
##          joint -e -sao yup -zso -oj xyz joint18;     
##       joint -e -sao yup -zso -oj xyz joint13;
##     joint -e -sao yup -zso -oj xyz joint12; 
##    joint -e -sao yup -zso -oj xyz joint1;   
##joint -e -sao yup -zso -oj xyz joint0;                       
##CreateLocator;
##setAttr "locator1.translateX" 12;
##setAttr "locator1.translateY" 0;
##setAttr "locator1.translateZ" 4;
##CreateLocator;
##setAttr "locator2.translateX" 12;
##setAttr "locator2.translateY" 0;
##setAttr "locator2.translateZ" -4; 
##CreateLocator;
##setAttr "locator3.translateX" 0;
##setAttr "locator3.translateY" 0;
##setAttr "locator3.translateZ" -16;
##CreateLocator;
##setAttr "locator4.translateX" -8;
##setAttr "locator4.translateY" 0;
##setAttr "locator4.translateZ" -16; 


import maya.cmds as cmds
import maya.OpenMaya as OpenMaya
import numpy as np
import math as math
##eulerRotAngle=[90, -90, 0]
##target3Axis=[ [1,0,0], [0,1,0], [0,0,1] ]
##result3Axis=target3Axis

def parentJoints(endEffectorName):
 old=endEffectorName
 parents=[]
 curParent=cmds.pickWalk(endEffectorName,d='up')
 assert( len(curParent) == 1)
 while( old!=curParent[0] ):
   parents.append(curParent[0])
   old=curParent[0]
   curParent=cmds.pickWalk(curParent[0],d='up')
   assert( len(curParent) == 1)
 cmds.select(cl=True)
 return parents



def specifyJoints_EffectiveJoints(targets,endEffectorNames):
  assert( len(targets)==len(endEffectorNames) )
  joints=[]
  effectiveJointNames=[]
  for i in range(len(endEffectorNames)):
    effectiveJointNames.append(parentJoints(endEffectorNames[i]))
    joints.extend(effectiveJointNames[i])
  joints=list(set(joints))
  return joints, effectiveJointNames



def transform3Axis(eulerRotAngle,target3Axis, result3Axis, isDegree=True):
  """transform the target3Axis according to eulerRotAngle in XYZ order"""
  targetX=OpenMaya.MVector(target3Axis[0][0],target3Axis[0][1],target3Axis[0][2])
  targetY=OpenMaya.MVector(target3Axis[1][0],target3Axis[1][1],target3Axis[1][2])
  targetZ=OpenMaya.MVector(target3Axis[2][0],target3Axis[2][1],target3Axis[2][2])  
  if(isDegree==True):
    PI=3.1415926535897932384626433832795
    for i in range(3):
      eulerRotAngle[i]=eulerRotAngle[i]*(PI/180.0)
  rot=OpenMaya.MEulerRotation(eulerRotAngle[0],eulerRotAngle[1],eulerRotAngle[2],OpenMaya.MEulerRotation.kXYZ)
  targetX=targetX.rotateBy(rot)
  targetY=targetY.rotateBy(rot)
  targetZ=targetZ.rotateBy(rot)
  for i in range(3):
    result3Axis[0][i]=targetX[i]
  for i in range(3):
    result3Axis[1][i]=targetY[i]
  for i in range(3):
    result3Axis[2][i]=targetZ[i]



def ClampMag(w,Dmax,D=0,i=0):
  wNorm=math.sqrt(sum(w*w))
  if(wNorm>Dmax):
    w=w*(Dmax/wNorm)
  return w


def computeDi(err,oldErr):
  assert(len(err)==len(oldErr))
  d=np.zeros(len(err))
  for i in range(len(err)):
    new=math.sqrt(sum(err[i]*err[i]))
    old=math.sqrt(sum(oldErr[i]*oldErr[i]))
    if(old-new>0):
     d[i]=old-new
    else:
     d[i]=0
  return d


def ClampMagDiscrete(w,Dmax,D,i):
  if(len(D)!=0):
    Dmaxi=D[i]+Dmax
    w=ClampMag(w,Dmaxi)
  return w  


def DLS(jointPositions, rotateAxes, targetPositions,sourcePositions,effectiveJoints1D,D,Clamping,Dmax=2.0,dampingConstant=3.0):
  """DiscreteJumpDLS"""
  #jointPositions have each joint for each rotateAxis in rotateAxes
  # jointPositions may have duplicate joints for x/y/z axes
  #effectiveJoints[i] have indices of effective joints for sourcePositions[i]
  # the indices can be used on jointPositions
  jointAngles=np.zeros( (len(jointPositions)) )
  J=np.zeros( (len(sourcePositions)*3,len(jointPositions)) )
  e=np.zeros( (len(sourcePositions)*3) )
  for srcIdx in range(len(sourcePositions)):
    posDisplace=targetPositions[srcIdx]-sourcePositions[srcIdx]
    posDisplace=Clamping(posDisplace,Dmax,D,srcIdx)
    e[3*srcIdx+0]=posDisplace[0]
    e[3*srcIdx+1]=posDisplace[1]
    e[3*srcIdx+2]=posDisplace[2]
  #print e
  for srcIdx in range(len(sourcePositions)):
    for j in range(len(jointPositions)):
      if(j in effectiveJoints1D[srcIdx]):
        dSdTheta=np.cross(rotateAxes[j],sourcePositions[srcIdx]-jointPositions[j])
        J[3*srcIdx+0][j]=dSdTheta[0]
        J[3*srcIdx+1][j]=dSdTheta[1]
        J[3*srcIdx+2][j]=dSdTheta[2]
      else:
        J[3*srcIdx+0][j]=0
        J[3*srcIdx+1][j]=0
        J[3*srcIdx+2][j]=0
  #print J
  part=np.dot(J,J.T)+(dampingConstant*dampingConstant)*np.identity(len(sourcePositions)*3,float)
  J_inverse=np.dot(J.T,np.linalg.inv(part))
  jointAngles=np.dot(J_inverse,e)
  return jointAngles


def IKSolveByGivenNames(IKSolver,Clamping,targets,endEffectorNames,joints,effectiveJointNames,prevErr,Dmax=2.0,dampingConstant=3.0):
  targetPositions=np.zeros( (len(targets),3) )
  for i in range(len(targets)):
    targetPositions[i]=np.array(cmds.xform(targets[i],q=True,ws=True,rp=True))
  numJoints=len(joints)
  endEffectorPositions=np.zeros( (len(targets),3) )
  for i in range(len(targets)):
    endEffectorPositions[i]=np.array(cmds.xform(endEffectorNames[i],q=True,ws=True,rp=True))
  jointPositions=np.zeros( (numJoints,3) )
  jointAngles=np.zeros( (numJoints,3) )
  for i in range(numJoints):
    jointAngles[i]=np.array(cmds.xform(joints[i],q=True,ws=True,ro=True))
    jointPositions[i]=np.array(cmds.xform(joints[i],q=True,ws=True,rp=True))
  jointLocalAxes=np.zeros( (numJoints,3,3) )
  target3Axis=[ [1,0,0], [0,1,0], [0,0,1] ]
  for i in range(numJoints):
    temp_axis=[ [1,0,0], [0,1,0], [0,0,1] ]
    transform3Axis(jointAngles[i],target3Axis,temp_axis,isDegree=True)
    jointLocalAxes[i]=np.array(temp_axis)
    #print jointLocalAxes[i]
  jointPositions1D=np.zeros( (numJoints*3,3) )
  rotateAxes=np.zeros( (numJoints*3,3) )
  overallIdx=0
  for jointIdx in range(numJoints):
    for axisIdx in range(3):
      rotateAxes[overallIdx]=np.array(jointLocalAxes[jointIdx][axisIdx])
      jointPositions1D[overallIdx]=jointPositions[jointIdx]
      overallIdx=overallIdx+1 #!! not ++overallIdx  
  effectiveJoints1D=[]
  #effectiveJoints[i] have indices of effective joints to endEffectorNames[i]
  for i in range(len(endEffectorNames)):
    effectiveJoints1D.append(set())
    for jointName in effectiveJointNames[i]:
      #print endEffectorNames[i], jointName
      print "joints.index(jointName) : ", joints.index(jointName)
      effectiveJoints1D[i].add(3*joints.index(jointName))
      effectiveJoints1D[i].add(3*joints.index(jointName)+1)
      effectiveJoints1D[i].add(3*joints.index(jointName)+2)
  #print joints
  print effectiveJoints1D
  jointAngles1D=np.zeros( (numJoints*3) )
  if( len(prevErr)==0):
    D=np.zeros(0)
  else:
    D=computeDi(targetPositions-endEffectorPositions,prevErr)
  jointAngles1D=IKSolver(jointPositions1D,rotateAxes,targetPositions,endEffectorPositions,effectiveJoints1D,D,Clamping,Dmax,dampingConstant)
  PI=3.1415926535897932384626433832795
  for jointIdx in range(numJoints):
    degX=jointAngles1D[3*jointIdx]*(180.0/PI)
    degY=jointAngles1D[3*jointIdx+1]*(180.0/PI)
    degZ=jointAngles1D[3*jointIdx+2]*(180.0/PI)
    cmds.xform(joints[jointIdx],eu=True,r=True,ro=[degX, degY, degZ])
  prevErr=targetPositions-endEffectorPositions 
  return prevErr

#Usage
targetsA=['locator1', 'locator2']
endEffectorNamesA=['joint4', 'joint7']
errA=np.zeros(0)
jointsA, effectiveJointNamesA=specifyJoints_EffectiveJoints(targetsA,endEffectorNamesA)
print "jointsA : ", jointsA
print "effectiveJoints : ", effectiveJointNamesA
errA=IKSolveByGivenNames(DLS,ClampMag,targetsA,endEffectorNamesA,jointsA,effectiveJointNamesA,errA,Dmax=2.0,dampingConstant=3)

