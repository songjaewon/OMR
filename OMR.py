import maya.cmds as mc
import maya.OpenMaya as OpenMaya
import numpy as np
import math as math

class OMR() :
    def __init__(self, targetList, endEffectorList):
        self.targetList = targetList
        self.endEffectorList = endEffectorList   
        self.jntList, self.effectiveJntDic = self.getJntListAndEffectiveJntDic(self.targetList, self.endEffectorList)
        self.PI = 3.1415926535897932384626433832795     
        return    
    
    def getLocalAxes(self, eulerRotAngle, isEulerAngle=True):        
        axisX = OpenMaya.MVector(1,0,0)
        axisY = OpenMaya.MVector(0,1,0)
        axisZ = OpenMaya.MVector(0,0,1)        
        if(isEulerAngle==True):            
            for i in range(3):
                eulerRotAngle[i]=eulerRotAngle[i]*(self.PI/180.0)
        rot = OpenMaya.MEulerRotation(eulerRotAngle[0],eulerRotAngle[1],eulerRotAngle[2],OpenMaya.MEulerRotation.kXYZ)
        axisX = axisX.rotateBy(rot)
        axisY = axisY.rotateBy(rot)
        axisZ = axisZ.rotateBy(rot)
        tempAxis = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
        for i in range(3):
            tempAxis[0][i] = axisX[i]
            tempAxis[1][i] = axisY[i]
            tempAxis[2][i] = axisZ[i]
        return tempAxis
            
    def parentJnts(self, endEffector):
        old = endEffector
        parentJntList = []
        currentParent = mc.pickWalk(endEffector, d='up')
        assert(len(currentParent)==1)
        while(old != currentParent[0]):
            parentJntList.append(currentParent[0])
            old=currentParent[0]
            currentParent = mc.pickWalk(currentParent[0], d='up')
            assert(len(currentParent)==1)
        mc.select(cl=True)
        return parentJntList
    
    def getJntListAndEffectiveJntDic(self, targetList, endEffectorList):
        assert(len(targetList)==len(endEffectorList))
        jntList=[]
        effectiveJntDic={}
        for i in range(len(endEffectorList)):
            effectiveJntDic[endEffectorList[i]] = (self.parentJnts(endEffectorList[i]))
            jntList.extend(self.parentJnts(endEffectorList[i]))
        jntList = list(set(jntList))
        return jntList, effectiveJntDic   
    
    def clampMagDiscrete(self):
        return
    
    def clampMag(self, posDisplace, Dmax):
        if (sum(posDisplace*posDisplace)**0.5 <= Dmax ) :
            result = posDisplace
        else :
            result = Dmax*(posDisplace/(sum(posDisplace*posDisplace)**0.5))            
        return result
    
    def DampedLeastSquare(self, jntPosList, jntLocalAxesMat, targetPosList, endEffectorPosList, D=0.0, Dmax=2.0, dampingConstant=3.0):
        J = np.zeros((len(self.endEffectorList)*3, len(self.jntList)*3))
        e = np.zeros((len(self.endEffectorList)*3))
        for i in range(len(self.endEffectorList)):
            posDisplace = targetPosList[i] - endEffectorPosList[i]
            posDisplace = self.clampMag(posDisplace, Dmax)
            e[3*i + 0] = posDisplace[0]
            e[3*i + 1] = posDisplace[1]
            e[3*i + 2] = posDisplace[2]        
        for i in range(len(self.endEffectorList)):
            for j in range(len(self.jntList)):
                if self.jntList[j] in self.effectiveJntDic[self.endEffectorList[i]]:
                    for k in range(3) :
                        element = np.cross(jntLocalAxesMat[j][k], endEffectorPosList[i] - jntPosList[j])
                        J[3*i+0][3*j+k] = element[0]
                        J[3*i+1][3*j+k] = element[1]
                        J[3*i+2][3*j+k] = element[2]
                else :
                    for k in range(3):
                        J[3*i+0][3*j+k] = 0
                        J[3*i+1][3*j+k] = 0
                        J[3*i+2][3*j+k] = 0        
        temp = np.dot(J,J.T)+(dampingConstant*dampingConstant)*np.identity(len(self.endEffectorList)*3,float)
        J_inverse = np.dot(J.T, np.linalg.inv(temp))
        jointAngles = np.dot(J_inverse, e)               
        return jointAngles
    
    def ikSolver(self):
        targetPosMat = np.zeros((len(self.targetList), 3))
        endEffectorPosMat = np.zeros((len(self.targetList), 3))
        jntPosMat = np.zeros((len(self.jntList), 3))
        jntAngleMat = np.zeros((len(self.jntList), 3))
        for i in range(len(self.targetList)):
            targetPosMat[i] = np.array(mc.xform(self.targetList[i], q=True, ws=True, t=True))
            endEffectorPosMat[i] = np.array(mc.xform(self.endEffectorList[i], q=True, ws=True, t=True))
        for i in range(len(self.jntList)):
            jntPosMat[i] = np.array(mc.xform(self.jntList[i], q=True, ws=True, t=True))
            jntAngleMat[i] = np.array(mc.xform(self.jntList[i], q=True, ws=True, ro=True))
        jntLocalAxesMat = np.zeros((len(self.jntList),3, 3))       
        for i in range(len(self.jntList)):
            jntLocalAxesMat[i] = self.getLocalAxes(jntAngleMat[i])
        jntAngles = self.DampedLeastSquare(jntPosMat, jntLocalAxesMat, targetPosMat, endEffectorPosMat) 
        for i in range(len(self.jntList)):
            degX = jntAngles[3*i]*(180/self.PI)
            degY = jntAngles[3*i+1]*(180/self.PI)
            degZ = jntAngles[3*i+2]*(180/self.PI)
            mc.xform(self.jntList[i], eu=True, r=True, ro=[degX, degY, degZ])


#a = OMR(['locator1', 'locator2', 'locator3'],['joint4', 'joint7', 'joint10'])
a = OMR(['locator1','locator2','locator3','locator4'],['joint5','joint9','joint13','joint17'])
#a.ikSolver(['pSphere1', 'pSphere2'], [])
a.ikSolver()

