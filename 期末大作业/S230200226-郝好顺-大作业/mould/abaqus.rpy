# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.14-4 replay file
# Internal Version: 2015_06_12-04.41.13 135079
# Run by haoha on Wed Jan 10 14:40:08 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=217.015625, 
    height=131.570602416992)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(14.0, 0.0))
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseShellExtrude(sketch=s, depth=100.0)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)
s1.rectangle(point1=(-30.0, 20.0), point2=(30.0, -20.0))
p = mdb.models['Model-1'].Part(name='Part-2', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-2']
p.BaseShell(sketch=s1)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Part-2']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models['Model-1'].sketches['__profile__']
p = mdb.models['Model-1'].parts['Part-2']
p.ReferencePoint(point=(0.0, 0.0, 0.0))
p1 = mdb.models['Model-1'].parts['Part-2']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
p = mdb.models['Model-1'].Part(name='Part-3', 
    objectToCopy=mdb.models['Model-1'].parts['Part-2'])
session.viewports['Viewport: 1'].setValues(displayedObject=p)
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
mdb.models['Model-1'].Material(name='Material-1')
mdb.models['Model-1'].materials['Material-1'].Elastic(table=((200000000000.0, 
    0.3), ))
mdb.models['Model-1'].HomogeneousShellSection(name='Section-1', 
    preIntegrate=OFF, material='Material-1', thicknessType=UNIFORM, 
    thickness=0.02, thicknessField='', idealization=NO_IDEALIZATION, 
    poissonDefinition=DEFAULT, thicknessModulus=None, temperature=GRADIENT, 
    useDensity=OFF, integrationRule=SIMPSON, numIntPts=3)
p = mdb.models['Model-1'].parts['Part-1']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
region = regionToolset.Region(faces=faces)
p = mdb.models['Model-1'].parts['Part-1']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
p = mdb.models['Model-1'].parts['Part-2']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['Part-2']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
region = regionToolset.Region(faces=faces)
p = mdb.models['Model-1'].parts['Part-2']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
p = mdb.models['Model-1'].parts['Part-3']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
p = mdb.models['Model-1'].parts['Part-3']
f = p.faces
faces = f.getSequenceFromMask(mask=('[#1 ]', ), )
region = regionToolset.Region(faces=faces)
p = mdb.models['Model-1'].parts['Part-3']
p.SectionAssignment(region=region, sectionName='Section-1', offset=0.0, 
    offsetType=MIDDLE_SURFACE, offsetField='', 
    thicknessAssignment=FROM_SECTION)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=OFF)
p = mdb.models['Model-1'].parts['Part-2']
a.Instance(name='Part-2-1', part=p, dependent=OFF)
p = mdb.models['Model-1'].parts['Part-3']
a.Instance(name='Part-3-1', part=p, dependent=OFF)
a = mdb.models['Model-1'].rootAssembly
a.translate(instanceList=('Part-3-1', ), vector=(0.0, 0.0, 100.0))
#: The instance Part-3-1 was translated by 0., 0., 100. (相对于装配坐标系)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
mdb.models['Model-1'].BuckleStep(name='Step-1', previous='Initial', 
    numEigen=10, vectors=18, maxIterations=50)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
mdb.models['Model-1'].fieldOutputRequests['F-Output-1'].setValues(variables=(
    'U', 'UT', 'UR', 'RBANG', 'RBROT'))
session.viewports['Viewport: 1'].assemblyDisplay.setValues(interactions=ON, 
    constraints=ON, connectors=ON, engineeringFeatures=ON, 
    adaptiveMeshConstraints=OFF)
mdb.models['Model-1'].ContactProperty('IntProp-1')
mdb.models['Model-1'].interactionProperties['IntProp-1'].TangentialBehavior(
    formulation=FRICTIONLESS)
mdb.models['Model-1'].interactionProperties['IntProp-1'].NormalBehavior(
    pressureOverclosure=HARD, allowSeparation=ON, 
    constraintEnforcementMethod=DEFAULT)
#: 相互作用属性 "IntProp-1" 已创建.
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Initial')
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Part-2-1'].faces
side2Faces1 = s1.getSequenceFromMask(mask=('[#1 ]', ), )
region1=regionToolset.Region(side2Faces=side2Faces1)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#2 ]', ), )
region2=regionToolset.Region(edges=edges1)
mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Int-1', 
    createStepName='Initial', master=region1, slave=region2, sliding=FINITE, 
    enforcement=NODE_TO_SURFACE, thickness=OFF, 
    interactionProperty='IntProp-1', surfaceSmoothing=NONE, adjustMethod=NONE, 
    smooth=0.2, initialClearance=OMIT, datumAxis=None, clearanceRegion=None)
#: 相互作用 "Int-1" 已创建.
a = mdb.models['Model-1'].rootAssembly
s1 = a.instances['Part-3-1'].faces
side1Faces1 = s1.getSequenceFromMask(mask=('[#1 ]', ), )
region1=regionToolset.Region(side1Faces=side1Faces1)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
edges1 = e1.getSequenceFromMask(mask=('[#1 ]', ), )
region2=regionToolset.Region(edges=edges1)
mdb.models['Model-1'].Tie(name='Constraint-1', master=region1, slave=region2, 
    positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)
session.viewports['Viewport: 1'].view.setValues(nearPlane=183.457, 
    farPlane=337.817, width=160.059, height=71.1316, cameraPosition=(71.6748, 
    188.617, -114.992), cameraUpVector=(-0.687567, -0.583842, -0.431718))
session.viewports['Viewport: 1'].view.setValues(nearPlane=175.3, 
    farPlane=345.973, width=259.01, height=115.106, viewOffsetX=9.78767, 
    viewOffsetY=-2.47856)
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-2-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
region2=regionToolset.Region(faces=faces1)
a = mdb.models['Model-1'].rootAssembly
r1 = a.instances['Part-2-1'].referencePoints
refPoints1=(r1[2], )
region1=regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].RigidBody(name='Constraint-2', refPointRegion=region1, 
    bodyRegion=region2)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, interactions=OFF, constraints=OFF, 
    engineeringFeatures=OFF)
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-2-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
region = regionToolset.Region(faces=faces1)
mdb.models['Model-1'].DisplacementBC(name='BC-1', createStepName='Initial', 
    region=region, u1=SET, u2=SET, u3=UNSET, ur1=SET, ur2=SET, ur3=SET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)
del mdb.models['Model-1'].boundaryConditions['BC-1']
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-3-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
region = regionToolset.Region(faces=faces1)
mdb.models['Model-1'].EncastreBC(name='BC-1', createStepName='Initial', 
    region=region, localCsys=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='Step-1')
a = mdb.models['Model-1'].rootAssembly
r1 = a.instances['Part-2-1'].referencePoints
refPoints1=(r1[2], )
region = regionToolset.Region(referencePoints=refPoints1)
mdb.models['Model-1'].ConcentratedForce(name='Load-1', createStepName='Step-1', 
    region=region, cf3=500.0, distributionType=UNIFORM, field='', 
    localCsys=None)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON, loads=OFF, 
    bcs=OFF, predefinedFields=OFF, connectors=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
e2 = a.instances['Part-2-1'].edges
e3 = a.instances['Part-3-1'].edges
pickedEdges = e1.getSequenceFromMask(mask=('[#3 ]', ), )+\
    e2.getSequenceFromMask(mask=('[#f ]', ), )+e3.getSequenceFromMask(mask=(
    '[#f ]', ), )
a.seedEdgeByNumber(edges=pickedEdges, number=25, constraint=FINER)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
e2 = a.instances['Part-2-1'].edges
e3 = a.instances['Part-3-1'].edges
pickedEdges = e1.getSequenceFromMask(mask=('[#3 ]', ), )+\
    e2.getSequenceFromMask(mask=('[#f ]', ), )+e3.getSequenceFromMask(mask=(
    '[#f ]', ), )
a.seedEdgeByNumber(edges=pickedEdges, number=5, constraint=FINER)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-2-1'].edges
e2 = a.instances['Part-3-1'].edges
pickedEdges = e1.getSequenceFromMask(mask=('[#f ]', ), )+\
    e2.getSequenceFromMask(mask=('[#f ]', ), )
a.seedEdgeByNumber(edges=pickedEdges, number=10, constraint=FINER)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
pickedEdges = e1.getSequenceFromMask(mask=('[#3 ]', ), )
a.seedEdgeByNumber(edges=pickedEdges, number=40, constraint=FINER)
elemType1 = mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD, 
    secondOrderAccuracy=OFF, hourglassControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
a = mdb.models['Model-1'].rootAssembly
f1 = a.instances['Part-1-1'].faces
faces1 = f1.getSequenceFromMask(mask=('[#1 ]', ), )
f2 = a.instances['Part-2-1'].faces
faces2 = f2.getSequenceFromMask(mask=('[#1 ]', ), )
f3 = a.instances['Part-3-1'].faces
faces3 = f3.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =((faces1+faces2+faces3), )
a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2))
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['Part-1-1'], a.instances['Part-2-1'], 
    a.instances['Part-3-1'], )
a.generateMesh(regions=partInstances)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=OFF)
mdb.Job(name='Job-1', model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
    numGPUs=0)
mdb.jobs['Job-1'].submit(consistencyChecking=OFF)
#: 作业输入文件 "Job-1.inp" 已提交分析.
#: Job Job-1: Analysis Input File Processor completed successfully.
#: Job Job-1: Abaqus/Standard completed successfully.
#: Job Job-1 completed successfully. 
o3 = session.openOdb(name='C:/Temp/Job-1.odb')
#: 模型: C:/Temp/Job-1.odb
#: 装配件个数:         1
#: 装配件实例个数: 0
#: 部件实例的个数:     3
#: 网格数:             3
#: 单元集合数:       1
#: 结点集合数:          2
#: 分析步的个数:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=2 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=3 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=4 )
session.viewports['Viewport: 1'].odbDisplay.setFrame(step=0, frame=5 )
mdb.saveAs(pathName='C:/Temp/finalwork')
#: 模型数据库已保存到 "C:\Temp\finalwork.cae".
