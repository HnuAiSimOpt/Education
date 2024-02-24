#%%
import numpy as np
import matplotlib.pyplot as plt
import time
import vtkmodules.all as vtk

#%%
#生成可视化vtk文件
class FEDataModel:
    """有限元数据模型类"""

    def __init__(self):
        self.nodes = []  # 节点几何坐标
        self.elements = []  # 单元拓扑信息
        self.s = []
        self.scalars = {}  # 节点标量属性
        self.vectors = {}  # 节点向量属性
        self.ugrid = vtk.vtkUnstructuredGrid()  # 用于VTK可视化的数据模型
        self.ugrid.Allocate(100)

    # 得到节点坐标单元节点编号
    def read_nodes_elements(self, node_file, element_file, s_file):
        with open(node_file) as f:
            for line in f.readlines():
                line = line.strip("\n")
                self.nodes.append(list(map(lambda x: float(x), line.split(",")))[1:])
            f.close()
        with open(element_file) as f:
            for line in f.readlines():
                line = line.strip("\n")
                self.elements.append(list(map(lambda x: int(x), line.split(",")))[1:])
            f.close()
        with open(s_file) as f:
            for line in f.readlines():
                line = line.strip("\n")
                self.s.append(list(map(lambda x: float(x), line.split(","))))
            f.close()

        nodes = vtk.vtkPoints()
        for i in range(0, len(self.nodes)):
            nodes.InsertPoint(i, self.nodes[i])

        for i in range(0, len(self.elements)):
            try:
                hexahedron = vtk.vtkHexahedron()
                for j in range(8):
                    hexahedron.GetPointIds().SetId(j, self.elements[i][j] - 1)
                self.ugrid.InsertNextCell(hexahedron.GetCellType(), hexahedron.GetPointIds())
            except Exception as err:
                print("FEDataModel构建中遇到错误单元类型！")
                print(err)
        self.ugrid.SetPoints(nodes)

    # 获得标量信息，应力、温度场等等
    def read_ntl(self):

        scalar = self.s
        # 存储标量值
        scalars = vtk.vtkFloatArray()
        scalars.SetName("S")
        for i in range(0, len(scalar)):
            scalars.InsertTuple1(i, scalar[i][0])
        # 设定每个节点的标量值
        self.ugrid.GetPointData().SetScalars(scalars)

    def display(self):
        renderer = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(renderer)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)

        colors = vtk.vtkNamedColors()
        ugridMapper = vtk.vtkDataSetMapper()
        ugridMapper.SetInputData(self.ugrid)

        ugridActor = vtk.vtkActor()
        ugridActor.SetMapper(ugridMapper)
        ugridActor.GetProperty().SetColor(colors.GetColor3d("AliceBlue"))
        ugridActor.GetProperty().EdgeVisibilityOn()

        renderer.AddActor(ugridActor)
        renderer.SetBackground(colors.GetColor3d("AliceBlue"))

        renderer.ResetCamera()
        renderer.GetActiveCamera().Elevation(60.0)
        renderer.GetActiveCamera().Azimuth(30.0)
        renderer.GetActiveCamera().Dolly(1.2)
        renWin.SetSize(640, 480)
        # Interact with the data.
        renWin.Render()
        iren.Start()

    def drawScalarField(self, scalar_mapper, scalarRange, title):
        # 定义颜色映射表
        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.5, 0.0)  # 色调范围从红色到蓝色
        lut.SetAlphaRange(1.0, 1.0)  # 透明度范围
        lut.SetValueRange(1.0, 1.0)
        lut.SetSaturationRange(0.5, 0.5)  # 颜色饱和度
        lut.SetNumberOfTableValues(16)
        lut.SetNumberOfColors(16)  # 颜色个数
        lut.SetRange(scalarRange)
        lut.Build()

        scalar_mapper.SetScalarRange(scalarRange)
        scalar_mapper.SetLookupTable(lut)
        scalar_actor = vtk.vtkActor()
        scalar_actor.SetMapper(scalar_mapper)
        self.renderer.AddActor(scalar_actor)
        # 色标带
        scalarBar = vtk.vtkScalarBarActor()
        scalarBar.SetLookupTable(scalar_mapper.GetLookupTable())  # 将颜色查找表传入窗口中的色标带
        scalarBar.SetTitle(title)
        scalarBar.SetNumberOfLabels(5)
        self.renderer.AddActor2D(scalarBar)

    def save_vtk(self, filename):
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(filename)
        writer.SetInputData(self.ugrid)
        writer.Write()

for i in ["1", "2", "3"]:
    model = FEDataModel()
    model.read_nodes_elements("Nodes.txt", "3D4C_Elements.txt", "U"+i+".txt")
    model.read_ntl()
    model.display()
    model.save_vtk("visualize/" + inp_file + "U"+i+".vtk")
print("Done")