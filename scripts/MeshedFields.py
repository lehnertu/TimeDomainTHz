"""
This module allows the creation and handling of time-domain electromagnetic fields
defined on  surfaces given by unstructured triangular meshes.
"""

import vtk
import pygmsh
import numpy as np
from scipy import constants
import h5py

def CircularMesh(R, ratio=1.0, lcar=0.1):
    """
    Create a circular mesh with radius R around the origin [0,0,0].
    The normal vector is (0,0,1) in z-direction, the disk lies in the x-y-plane.
    
    The density of the grid is increased by ratio in the x-direction.
    The average cell dimensions are given by lcar. The average cell dimension in
    x-direction, thus, is lcar/ratio.
    
    returns:
    - an array of point corrdinates [3 float]
    - a connection list of triangles [3 int references into the array of points]
    """
    geom = pygmsh.built_in.Geometry()
    # we create the initial geometry as a streched ellipse to create
    # different scaling lengths (cell sizes) along the different axes
    p1 = geom.add_point([ratio, 0.0, 0.0], lcar)
    p2 = geom.add_point([0.0, 1.0, 0.0], lcar)
    p3 = geom.add_point([-ratio, 0.0, 0.0], lcar)
    p4 = geom.add_point([0.0, -1.0, 0.0], lcar)
    pc = geom.add_point([0.0, 0.0, 0.0], lcar)
    pa = geom.add_point([1.0, 0.0, 0.0], lcar)
    # the mesh is circumscribed with four elliptic arcs
    e1 = geom.add_ellipse_arc(start=p1, center=pc, point_on_major_axis=pa, end=p2)
    e2 = geom.add_ellipse_arc(start=p2, center=pc, point_on_major_axis=pa, end=p3)
    e3 = geom.add_ellipse_arc(start=p3, center=pc, point_on_major_axis=pa, end=p4)
    e4 = geom.add_ellipse_arc(start=p4, center=pc, point_on_major_axis=pa, end=p1)
    # these are combined into a line loop
    ll = geom.add_line_loop([e1,e2,e3,e4])
    geom.add_plane_surface(ll)
    # now we can create the mesh
    mesh = pygmsh.generate_mesh(geom, dim=2, verbose=False)
    # we reverse the streching by scaling the coordinates accordingly
    points = np.array([ p*R*[1.0/ratio,1,1] for p in mesh.points ])
    triangles = mesh.cells['triangle']
    return points, triangles

def MeshArea(points, triangles):
    """
    Compute the area of all mesh cells
    """
    area=[]
    for i, t in enumerate(triangles):
        p1 = points[t[0]]
        p2 = points[t[1]]
        p3 = points[t[2]]
        r1 = p2-p1
        r2 = p3-p1
        area.append(0.5*np.linalg.norm(np.cross(r1,r2)))
    return np.array(area)

def MeshNormals(points, triangles):
    """
    Compute the normal vectors of all mesh cells
    """
    normals=[]
    for i, t in enumerate(triangles):
        p1 = points[t[0]]
        p2 = points[t[1]]
        p3 = points[t[2]]
        r1 = p2-p1
        r2 = p3-p1
        n = np.cross(r1,r2)
        normals.append(n / np.linalg.norm(n))
    return np.array(normals)

def EnergyFlowDensity(A,dt):
    """
    Compute the Poynting vector (energy flow density) of the field
    integrated over all time with time step dt
    """
    Energy = []
    for trace in A:
        EVec = trace[:,0:3]
        BVec = trace[:,3:6]
        SVec = np.cross(EVec, BVec) / constants.mu_0
        Energy.append((SVec.sum(axis=0))*dt)
    return np.array(Energy)

def powerLUT():
    lut = vtk.vtkLookupTable()
    nc = 256
    ctf = vtk.vtkColorTransferFunction()
    # black-red-yellow-white
    # ctf.SetColorSpaceToDiverging()
    ctf.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
    ctf.AddRGBPoint(0.5, 0.7, 0.0, 0.0)
    ctf.AddRGBPoint(0.8, 0.7, 0.7, 0.0)
    ctf.AddRGBPoint(1.0, 1.0, 1.0, 1.0)
    lut.SetNumberOfTableValues(nc)
    lut.Build()
    for i in range(0, nc):
        rgb = list(ctf.GetColor(float(i) / nc))
        rgb.append(1.0)
        lut.SetTableValue(i, *rgb)
    return lut

class MyInteractor(vtk.vtkInteractorStyleTrackballCamera):

    def __init__(self, textMapper, parent=None):
        self.AddObserver("LeftButtonPressEvent", self.leftButtonPressEvent)
        selectedMapper = vtk.vtkDataSetMapper()
        selectedActor = vtk.vtkActor()
        self.text = textMapper

    def leftButtonPressEvent(self, obj, event):
        clickPos = self.GetInteractor().GetEventPosition()
        picker = vtk.vtkCellPicker()
        picker.Pick(clickPos[0], clickPos[1], 0, self.GetDefaultRenderer())
        id = picker.GetCellId()
        self.text.SetInput("Cell index:%d\npos:\nval=" % id)
        self.OnLeftButtonDown()
        return

def ShowMeshedField(points, triangles, centers=[], scalars=[], scalarTitle="", showAxes=False):
    """
    Render a display of a mesh geometry using VTK.
    If given the center points of the triangles are rendered.
    If given a scalar field is used to color the triangles accordingly.
    """
    # create a dataset for the triangle mesh
    pts = vtk.vtkPoints()
    for p in points:
        pts.InsertNextPoint(p)
    cells = vtk.vtkCellArray()
    Ncells = len(triangles)
    for t in triangles:
        cells.InsertNextCell(3, t)
    meshData = vtk.vtkPolyData()
    meshData.SetPoints(pts)
    meshData.SetPolys(cells)
    # map the triangle mesh into the scene
    meshMapper = vtk.vtkPolyDataMapper()
    meshMapper.SetInputData(meshData)
    # color triangles by given scalar value
    if len(scalars)>0:
        if len(scalars)==Ncells:
            scal = vtk.vtkFloatArray()
            scal.SetNumberOfValues(Ncells)
            for i,val in enumerate(scalars):
                scal.SetValue(i,val)
            lut = powerLUT()
            meshData.GetCellData().SetScalars(scal)
            meshMapper.SetLookupTable(lut)
            meshMapper.SetScalarRange(0.0,np.max(scalars))
            # create a color scale bar
            sbar = vtk.vtkScalarBarActor()
            sbar.SetLookupTable(lut)
            sbar.SetTitle(scalarTitle)
            sbar.SetPosition(0.05,0.4)
            sbar.SetPosition2(0.1,0.5)
        else:
            print("number of scalars differs from number of triangles - no coloring")
    # create a dataset for the center points
    if len(centers)>0:
        cpts = vtk.vtkPoints()
        for p in centers:
            cpts.InsertNextPoint(p)
        ccells = vtk.vtkCellArray()
        for i in range(len(centers)):
            ccells.InsertNextCell(1, [i])
        centerData = vtk.vtkPolyData()
        centerData.SetPoints(cpts)
        centerData.SetVerts(ccells)
    colors = vtk.vtkNamedColors()
    # add the triangle mesh actor to the scene
    meshActor = vtk.vtkActor()
    meshActor.SetMapper(meshMapper)
    meshActor.GetProperty().SetPointSize(5)
    meshActor.GetProperty().SetColor(colors.GetColor3d("Red"))
    meshActor.GetProperty().EdgeVisibilityOn()
    # map the center points into the scene
    if len(centers)>0:
        centerMapper = vtk.vtkPolyDataMapper()
        centerMapper.SetInputData(centerData)
        centerActor = vtk.vtkActor()
        centerActor.SetMapper(centerMapper)
        centerActor.GetProperty().SetPointSize(3)
        centerActor.GetProperty().SetColor(colors.GetColor3d("Blue"))
    # add some text to annotate the selected cell
    textMapper = vtk.vtkTextMapper()
    textMapper.SetInput("Cell index:\npos:\nval=")
    tprop = textMapper.GetTextProperty()
    tprop.SetJustificationToLeft()
    tprop.SetColor(colors.GetColor3d("LightBlue"))
    tprop.SetFontSize(20)
    textActor = vtk.vtkActor2D()
    textActor.SetMapper(textMapper)
    textActor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
    textActor.GetPositionCoordinate().SetValue(0.05, 0.2)
    # create a render window
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(colors.GetColor3d("SlateGray"))
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetSize(800,600)
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()
    # style = vtk.vtkInteractorStyleTrackballCamera()
    style = MyInteractor(textMapper)
    style.SetDefaultRenderer(renderer)
    renderWindowInteractor.SetInteractorStyle(style)
    # add the actors to the scene
    renderer.AddActor(meshActor)
    if len(centers)>0: renderer.AddActor(centerActor)
    # show the scalar bar
    if len(scalars)==Ncells:
        renderer.AddActor(sbar)
    # visualize the coordinate system
    if showAxes:
        axesActor = vtk.vtkAxesActor()
        renderer.AddActor(axesActor)
    renderer.AddActor(textActor)
    # render and interact
    renderWindow.Render()
    renderWindowInteractor.Start()
    # now the interaction is running until we close the window
    # cleanup after closing the window
    del renderWindow
    del renderWindowInteractor

def WriteMesh(filename, points, triangles, pos):
    """
    Write the meshed geometry to an HDF5 file.
    No time or field datasets are created in the file.
    """
    hf = h5py.File(filename, 'w')
    h5p = hf.create_dataset('MeshCornerPoints', data=points, dtype='f8')
    h5p.attrs['Ncp'] = len(points)
    h5p = hf.create_dataset('MeshTriangles', data=triangles, dtype='u4')
    h5p.attrs['Ntri'] = len(triangles)
    h5p = hf.create_dataset('ObservationPosition', data=pos)
    h5p.attrs['Np'] = len(pos)
    hf.close()

def WriteMeshedField(filename, points, triangles, pos, t0, dt, A):
    """
    Write the fields on a meshed geometry to an HDF5 file.
    """
    hf = h5py.File(filename, 'w')
    h5p = hf.create_dataset('MeshCornerPoints', data=points, dtype='f8')
    h5p.attrs['Ncp'] = len(points)
    h5p = hf.create_dataset('MeshTriangles', data=triangles, dtype='u4')
    h5p.attrs['Ntri'] = len(triangles)
    h5p = hf.create_dataset('ObservationPosition', data=pos, dtype='f8')
    h5p.attrs['Np'] = len(pos)
    h5p = hf.create_dataset('ObservationTime',data=t0, dtype='f8')
    h5p.attrs['Nt'] = A.shape[1]
    h5p.attrs['dt'] = dt
    h5p = hf.create_dataset('ElMagField', data=A, dtype='f8')
    hf.close()

def ReadMesh(filename):
    """
    Read the meshed geometry from an HDF5 file.
    No time or field datasets are read.
    Return arrays for mesh corner points, triangle associations and center (trace) positions.
    """
    hdf = h5py.File(filename, "r")
    p = hdf['MeshCornerPoints']
    points = np.array(p)
    p = hdf['MeshTriangles']
    triangles = np.array(p)
    p = hdf['ObservationPosition']
    pos = np.array(p)
    hdf.close()
    return points, triangles, pos

def ReadMeshedField(filename):
    """
    Read the meshed geometry from an HDF5 file.
    Return values:
        - arrayof mesh corner points
        - array of triangle associations
        - array of center (trace) positions
        - array of trace starting times
        - time-step
        - array of fields
    """
    hdf = h5py.File(filename, "r")
    dataset = hdf['MeshCornerPoints']
    points = np.array(dataset)
    dataset = hdf['MeshTriangles']
    triangles = np.array(dataset)
    dataset = hdf['ObservationPosition']
    pos = np.array(dataset)
    dataset = hdf['ObservationTime']
    dt = dataset.attrs.get('dt')
    t0 = np.array(dataset)
    dataset = hdf['ElMagField']
    A = np.array(dataset)
    hdf.close()
    return points, triangles, pos, t0, dt, A
