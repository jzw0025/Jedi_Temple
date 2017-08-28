from vtk import *
# The source file
file_name = "/Users/junchaowei/Desktop/Jedi_Temple/CA3_template/solution.vtk"

# Read the source file.
reader = vtkUnstructuredGridReader()
reader.SetFileName(file_name)
reader.Update() # Needed because of GetScalarRange
output = reader.GetOutput()
scalar_range = output.GetScalarRange()

# Create the mapper that corresponds the objects of the vtk file
# into graphics elements
mapper = vtkDataSetMapper()
mapper.SetInputData(output)
mapper.SetScalarRange(scalar_range)

# Create the Actor
actor = vtkActor()
actor.SetMapper(mapper)

# Create the Renderer
renderer = vtkRenderer()
renderer.AddActor(actor)
renderer.SetBackground(1, 1, 1) # Set background to white

# Create the RendererWindow
renderer_window = vtkRenderWindow()
renderer_window.AddRenderer(renderer)

# Create the RendererWindowInteractor and display the vtk_file
interactor = vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderer_window)
interactor.Initialize()
interactor.Start()
