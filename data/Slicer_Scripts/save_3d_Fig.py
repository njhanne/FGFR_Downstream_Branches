# copy paste this into the interpreter in Slicer
# it will change the background to black and save it as a .png without any background

### Don't foget to change the filename and path! ###

# Set background to black (required for transparent background)
view = slicer.app.layoutManager().threeDWidget(0).threeDView()
view.mrmlViewNode().SetBackgroundColor(0,0,0)
view.mrmlViewNode().SetBackgroundColor2(0,0,0)
# Capture RGBA image
renderWindow = view.renderWindow()
renderWindow.SetAlphaBitPlanes(1)
wti = vtk.vtkWindowToImageFilter()
wti.SetInputBufferTypeToRGBA()
wti.SetInput(renderWindow)
writer = vtk.vtkPNGWriter()
writer.SetFileName("C:/Users/njhan/Desktop/screenshot.png")
writer.SetInputConnection(wti.GetOutputPort())
writer.Write()