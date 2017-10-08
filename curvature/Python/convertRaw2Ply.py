import vtk
import os
import sys
import re


if not len(sys.argv) == 2:
    print ">>>>Error!! Call the script as follows: 'python convertRaw2Ply.py <FILE_x_y_z.raw>' "
    exit(1)
    
src = sys.argv[1]
splittet_dir = os.path.split(src)
resolution = re.split('_', splittet_dir[1])
xx = int(resolution[1])-1
yy = int(resolution[2])-1
zz = int(resolution[3][:-4])-1

print xx, yy, zz

'''
(1) Input datafile
'''
  
reader = vtk.vtkImageReader2()
reader.SetFileName(src)
reader.SetDataScalarTypeToUnsignedChar()
reader.SetDataByteOrderToBigEndian()
reader.SetHeaderSize(0)
reader.SetNumberOfScalarComponents(1)
reader.SetFileDimensionality(3)
reader.SetDataOrigin(0,0,0)
reader.SetDataSpacing(1,1,1)
reader.SetDataExtent(0, xx, 0, yy, 0, zz)
reader.FileLowerLeftOn()
reader.Update()
 
datas=reader.GetOutput()
print "...reading done..."
  
'''
(2) filter
'''
npoints=datas.GetNumberOfPoints()
print npoints
 
# Contour
#contour = vtk.vtkContourFilter()
contour = vtk.vtkMarchingContourFilter()
contour.SetInputData(datas)
#contour.GenerateTrianglesOn()
# contour.ComputeNormalsOn()
contour.SetValue( 0, 128 )
contour.Update()
print contour.GetNumberOfContours()
print "...contour done..."
  
'''
(3) Output PLY
'''
filename = splittet_dir[0]+'/'+resolution[0]+'.ply'

plyWriter = vtk.vtkPLYWriter()
plyWriter.SetFileName(filename)
plyWriter.SetInputConnection(contour.GetOutputPort())
plyWriter.SetFileTypeToBinary()
# plyWriter.SetFileTypeToASCII()
plyWriter.Write()
print "...output written..."