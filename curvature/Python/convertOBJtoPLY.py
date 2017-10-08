import vtk
import os
import sys

def convertOBJtoPLY(obj_file, ply_file):

    print obj_file
    print ply_file
    
    '''
    (1) Input datafile
    '''
    objReader = vtk.vtkOBJReader()
    objReader.SetFileName(obj_file)
    objReader.Update()
    print "...reading done..."
    
    
    '''
    (2) Output PLY
    '''
    plyWriter = vtk.vtkPLYWriter()
    plyWriter.SetFileName(ply_file)
    plyWriter.SetInputConnection(objReader.GetOutputPort())
    plyWriter.SetFileTypeToBinary()
    plyWriter.SetFileTypeToASCII()
    plyWriter.Write()
    print "...output written..."


if __name__ == "__main__":

    if not len(sys.argv) == 3:
        print ">>>>Error!! Call the script as follows: 'python convertOBJtoPLY.py <OBJ-FILE.PLY> <PLY-FILE.PLY>' "
        exit(1)
        
    obj_file = sys.argv[1]
    ply_file = sys.argv[2]
    
    convertOBJtoPLY(obj_file, ply_file)