import vtk
import csv
import sys
import os.path
from array import *
from addCurvValues2PLY import *
from convertOBJtoPLY import *


def processMesh(ply_source,csv_source,csv_dest): 
    '''
    This function loads a PLY-file and the accompanying CSV-file with the curvature information
    and maps cell's (facet) areas to the vertices and creates a new CSV-file with structure:
    kappa_1 -- kappa_2 -- area -- region (from ISD-plot)
    '''

    '''
    (1) Load CSV-file with vertex-curvature-data
    '''
    numarray = array('f', [])
    
    with open(csv_source, 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in spamreader:
            # print row
            numarray.extend((float(row[3]), float(row[4])))
    
    print "...loaded CSV-file and created curvature-array..."
    
    
    '''
    (2) Input PLY-datafile
    '''
    plyReader = vtk.vtkPLYReader()
    plyReader.SetFileName(ply_source)
    plyReader.Update()
    print "...loaded PLY-file..."
    
    mesh_temp = plyReader.GetOutput()  # .GetPointData().AddArray(curv_array);
    
    
    '''
    (2) Initialize and fill "Area" POINT array
    '''
    zeros_array = [0] * mesh_temp.GetNumberOfPoints()
    zeros_array = array('f', zeros_array)
    
    area_array = vtk.vtkFloatArray()
    area_array.SetVoidArray(zeros_array, len(numarray), 1)
    area_array.SetNumberOfComponents(1)
    area_array.SetName('Area')

    mesh_temp.GetPointData().SetScalars(area_array)
    
    print "...created zero-array..."
        
    
    '''
    (3) Mesh quality for retrieving cell areas
    '''
    measure = vtk.vtkMeshQuality()
    measure.SetInputConnection(plyReader.GetOutputPort())
    measure.SetTriangleQualityMeasureToArea()
    measure.Update()
     
    cell_area = measure.GetOutput().GetCellData().GetScalars()
    cell_area.SetName('Facet areas')
    mesh_temp.GetCellData().SetScalars(cell_area)
     
    print "...added cell-areas to cell scalars...now iterating through faces..."

    
    '''
    (4) Iterating through Polygon faces and add third of face area to neighboring vertices
    '''    
    for kk in range(mesh_temp.GetNumberOfCells()):
        area_tmp = mesh_temp.GetCellData().GetScalars().GetTuple(kk)[0]
        for ii in range(3):
            point = mesh_temp.GetCell(kk).GetPointId(ii)
            point_area = float(mesh_temp.GetPointData().GetScalars().GetTuple(point)[0])
            point_area += float(area_tmp)/3
            area_array.SetValue(point, point_area)
            
    print "...mapped all cell areas to vertices areas..."
    
    
    '''
    (5) Write to CSV-file
    '''    
    with open(csv_dest, 'wb') as f:
        writer = csv.writer(f)
       
        for ii in range(mesh_temp.GetNumberOfPoints()):
            writer.writerow([#mesh_temp.GetPoint(ii)[0], \  ## vertex X-coordinate
                             #mesh_temp.GetPoint(ii)[1], \  ## vertex Y-coordinate
                             #mesh_temp.GetPoint(ii)[2], \  ## vertex Z-coordinate
                             numarray[2*ii], \
                             numarray[2*ii+1], \
                             mesh_temp.GetPointData().GetScalars().GetTuple(ii)[0], \
                             writeClassCurv(numarray[2*ii], numarray[2*ii+1]) \
                             ])

    print "...written to CSV..."


if __name__ == "__main__":
    '''
    This whole scirpt loads a OBJ-mesh, transforms it into a PLY-file, then maps the cell (facet)
    areas to vertices and creates a "_NEW.csv" CSV-File (for processing in R) and a new "_curv.ply"
    for plotting in Paraview.
    '''

    if not len(sys.argv) == 3:
        print ">>>>Error!! Call the script as follows: 'python mapCellArea2VertexArea.py <FILE.OBJ> <ORIG-CURV.csv>' "
        exit(1)

    obj_source = sys.argv[1]
    csv_source = sys.argv[2]
    
    ## (1) Transform OBJ to PLY
    print '...STEP 1: Transforming OBJ to PLY...'
    ply_source = obj_source[:-4] + '_fromOBJ.ply'
    if not os.path.isfile(ply_source):
        convertOBJtoPLY(obj_source, ply_source)
    
    ## (2) map CellArea to Vertex Area
    print '...STEP 2: Mapping Cell to Vertex Areas...'
    csv_dest = csv_source[:-4] + '_NEW.csv'
    processMesh(ply_source,csv_source,csv_dest)
    
    ## (3) add curvature values to PLY
    print '...STEP 3: Adding curvatures values to new PLY...'
    target_file = csv_source[:-4] + '_curv.ply'
    addCurvValues2PLY(ply_source, csv_source, target_file)
    
    