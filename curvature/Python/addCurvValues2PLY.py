import os
import sys
import numpy
import csv

def addCurvValues2PLY(ply_file, csv_file, target_file):
    '''
    (0) Read CSV-file
    '''
    # csv_lines = tuple(open(csv_file, 'r'))
    # numpy.loadtxt(open(csv_file,"rb"),delimiter=",",skiprows=0)
    csv_data = csv.reader(open(csv_file, "rb"), delimiter=',')
    print "...loaded CSV-data..."
#     x = list(csv_data)
#     print "...casted CSV to LIST..."
#     curv_data = numpy.array(x).astype('float')
#     print "...assigned to NUMPY array..."


    '''
    (1) Read PLY-header
    '''
    ply_lines = tuple(open(ply_file, 'r'))

    kk = 0
    header = []

    while ply_lines[kk] != 'end_header\n':
        header.append(ply_lines[kk])
        kk = kk + 1
    print "...loaded PLY-header..."


    '''
    (2) Get constants from PLY-file
    '''
    for kk, line in enumerate(header):
        if line[0:14] == 'element vertex':
            n_vertices = line[15:-1]
            pos_vertices = kk
        elif line[0:12] == 'element face':
            n_faces = line[13:-1]
            pos_faces = kk
    print "...loaded constants from PLY-header..."


    '''
    (3) Add new properties to header
    '''
    header.insert(pos_faces, 'property float nz\n')
    header.insert(pos_faces, 'property float ny\n')
    header.insert(pos_faces, 'property float nx\n')


    '''
    (4) Write new Header to file
    '''
    with open(target_file, 'w') as new_ply:
        for line in header:
            new_ply.write(line)
        new_ply.write('end_header\n')
    print "...wrote new PLY-header to new PLY-file...now iterating through vertices..."


    '''
    (5) Write vertex data
    '''
    with open(target_file, 'a') as new_ply:
        for kk in range(int(n_vertices)):
            newline = ply_lines[len(header) + kk - 2].strip()  # # -2 because of new header properties
            curv_data = csv_data.next()
            # (1) Write curvature classifiers:
            newline += ' ' + writeClassCurv(curv_data[3], curv_data[4])
            newline += ' ' + curv_data[3]
            newline += ' ' + curv_data[4] + '\n'

            # (2) Write 3 curvatures (mean, min, max):
            # newline += writeAllThreeCurv(curv_data[3], curv_data[4])

            new_ply.write(newline)
    print "...vertices written...now writing rest of file..."


    '''
    (6) Write rest of data
    '''
    rest_ply = ply_lines[len(header) - 3 + int(n_vertices) + 1:]
    with open(target_file, 'a') as new_ply:
        for line in rest_ply:
            new_ply.write(line)

    print "DONEEEE"
    # print ply_lines[kk+1]
    # print header


def writeClassCurv(k_min, k_max):
    '''
    This function writes curvatures according to the ISD-regions, i.e.
    the curvature classifiers.
    '''
    if float(k_min) >= 0.0 and float(k_max) >= 0.0:
        newline_tmp = '0.0'  # + str((float(curv_data[kk, 3]) + float(curv_data[kk, 4])) / 2.0)
    elif float(k_min) < 0.0 and float(k_max) >= 0.0 \
    and -float(k_min) <= float(k_max):
        newline_tmp = '0.25'  # + str((float(curv_data[kk, 3]) + float(curv_data[kk, 4])) / 2.0)
    elif float(k_min) < 0.0 and float(k_max) >= 0.0 \
    and -float(k_min) > float(k_max):
        newline_tmp = '0.5'  # + str((float(curv_data[kk, 3]) + float(curv_data[kk, 4])) / 2.0)
    elif float(k_min) < 0.0 and float(k_max) < 0.0:
        newline_tmp = '0.75'  # + str((float(curv_data[kk, 3]) + float(curv_data[kk, 4])) / 2.0)
    else:
        print "apparently missed a vertex"

    return newline_tmp


def writeAllThreeCurv(k_min, k_max):
    '''
    This function writes all the 3 curvatures standardly
    '''
    newline_tmp = ' ' + str((float(k_min) + float(k_max)) / 2.0)
    newline_tmp += ' ' + str(float(k_min))
    newline_tmp += ' ' + str(float(k_max)) + '\n'

    return newline_tmp



if __name__ == "__main__":

    if not len(sys.argv) == 4:
        print ">>>>Error!! Call the script as follows: 'python addCurvValues2PLY.py <FILE.PLY> <CURV.csv> <TARGET-FILE.PLY>' "
        exit(1)

    ply_file = sys.argv[1]
    csv_file = sys.argv[2]
    target_file = sys.argv[3]

    addCurvValues2PLY(ply_file, csv_file, target_file)
