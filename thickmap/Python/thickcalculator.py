import subprocess
import os
import sys
import psutil
from time import sleep


class LocalThickness(object):
    '''
    This is the class for calculating Local Thicknesses with Fiji.
    It stores variables that are necessary for communicating with
    Fiji over the bash environment. 
    '''
    def __init__(self,img_path):
        self.img_path = os.path.join(str(img_path), '')  # make sure "/" is at end
        
        ## Get name of first file
        tif_list = [name for name in os.listdir(self.img_path)
                    if name.lower().endswith('.tif') and
                    not name.startswith('.')]
        tif_list.sort()
        self.tif_file = tif_list[0]
        
        ## Get folder name
        parts = self.splitOsPath(self.img_path)
        self.dir_name = parts[-2]
        
        ## Command line string
        self.cmd_string = ''
        
        ## Fiji process that is the one launched in the beginning
        self.parent_proc = []
        
    
    def loadImageSequenceInFiji(self):        
        #run_str = '\"run(\\"Image Sequence...\\", \\"open=' + self.img_path + self.tif_file + ' sort\\");\"'
        run_str = 'run(\\"Image Sequence...\\", \\"open=' + self.img_path + self.tif_file + ' sort\\");'
        self.cmd_string = self.cmd_string + run_str
                
        
    def runLocalThicknessCalculation(self):
        #run_str = '\"run(\\"Local Thickness (complete process)\\", \\"threshold=128\\");\"'
        run_str = 'run(\\"Local Thickness (complete process)\\", \\"threshold=128\\");'
        
        self.cmd_string = self.cmd_string + run_str
        #self.selectFijiWindow(self.dir_name + '_EDT_DR_LT')
                
        
    def selectFijiWindow(self,window_name):
        run_str = 'selectWindow(\\"' + window_name + '\\");'
        self.cmd_string = self.cmd_string + run_str
        
        
    def saveRAWData(self):
        path_list = self.splitOsPath(self.img_path)
        parent_path = self.glueOsPath(path_list[:-2])
                
        new_RAWfile = os.path.join(str(parent_path), self.dir_name+'_LocThk.raw')
        
        run_str = 'saveAs(\\"Raw Data\\", \\"' + new_RAWfile + '\\");'
        self.cmd_string += run_str
        
    
    def evaluateFijiCMD(self):
        cmd = 'fiji --ij2 -eval ' + '\"' + self.cmd_string + '\"' ##run(\\"Quit\\");
        print cmd
        
        #return
        proc = psutil.Popen(cmd, shell=True, executable='/bin/bash')
        
        if not self.parent_proc:
            self.parent_proc = proc
            print self.parent_proc.pid
        #proc.communicate()
        #p = psutil.Process(proc.pid)
        
        # First the command is sent to Fiji, we wait the CPU to power up
        cpu_percent = proc.cpu_percent()
        while cpu_percent < 2.0:
            cpu_percent = proc.cpu_percent()
        else:
            print "...Fiji process started..."
            
        
        # CPU should now be more than 1 %. If it is less than 1 % for more than 2s,
        # we end the process...Additionally we also monitor the parent process from
        # the beginning, as it seems to run during the "Local Thickness Calculation"
        cpu_percent_parent = 10.0 #self.parent_proc.cpu_percent()
        
        while not cpu_percent < 4.0 or not cpu_percent_parent < 4.0:
            cpu_percent = 0.0
            cpu_percent_parent = 0.0
            for x in range(4):
                sleep(2)
                #print x
                cpu_percent += proc.cpu_percent()
                cpu_percent_parent += self.parent_proc.cpu_percent()
            #print cpu_percent
        else:
            print "...Fiji process ended..."
        
    
    def splitOsPath(self,path):
        '''
        This method splits any OS path to a list of folders and should
        be compatible on all platforms.
        Taken from Grecoman: https://github.com/gnudo/grecoman
        '''
        parts = []
        while True:
            path_tmp, rest = os.path.split(path)
            if path_tmp == path:
                assert not rest
                if path:
                    parts.append(path)
                break
            parts.append(rest)
            path = path_tmp
        parts.reverse()
        return parts
    
    def glueOsPath(self, pathlist):
        '''
        This method is the reverse to the "splitOsPath"
        Taken from Grecoman: https://github.com/gnudo/grecoman
        '''
        newpath = pathlist[0]
        for item in pathlist[1:]:
            newpath = os.path.join(str(newpath), str(item))
        return newpath


if __name__ == '__main__':
    if not len(sys.argv) == 2:
        print ">>>>Error!! Call the script as follows: 'python thickcalculator.py <FOLDER-OF-semented-tif-imgs>' "
        exit(1)
    
    src = sys.argv[1]
    
    thick = LocalThickness(src)
    
    thick.loadImageSequenceInFiji()
    thick.evaluateFijiCMD()
    
    thick.cmd_string = 'run(\\"Geometry to Distance Map\\", \\"threshold=128\\");'
    thick.evaluateFijiCMD()
    
    thick.cmd_string = ''
    thick.selectFijiWindow(thick.dir_name)
    thick.evaluateFijiCMD()
    
    thick.cmd_string = 'close();'
    thick.evaluateFijiCMD()
    
    thick.cmd_string = ''
    thick.selectFijiWindow(thick.dir_name + '_EDT')
    thick.evaluateFijiCMD()
    
    thick.cmd_string = 'run(\\"Distance Map to Distance Ridge\\");'
    thick.evaluateFijiCMD()
    
    thick.cmd_string = ''
    thick.selectFijiWindow(thick.dir_name + '_EDT')
    thick.evaluateFijiCMD()
    
    thick.cmd_string = 'close();'
    thick.evaluateFijiCMD()
    
    thick.cmd_string = ''
    thick.selectFijiWindow(thick.dir_name + '_EDT_DR')
    thick.evaluateFijiCMD()
    
    thick.cmd_string = 'run(\\"Distance Ridge to Local Thickness\\");'
    thick.evaluateFijiCMD()
    
    thick.cmd_string = ''
    thick.selectFijiWindow(thick.dir_name + '_EDT_DR')
    thick.evaluateFijiCMD()
    
    thick.cmd_string = 'close();'
    thick.evaluateFijiCMD()
    
    thick.cmd_string = ''
    thick.selectFijiWindow(thick.dir_name + '_EDT_DR_LT')
    thick.evaluateFijiCMD()
    
    thick.cmd_string = 'run(\\"Local Thickness to Cleaned-Up Local Thickness\\");'
    thick.evaluateFijiCMD()
    
    thick.cmd_string = ''
    thick.selectFijiWindow(thick.dir_name + '_EDT_DR_LT')
    thick.evaluateFijiCMD()
    
    thick.cmd_string = 'close();'
    thick.evaluateFijiCMD()
    
    thick.cmd_string = ''
    thick.selectFijiWindow(thick.dir_name + '_EDT_DR_LT_CL')
    thick.evaluateFijiCMD()
    
    thick.cmd_string = ''
    thick.saveRAWData()
    thick.evaluateFijiCMD()
    
    thick.cmd_string = 'close();'
    thick.evaluateFijiCMD()
    
    thick.cmd_string = 'run(\\"Quit\\");'
    thick.evaluateFijiCMD()
    
#     thick = LocalThickness(src)
#     
#     # (1) Load Image sequence
#     thick.loadImageSequenceInFiji()
#     
#     # (2) Calculate local thickness
#     thick.runLocalThicknessCalculation()
#     
#     # (3) Evaluate Fiji Command
#     thick.evaluateFijiCMD()
#     
#     # (4) Save RAW file
#     print "now quitting"
#     thick.cmd_string = ''
#     thick.saveRAWData()
#     #thick.cmd_string += 'run(\\"Quit\\");'
#     thick.evaluateFijiCMD()
    
    
    
    
    