import os
import sys
import subprocess
import Image


if len(sys.argv) < 2:
    print ">>>>Error!! Call the script as follows: 'python createRaw.py /<Folder-to-TIFs>/' [ /<DESTINATION-DIR>/ ]"
    exit(1)
src = sys.argv[1]

# List all tif-images from directory
imgs = [name for name in os.listdir(src)
                  if name.lower().endswith('.tif') and not name.startswith('.')]
imgs.sort()
len_str = str(len(imgs))

# Get image sizes
im = Image.open(src+imgs[0])
xx = str(im.size[0])
yy = str(im.size[1])
zz = len_str

# Split to get parent directory and directory base
sample_name = os.path.split(src[:-1])[1]
if not len(sys.argv) == 3:
    parent_dir = os.path.split(src[:-1])[0]
else:
    parent_dir = sys.argv[2]
    parent_dir = parent_dir[:-1]
    
print parent_dir

# (1) convert with Fiji
cmd1 = 'open='+src+imgs[0]+' number='+len_str+' starting=1 increment=1 scale=100 file=[] sort'
cmd2 = parent_dir+'/'+sample_name+'_'+xx+'_'+yy+'_'+zz+'.raw'
cmd = 'fiji -eval \"run(\\"Image Sequence...\\", \\"'+cmd1+'\\");'
cmd += 'saveAs(\\"Raw Data\\",\\"'+cmd2+'\\");'
cmd += 'run(\\"Quit\\");'  # close Fiji again
cmd += '\"'
print cmd

# (2) close Fiji
subprocess.check_call('fiji -eval \"run(\\"Quit\\");\"', shell=True)

# (3) OPEN + SAVE + CLOSE Fiji     
subprocess.check_call(cmd, shell=True)
print 'finish'
