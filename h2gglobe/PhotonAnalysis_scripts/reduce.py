import ROOT
import sys, os, commands
#from readJSON import ReadJSON # if python version >=2.6
from evalJSON import ReadJSON # if python version < 2.6

##~~~~~~~~ User Defined Things
dir = "../samples/"
##~~~~~~~~

ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("libRooFit.so");
ROOT.gSystem.Load("libLoopAll.so");

def chunks(l, n):
  for i in xrange(0, len(l), n):
    yield l[i:i+n]
def prompt_file_overwrite(path):
  if os.path.exists(path):
    print "File %s exists. Overwriting..." % path
    #print "File %s exists. Do you want to overwrite (default n)" % path
    return True#(raw_input("\ty/n:").lstrip().rstrip()=="y")
    #return (raw_input("\ty/n:").lstrip().rstrip()=="y")
  return True

##~~~~~~~~ Defaults
doBatch = False
nBatch = 10
queue = 'hepmedium.q'
isJob = False
nJob = -1
##~~~~~~~~

args = sys.argv[1:]#0 is filename

if len(args)<1 :
  raise NameError("Incorrect Usage (-help for help)")

for i, arg in enumerate(args) : 
  if "-b" in arg:
    doBatch = True
  elif "-n" == arg:
    nBatch = int(args[i+1])
  elif "-n" in arg:
    nBatch = int(arg[2:])
  elif "-q" == arg:
    queue = args[i+1]
  elif "-q" in arg:
    queue= arg[2:]
  elif "-j" == arg:
    isJob = True
    nJob = int(args[i+1])
  elif "-j" in arg:
    isJob = True
    nJob = int(arg[2:])

jsonFile = str(args[-1])
json = ReadJSON(dir+jsonFile)

if len(json) < 1:
  raise NameError("Error with JSON file.")


if doBatch :##SUBMIT TO BATCH QUEUE AND EXIT
  print "Setting up batch analysis."
  jobTemplate = '''#!/bin/sh
cd H2GDIR
source setup.sh
cd THISDIR
python ./THISFILE -j${SGE_TASK_ID} -nNUMBATCH  SAMPLE
  '''
  output = jobTemplate
  output = output.replace("H2GDIR",os.environ.get('H2G_WORKING_SW_DIR'))
  output = output.replace("THISDIR",sys.path[0])
  output = output.replace("THISFILE",sys.argv[0])
  output = output.replace("SAMPLE",jsonFile)
  output = output.replace("NUMBATCH",str(nBatch))
  jobFile = sys.argv[0].split(".")[0]+"_"+jsonFile.split(".")[0]+".sh"

  if prompt_file_overwrite(jobFile):
    print "Making job file: %s" % jobFile
  else : sys.exit(4)
  f=open(jobFile,"w")
  f.write(output)
  f.close()
  nJobs = len(list(chunks(json['files'],nBatch)))
  print "Submitting "+str(nJobs)+" jobs to batch"

  submit = commands.getstatusoutput("qsub -q %s -t 1-%i:1  %s " % (queue, nJobs, jobFile))
  #Succesful?
  if submit[0] != 0 :
    print "Error occured when submitting jobs:", submit[0]
    print submit[1]     
    sys.exit(3)
  else :
     print submit[1]
  sys.exit(0)

elif isJob: #RUN OVER A SUBSET OF FILES (ie this is an individual job in a parametric job)
  splitFiles = list(chunks(json['files'],nBatch))
  fileList = splitFiles[nJob-1]#only run over given files TODO?
  outputName = json['name']+"_"+str(nJob-1)+".root"

else :     ##RUN OVER ALL FILES
  fileList = json['files']
  outputName =   json['name']+".root"

print "Output file: ", outputName
ROOT.gBenchmark.Start("Reduction")
ut = ROOT.Util()
ut.SetTypeRun(1, outputName)
print "Files:"
for file in fileList :
  print file
  ut.AddFile(str(file),1)
print "Starting analysis."
ut.LoopAndFillHistos();
ROOT.gBenchmark.Show("Reduction")
print "Writing output."
ut.WriteHist();  
print "Done." 

