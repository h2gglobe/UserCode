import ROOT
import array,os,sys
import datetime
from optparse import OptionParser

#-------------------------------------------------------------------------
# UserInput
parser=OptionParser()
parser.add_option("-i","--input",dest="fileName",help="Input file name")
parser.add_option("-o","--output",dest="outputDir",default="h2g",help="Output directory for plots")
parser.add_option("-d","--diagnostics",action="store_false",default=True,help="Turn diagnostic plot making off")
parser.add_option("-w","--www",action="store_false",default=True,help="Turn webpage making off")
parser.add_option("-W","--htmlOnly",action="store_true",default=False,help="Run only html making")
parser.add_option("-B","--backgroundOnly",action="store_true",default=False,help="Run correction to background model only")
parser.add_option("-I","--sigInterpOnly",action="store_true",default=False,help="Run signal interpolation only")
parser.add_option("-D","--datacardsOnly",action="store_true",default=False,help="Run datacard creation only")
parser.add_option("","--atCERN",action="store_true",default=False,help="need to source setupROOT csh if at CERN")
(options,args)=parser.parse_args()
  
interpFileName = options.fileName+"_interpolated.root"
#-------------------------------------------------------------------------
if options.atCERN :os.system("source setupROOT.csh")
else:os.system("source setupROOT.sh")
ROOT.gROOT.ProcessLine(".L createCorrectedBackgroundModel.C+g")
os.system("cmsenv")
ROOT.gROOT.ProcessLine(".L BDTInterpolation.C+g")

from ROOT import createCorrectedBackgroundModel
from ROOT import BDTInterpolation

ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)

nSidebands=6
datacardFile="mva-datacards"
mvaPlotFile="mva-plots"


if not options.htmlOnly:
  if not options.sigInterpOnly and not options.datacardsOnly:
    # First correct background model
    print '----------------------------------------------------------'
    print '   Corrected Background Model '
    print '----------------------------------------------------------'
    createCorrectedBackgroundModel(options.fileName,nSidebands,options.diagnostics)

  if not options.backgroundOnly and not options.datacardsOnly:
    # Second interpolate intermediate signal distributions
    print '----------------------------------------------------------'
    print '   Interpolating Signal Model to intermediate masses '
    print '----------------------------------------------------------'
    BDTInterpolation(options.fileName,options.diagnostics,True,True)

  if not options.backgroundOnly and not options.sigInterpOnly:
    # Third write datacards
    print '----------------------------------------------------------'
    print '   Writing datacards for limit '
    print '----------------------------------------------------------'
    bashCommand = "python writeBinnedMvaCard.py -i "+interpFileName +" --makePlot -t"
    os.system(bashCommand+" ada")
    os.system(bashCommand+" grad")
# Now make html pages
if options.www or options.htmlOnly:
  print '----------------------------------------------------------'
  print '   Publishing diagnostics to web '
  print '----------------------------------------------------------'
  now = datetime.datetime.now()
  temp = now.strftime("%Y-%m-%d %H:%M")
  timestring=temp.replace(' ','_')
  os.system("whoami > temp.txt")
  tempF = open('temp.txt')
  user = tempF.readlines()[0]

  os.system("cp mva-plots-ada/model* BMplots/ada")
  os.system("cp mva-plots-grad/model* BMplots/grad")
  os.system("cp mva-plots-ada/diff_model* BMplots/ada")
  os.system("cp mva-plots-grad/diff_model* BMplots/grad")

  os.system("python make_html.py "+interpFileName)
  os.system("python make_bkg_html.py "+interpFileName)

  outputLocation = "~/public_html/"+options.outputDir+"/"+timestring
  os.system("mkdir -p "+outputLocation)
  os.system("cp -r plots "+outputLocation)
  os.system("cp -r BMplots "+outputLocation)
  os.system("ln -s "+outputLocation+" ~/public_html/"+options.outputDir+"/Latest")
  print "Plots available to view in ~/public_html/"+options.outputDir+"/Latest/"
  print "If working on /vols/ at IC plots can be view at: \n \t \t \t www.hep.ph.ic.ac.uk/~"+user+"/"+options.outputDir+"/Latest/plots/plots.html"

