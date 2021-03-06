#!/usr/bin/env python
# 
# --------------------------------------------
# Standard python import
from optparse import OptionParser, make_option
import fnmatch, glob, os, sys, json, itertools

def mkReader(chain,redo=False,name="ntuple_reader"):
    from ROOT import gROOT
    if not os.path.isfile("%s.C" % name) or redo:
        chain.MakeClass(name)
    gROOT.LoadMacro("%s.C+" % name)
    
    import ROOT
    return getattr(ROOT,name)(chain)

## ------------------------------------------------------------------------------------------------------------------------------------------------------
def mkChain(files, name="vtxOptTree"):
    from ROOT import TChain
    chain = TChain(name)
    for f in files:
        print "Adding file %s" % f
        chain.AddFile(f)
    return chain

## ------------------------------------------------------------------------------------------------------------------------------------------------------
def getListOfFiles( baseDir, filePattern, expr="{baseDir}/{filePattern}"):
    if baseDir.startswith("/store"):
        return [ "root://eoscms/%s" % f for f in eostools.listFiles( expr.format( baseDir=baseDir, filePattern="" ) )
                 if os.path.basename(f) == filePattern or fnmatch.fnmatch(f,filePattern) ]
    else:
        return glob.glob(expr.format( baseDir=baseDir, filePattern=filePattern))

## ------------------------------------------------------------------------------------------------------------------------------------------------------
def getTMVASettings(cfgs,dest):
    for cfg in cfgs.split(","):
        cf = open(cfg)
        settings = json.load(cf)
        for k,v in settings.iteritems():
            setattr(dest,k,v)
        cf.close()

# Main routine
def main(o,args):

    # Import TMVA classes from ROOT
    from ROOT import TMVA, TFile, TCut 

    print o

    # Output file
    outputFile = TFile( o.outfile % { "label" : o.label }, 'RECREATE' )
    
    factory = TMVA.Factory( "TMVAClassification", outputFile, 
                            "!V:!Silent:!Color:!DrawProgressBar:Transformations=I:AnalysisType=Classification" )

    # Set verbosity
    factory.SetVerbose( o.verbose )
    
    TMVA.Config.Instance().GetIONames().fWeightFileDir = o.weightsdir

    # variables 
    if type(o.variables) == str:
        o.variables = [ v.lstrip().rstrip() for v in o.variables.split(":") if v != "" ]
    allvars = ""
    for v in o.variables:
        factory.AddVariable( str(v),   str(v), "" )
        if allvars != "": allvars += ":"
        allvars += v
    print "variables %s" % allvars

    print o.spectators
    for s in o.spectators:
        if not s in o.variables:
            factory.AddSpectator( str(s),   str(s), "" )

    # categories and sub categories   
    categories = []
    subcategories = [] 
    if hasattr(o,"subcategories") and len(o.subcategories) > 0:
        subcategories = o.subcategories[0]  
        for sc in o.subcategories[1:]:
            subcategories = map( lambda x : ( TCut(x[0][0])*TCut(x[1][0]), "%s_%s" % (x[0][1],x[1][1]) ), itertools.product(subcategories,sc) )
        
    for cut,name,vars in o.categories:
        myvars = allvars
        if vars != "":
            for v in vars.split(":"):
                myvars = myvars.replace(v,"").replace("::",":")
            myvars = myvars.rstrip(":")
            
        vars = str(myvars)
        print vars
        
        if len(subcategories) > 0:
            for subcut,subname in subcategories:
                if subname == "":
                    subname = subname.replace(" ","").replace(">","_gt_").replace("<","_lt_").replace("=","_eq_").replace("&","_and_")
                fullname = "%s_%s" % (name,subname)
                categories.append( ( TCut(cut) * TCut(subcut),str(fullname),vars) )
        else:
            categories.append( (TCut(cut),str(name),vars) )

    # load tree
    background = mkChain(getListOfFiles(o.indir,o.filePattern), o.treename)
    signal = background
    
    # Global event weights (see below for setting event-wise weights)
    signalWeight     = 1.0
    backgroundWeight = 1.0

    factory.AddSignalTree( signal, signalWeight )
    factory.AddBackgroundTree( background, backgroundWeight )

    # sig and bkg definition
    mycutSig = TCut( o.sigcut ) 
    mycutBkg = TCut( o.bkgcut ) 

    # weights
    if hasattr(o,"sigwei"):
        print "reweighing signal"
        factory.AddSpectator( str("sigwei := %s" % o.sigwei) )
        factory.SetSignalWeightExpression    (str(o.sigwei))
    else:
        factory.SetSignalWeightExpression    ("1.")
    if hasattr(o,"bkgwei"):
        print "reweighing background"
        factory.AddSpectator( str("bkgwei := %s" % o.bkgwei) )
        factory.SetBackgroundWeightExpression    (str(o.bkgwei))
    else:
        factory.SetBackgroundWeightExpression    ("1.")

    # "SplitMode=Random" means that the input events are randomly shuffled before
    # splitting them into training and test samples
    factory.PrepareTrainingAndTestTree( mycutSig, mycutBkg,
                                        "nTrain_Signal=0:nTrain_Background=0"
                                        ":NormMode=NumEvents:!V" )
    
    # --------------------------------------------------------------------------------------------------
    # Fisher discriminant (same as LD)
    if "FisherD" in o.methods:
        mname =  "FisherD%s" % o.label
        fcats = factory.BookMethod( TMVA.Types.kCategory, mname )
        
        for cut,name,vars in categories:        
            print "booking sub-category classifier : %s %s %s" % ( cut, name, vars )
            fcats.AddMethod(cut,
                            vars,TMVA.Types.kFisher,"%s_%s" % (mname,name),
                            "!H:!V:Fisher:!CreateMVAPdfs:VarTransform=D"
                            )

    if "Fisher" in o.methods:
        mname =  "Fisher%s" % o.label
        fcats = factory.BookMethod( TMVA.Types.kCategory, mname )
        
        for cut,name,vars in categories:        
            print "booking sub-category classifier : %s %s %s" % ( cut, name, vars )
            fcats.AddMethod(cut,
                            vars,TMVA.Types.kFisher,"%s_%s" % (mname,name),
                            "!H:!V:Fisher:!CreateMVAPdfs"
                            )

    if "Likelihood" in o.methods:
        mname =  "Likelihood%s" % o.label
        fcats = factory.BookMethod( TMVA.Types.kCategory, mname )
        
        for cut,name,vars in categories:        
            print "booking sub-category classifier : %s %s %s" % ( cut, name, vars )
            fcats.AddMethod(cut,
                            vars,TMVA.Types.kLikelihood,"%s_%s" % (mname,name),
                            "!H:!V:!CreateMVAPdfs:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=150"
                            )

    if "LikelihoodD" in o.methods:
        mname =  "LikelihoodD%s" % o.label
        fcats = factory.BookMethod( TMVA.Types.kCategory, mname )
        
        for cut,name,vars in categories:        
            print "booking sub-category classifier : %s %s %s" % ( cut, name, vars )
            fcats.AddMethod(cut,
                            vars,TMVA.Types.kLikelihood,"%s_%s" % (mname,name),
                            "!H:!V:!CreateMVAPdfs:!TransformOutput:VarTransform=D:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=150"
                            )

    if "BDT" in o.methods:
        mname =  "BDT%s" % o.label 
        if len(categories) == 0:
            cats = factory.BookMethod(TMVA.Types.kBDT,mname,
                                      "!H:!V:!CreateMVAPdfs:BoostType=Grad:UseBaggedGrad"
                                      ":GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5"
                                      ":Shrinkage=0.3:NTrees=1000"
                                      )
        else:
            cats = factory.BookMethod( TMVA.Types.kCategory, mname)
            
            for cut,name,vars in categories:
                print "booking sub-category classifier : %s %s %s" % ( cut, name, vars )
                cats.AddMethod(cut,
                               vars,TMVA.Types.kBDT,"%s_%s" % (mname,name),
                               "!H:!V:!CreateMVAPdfs:BoostType=Grad:UseBaggedGrad"
                               ":GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=5"
                               ":Shrinkage=0.3:NTrees=1000"
                               )
    
    # ---- Now you can tell the factory to train, test, and evaluate the MVAs.
    if o.optimize:
        print "Optimizing?"
        factory.OptimizeAllMethods()
        
    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()    
    
    # Save the output.
    outputFile.Close()

## ------------------------------------------------------------------------------------------------------------------------------------------------------    
if __name__ == "__main__":
    parser = OptionParser(option_list=[
            make_option("-i", "--indir",
                        action="store", type="string", dest="indir",
                        default="./",
                        help="input directory", metavar="DIR"
                        ),
            make_option("-f", "--files",
                        action="store", type="string", dest="filePattern",
                        default="vtxOptReduction.root",
                        help="pattern of files to be read", metavar="PATTERN"
                        ), 
            make_option("-t", "--treeName",
                        action="store", type="string", dest="treename",
                        default="vtxOptTree",
                        help="TTree name", metavar="TREENAME"
                        ),
            make_option("-o", "--outfile",
                        action="store", type="string", dest="outfile",
                        default="tmva%(label)s.root",
                        help="outputfile", metavar="FILE"
                        ),
            make_option("-l", "--label",
                        action="store", type="string", dest="label",
                        default="",
                        help="label", metavar="LABEL"
                        ),
            make_option("-w", "--weightsdir",
                        action="store", type="string", dest="weightsdir",
                        default="weights",
                        help="tmva weights dir", metavar="DIR"
                        ),
            make_option("-V", "--variables",
                        action="store", dest="variables", type="string",
                        default="",
                        help="list of variables"
                        ),
            make_option("-T", "--tmvaSettings",
                        action="store", dest="tmvaSettings", type="string",
                        default="pervtx.json",
                        help="settings for the TMVA training"
                        ),
            make_option("-v", "--verbose",
                        action="store_true", dest="verbose",
                        default=False,
                        ),
            make_option("-O", "--optimize",
                        action="store_true", dest="optimize",
                        default=False,
                        ),
            ])

    (options, args) = parser.parse_args()
    getTMVASettings(options.tmvaSettings, options)
        
    main(options, args)
