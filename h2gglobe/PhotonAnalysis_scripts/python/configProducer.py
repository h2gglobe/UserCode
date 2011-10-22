# Python Configuration Handler for Util
# Original Author - Nicholas Wardle

PYDEBUG = 0

# System python imports
import sys,os, json
import array
import ROOT

# Local python imports
from datBlocks import *
from makeFilelist import *
from getTreeEntry import *


def getJson(fname):
    try:
        jstring = open(fname).read()
    except IOError:
        jstring = fname
    return json.loads( jstring )

def defineJsonFilter(fname,dataset):
    goodlumis = getJson( fname )
    
    for run in goodlumis.keys():
        for lumis in goodlumis[run]:
            dataset.addGoodLumi(int(run), int(lumis[0]), int(lumis[1]) )
  

def defineEvList(fname,dataset):
    evlist = getJson( fname )
    
    for run in evlist.keys():
        for ev in evlist[run]:
            dataset.addEventToList(int(run), int(ev[0]), int(ev[1]) )


class configProducer:

  def __init__(self,Ut,conf_filename,Type,njobs=-1,jobId=0):

    print "h2gglobe: step %d, with Config %s. Number of jobs %d. Running job %d" %(Type,conf_filename,njobs,jobId)

    self.ut_   = Ut;
    self.type_ = Type;
    self.is_data_ = True

    self.njobs_ = njobs
    self.jobId_ = jobId
    self.nf_ 	= [0]

    self.sample_weights_file_ = 0
    self.file_processed_events_ = {}

    self.conf_filename = str(conf_filename)
    self.lines_ = []

    self.conf_    = confBlock()
    self.plotvar_ = datBlock()
  
    if self.type_ == 0 or self.type_ == 2:
      self.init_loop()
      self.init_cuts()
      self.init_histos()
      self.init_counters()

    elif self.type_ == 1:
      self.init_reduce()
      self.init_cuts()

    elif self.type_ == 3:
      self.init_loop()

    else: 
      sys.exit("No Such Type As: %d"%self.type_)
      
  def read_weights_file(self):
    weights_lines=[];
    self.read_file(self.sample_weights_file_,weights_lines) 
    for line in weights_lines:
	file,pEvents = line.split("=")
	self.file_processed_events_[file]=int(pEvents)
 
  def read_file(self,conf_filename,lines=None):
    if lines == None:
      self.lines_ = [ ]
      lines = self.lines_
    if not os.path.isfile(conf_filename):
      sys.exit("No Configuration file named - %s"%conf_filename)
    else:
      cf = open(conf_filename,"r")
      comment_status = False
      for l in cf.read().split("\n"):
        line = l.strip(" ").lstrip(" ")
        if len(line) < 2:
          continue
        if line[0:2] == '->':	
          if comment_status:
            comment_status = False
          else:
            comment_status = True
        elif not comment_status:
          lines.append(line)
        else:
          self.conf_.comments+=line

  def init_cuts(self):
    self.read_dat_cuts('cuts.dat')
    for dum in self.plotvar_.vardef:
       self.ut_.AddCut(dum['cutname'],dum['ncat'],dum['dir'],dum['fin'],dum['cutValuel'],dum['cutValueh'])
       
  def init_counters(self):
    self.read_dat_counters('counters.dat')
    self.ut_.InitCounters()
    for dum in self.plotvar_.vardef:
      self.ut_.AddCounter(dum['ncat'] ,dum['countername'], dum['denomname1'], dum['denomname2'], dum['denomname3'])
      
  def init_histos(self):
    self.read_dat_plotvariables('plotvariables.dat')
    self.ut_.InitHistos()
    for dum in self.plotvar_.vardef:
      self.ut_.BookHisto(dum['htyp'],dum['plot'],dum['default'],dum['ncat'],dum['xbins'],dum['ybins'],dum['xmin'],dum['xmax'],dum['ymin'],dum['ymax'],dum['name'])
      
  def init_reduce(self):
    self.read_config_reduce(self.conf_filename)
    if PYDEBUG:
      self.conf_.print_conf()
    self.add_files()
    self.ut_.SetTypeRun(self.type_,self.conf_.outfile)
    # self.ut_.ReadInput(self.type_)
    ## self.read_struct_from_file(self.ut_.vtxAlgoParams,)

  def init_loop(self):
    self.read_config_loop(self.conf_filename)
    if PYDEBUG:
      self.conf_.print_conf()
    self.add_files()
    self.ut_.SetTypeRun(self.type_,self.conf_.histfile)
    for dum in self.conf_.confs:
      dataContainer = self.ut_.DefineSamples(dum['Nam'],dum['typ'],dum['ind'],dum['draw'],dum['red'],dum['tot'],dum['intL'],dum['lum'],dum['xsec'],dum['kfac'],dum['scal'],dum['addnevents'])
      if("json" in dum and dum["json"] != ""):
        defineJsonFilter(dum["json"], dataContainer)
      if("evlist" in dum and dum["evlist"] != ""):
        defineEvList(dum["evlist"], dataContainer)
        
  def add_files(self):
    for t_f in self.conf_.files:
      if (t_f[0]):  # only adding files which aren't Null
        self.ut_.AddFile(t_f[0],t_f[1])
      
  # File Parsers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def read_dat_plotvariables(self,f): 
    print "Parsing of the plotting variables dat"
    self.plotvar_.clear()
    self.read_file(f)
    map_dict = { "htyp" : int, "plot": int, "ncat": int, "xbins" : int,  "ybins": int, "xmin": float, "xmax": float, "ymin": float, "ymax": float, "name": str };
    map_c    = { }
    default = 0
    for line in self.lines_:       
      if len(line) < 2:
        continue
      if (len(line.split()) < 1):
        continue

      # Decide whether this is a define line or a file line:
      if "default" in line:
        split_line = line.split()
        # We have the definiion line
        for sp in split_line:
          val = sp.split("=")
          if val[0] == "default":
            default = int(val[1])
          else:
            sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '"%(val[0],line))

      # We have one of the file def lines
      elif "name" in line:
        split_line = line.split()
        for sp in split_line:
          name, val = sp.split("=")
          if name in  map_dict :
            map_c[name] = map_dict[name](val)
          else:
            sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '"%(name,line))
        self.plotvar_.vardef.append(map_c.copy())
      else:
        sys.exit("Config Line Unrecognised:\n ' %s '"%line)
      for cc in self.plotvar_.vardef:
        cc["default"] = default

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_dat_counters(self,f): 
     print "Parsing of the counters dat"
     self.plotvar_.clear()
     self.read_file(f)
     map_dict = { "ncat":int, "countername": str, "denomname1": str, "denomname2": str, "denomname3": str }
     map_c   = {}
     for line in self.lines_:       
      if len(line) < 2: continue
      if (len(line.split()) < 1): continue

      if "countername" in line:
        # We have one of the file def lines
        split_line = line.split()
        for sp in split_line:
          name, val = sp.split("=")
          if name in  map_dict :
            map_c[name] = map_dict[name](val)
          else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '" %(name,line))

        self.plotvar_.vardef.append(map_c.copy())
      else: sys.exit("Config Line Unrecognised:\n ' %s '"%line)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_dat_cuts(self,f): 
     "Parsing of the cuts dat"
     self.plotvar_.clear()
     self.read_file(f)
     map_c   = {}
     map_dict = { "cutname": str,  "ncat": int, "dir": int, "fin": int }
     for line in self.lines_:       
      if len(line) < 2: continue
      if (len(line.split()) < 1): continue

      if "cutname" in line:
        # We have one of the file def lines
        split_line = line.split()

        val_arr    = []
        minval_arr = []
        maxval_arr = []

        for sp in split_line:
          name, val = sp.split("=")
          if name in  map_dict :
            map_c[name] = map_dict[name](val)
          elif name == "minval":
            minval_arr.append(float(val))
          elif name == "maxval":
            maxval_arr.append(float(val))
          elif name == "val":
            val_arr.append(float(val))
            
          else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '" (val[0],line))
          
        if map_c["dir"] == 2 : 
	   map_c["cutValuel"] = array.array('f',minval_arr)
	   map_c["cutValueh"] = array.array('f',maxval_arr)
	else: 
	   map_c["cutValuel"] = array.array('f',val_arr)
	   map_c["cutValueh"] = array.array('f',[float(-9999.)\
					         for v in val_arr])
	
        self.plotvar_.vardef.append(map_c.copy())
      else: sys.exit("Config Line Unrecognised:\n ' %s '"%line)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_config_reduce(self,f):
     "Parsing of the reduction configuration"        
     self.read_file(f)
     comment_status = False 

     for line in self.lines_:
       # Decide whether this is a define line or a file line:
       if "output=" in line:
         self.read_output_file(line)

       elif "Fil=" in line or "Dir=" in line:
         self.read_input_files_reduce(line)

       # load analyzer
       elif line.startswith("analyzer "):
         self.read_analyzer( line )  
         
       # input and output branches  
       elif line.startswith("inputBranches "):
         self.read_input_branches(line)

       # input and output branches
       elif line.startswith("outputBranches "):
         self.read_output_branches(line)

       # Read a generic member of the LoopAll class
       else:
         self.read_struct_line(line,self.ut_)
           
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  def generate_weights_file(self,filename):
     if len(self.file_processed_events_):
       wfile = open(filename,"w")
       for samplefile in self.file_processed_events_:
	 wfile.write("%s=%d\n"%(samplefile,self.file_processed_events_[samplefile]) )
       wfile.close()
       print "Written Number Of Processed Events file -- ",filename

  def read_config_loop(self,f):
     "Parsing of the looper configuration"        
     self.read_file(f)
     self.intL	    = 0
     # look for a file called f.weights
     if os.path.isfile(f+".pevents") : 
	self.sample_weights_file_= f+".pevents"
	self.read_weights_file()

     for line in self.lines_:
       # Decide whether this is a define line or a file line:
       if "histfile" in line:  
         self.read_histfile(line)

       elif "Fil=" in line or "Dir=" in line:
         self.read_input_files_loop(line)

       # load analyzer
       elif line.startswith("analyzer "):
         self.read_analyzer( line )  
         
       # input and output branches  
       elif line.startswith("inputBranches "):
         self.read_input_branches(line)
         
       # input and output branches      
       elif line.startswith("outputBranches "):
         self.read_output_branches(line)
         
       # Read a generic member of the LoopAll class
       else:
         self.read_struct_line(line,self.ut_)
     
     for cc in  self.conf_.confs:
       cc["intL"] = self.intL 

     if  self.sample_weights_file_==0:
	 if self.njobs_==-1 or self.jobId_==0:
	    self.generate_weights_file(f+".pevents")
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_struct_line(self,line,struct):
    split_line = line.split()
    struct_name = split_line.pop(0)

    # check if the line is a simple parameter
    if "=" in struct_name:
      name,val = struct_name.split("=")
      t = type( struct.__getattribute__(name) )
      struct.__setattr__(name, t(val) )
      print "%s.%s = %s" % ( struct.__name__, name, str(val) )
      return

    # otherwise it is complex structure
    struct = struct.__getattribute__(struct_name)
    try:
      if os.path.isfile(split_line[0]):
        self.read_struct_from_file(split_line[0],struct)
      else:
        self.read_struct(split_line,struct)
    except Exception, e:
      sys.exit("Unkown data member %s at line:\n\n%s\n%s"% (struct_name, line, e) )


  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_struct(self,lines,struct):
    """Read a structure from a file
    """
    if PYDEBUG:
      print "Reading structure %s" % ( str(struct) )
    for line in lines:
      split_line = line.split()

      # check if this line specify a complex structure
      if not "=" in split_line[0]:
        self.read_struct_line(line,struct)
        
      # otherwise read all parameters
      else:
        for sp in split_line:
          name,val = [ s.lstrip(" ").rstrip(" ") for s in sp.split("=") ]
          try:
	    if ".root" in val and not "/castor" in val and not os.path.isfile(val): sys.exit("No File found - %s, check the line %s"%(val,line))
            t = type( struct.__getattribute__(name) )
            struct.__setattr__(name, t(val) )
            print "%s = %s" % ( name, str(struct.__getattribute__(name)) )
          except AttributeError:
            sys.exit( "Error: unkown attribute: %s\nline: %s\n%s" % ( name, line, struct ) )
          except ValueError, e:
            sys.exit( "Syntax error in line: %s\n%s" % ( line, e ) )
    

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_struct_from_file(self,f,struct):
    lines = []
    if PYDEBUG:
      print "Reading structure from file %s %s " % ( str(struct), f )
    self.read_file(f,lines)
    
    self.read_struct(lines,struct)
   

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_histfile(self,line):
    # We have the definition line           
    split_line = line.split()
    for sp in split_line:
      val = sp.split("=")
      
      if val[0] == "output":
        self.conf_.outfile = str(val[1])
        if self.njobs_ > 0:
            name,ext=self.conf_.outfile.rsplit(".",1)
            self.conf_.outfile = "%s_%d.%s" % ( name, self.jobId_, ext )

      elif val[0] == 'histfile':
        self.conf_.histfile = str(val[1])
        if self.njobs_ > 0:
            name,ext=self.conf_.histfile.rsplit(".",1)
            self.conf_.histfile = "%s_%d.%s" % ( name, self.jobId_, ext )
        
      elif val[0] == 'intL':
        self.intL = float(val[1])
      else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '" %(val[0],line))
      
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_analyzer(self,line):
    ## a, name, config = line.split()
    sl = line.split()
    sl.pop(0)
    name = sl.pop(0)
    analyzer = ROOT.__getattr__(name)()
    print "Loading analyzer %s " % name
    self.ut_.AddAnalysis( analyzer )
    for config in sl:
        try:
            if os.path.isfile(config):
                self.read_struct_from_file(config,analyzer)
            else:
                self.read_struct(config,analyzer)
        except Exception, e:
            sys.exit("Unable to read analyzer %s at line:\n\n%s\n%s"% (name, line, e) )
    print "Loaded analyzer %s " % name
    
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_output_file(self,line):
    split_line = line.split()
    # We have the definition line
    for sp in split_line:
      val = sp.split("=")
      if val[0] == "output":
        outfile = str(val[1])
        if self.njobs_ > 0:
            name,ext=outfile.rsplit(".",1)
            outfile = "%s_%d.%s" % ( name, self.jobId_, ext )
        print "Outfile: %s" % outfile
        self.conf_.outfile = outfile
        outdir=outfile.rsplit("/",1)[0]
        if outdir != outfile:
          try:
            os.system( "mkdir -p %s" % outdir )
          except Exception, e:
            print e
      else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '" %(val[0],line))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_input_files_reduce(self,line):
    values = { "CaDir" : "","DcDir" : "", "Dir" : "", "typ" : -1, "Fil" : ""  }; 
    # We have one of the file def lines
    split_line = line.split()
    for sp in split_line:
      name,val = sp.split("=")
      if name in values:
        values[name] = type(values[name])(val)
      else:
        sys.exit( "Unrecognised Argument:\n ' %s ' in line:\n ' %s '"% (name,line) )
    directory = values["Dir"]
    cas_directory = values["CaDir"]
    dcs_directory = values["DcDir"]
    fi_name   = values["Fil"]
    fi_type   = values["typ"]
    if fi_type != 0:
      self.is_data_ = False
    
    if fi_name != '':
      if not fi_name.startswith("rfio") and not os.path.isfile(fi_name): 
        sys.exit("No Input File Named: %s"%fi_name)
      tuple_n = fi_name, fi_type
      self.conf_.files.append(tuple_n)
        
    if cas_directory != '':
      ca_files = makeCaFiles(cas_directory,self.njobs_,self.jobId_)
      for file_s in ca_files:
	if file_s[1]:self.conf_.files.append((file_s[0],fi_type))
	else	    :self.conf_.files.append((None,fi_type))

    if dcs_directory != '':
      dc_files = makeDcFiles(dcs_directory,self.njobs_,self.jobId_)
      for file_s in dc_files:
	if file_s[1]:self.conf_.files.append((file_s[0],fi_type))
	else	    :self.conf_.files.append((None,fi_type))

    if directory != '':
        files = makeFiles(directory,self.njobs_,self.jobId_)
        for file_s in files:
	  if file_s[1]:self.conf_.files.append((file_s[0],fi_type))
	  else	    :self.conf_.files.append((None,fi_type))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_input_files_loop(self,line):
    print "read_input_files_loop"
    map_c = {"typ":-99999,"Nam":"default","draw":-999,"ind":-999,"tot":0,"red":-999,"lum":1.0,"xsec":1.0,"kfac":1.0,"scal":1.0,"json":"","evlist":""}
    #map_c["tot"]=-1
    map_c["addnevents"]=0
    directory = ''
    cas_directory = ''
    dcs_directory = ''
    fi_name   = ''
    fi_type   = -99999
    # We have one of the file def lines
    split_line = line.split()
    for sp in split_line:
      val = sp.split("=")
      if val[0] == "Fil":
        fi_name = str(val[1])
      elif val[0]== "Dir":
        directory=str(val[1])
      elif val[0]== "CaDir":
        cas_directory=str(val[1])
      elif val[0]== "DcDir":
        dcs_directory=str(val[1])
      elif val[0] == "typ":
        fi_type = int(val[1])
        map_c["typ"] = int(val[1])
      elif val[0] in map_c:
        map_c[val[0]] = type(map_c[val[0]])(val[1])
      else: sys.exit("Unrecognised Argument:\n ' %s ' in line:\n ' %s '"
                     %(val[0],line))
    if map_c["typ"] != 0:
      self.is_data_ = False
      
    if fi_name != '':
      if not fi_name.startswith("rfio") and not os.path.isfile(fi_name): 
        sys.exit("No Input File Named: %s"%fi_name)
      tuple_n = fi_name, fi_type
      self.nf_[0]+=1
      if not ( (self.njobs_>0) and (self.nf_[0] % self.njobs_ != self.jobId_) ):
        self.conf_.files.append(tuple_n)
      else: self.conf_.files.append((None,fi_type));
      if fi_type!=0 and fi_type!=-99999 and map_c["tot"] == 0:
	if self.sample_weights_file_==0 :
	  nEventsInFile = getTreeEntry(fi_name,"global_variables","processedEvents")
	  self.file_processed_events_[fi_name] = nEventsInFile
          map_c["tot"] = nEventsInFile;
	  
	else:  
	  if fi_name in self.file_processed_events_:map_c["tot"] = self.file_processed_events_[fi_name]
	  else: (sys.exit("No Entry for %s found in %s, Please Delete %s and re-run with option --dryRun to regenerate it"%(fi_name,self.sample_weights_file_,self.sample_weights_file_)))

	map_c["addnevents"] = int(1)
      self.conf_.confs.append(map_c.copy())

    mkFiles = None
    dir = None
    if cas_directory != '':
        mkFiles = makeCaFiles
        dir = cas_directory  
    if dcs_directory != '':
        mkFiles = makeDcFiles
        dir = dcs_directory  
    if directory != '':
        mkFiles = makeFiles
        dir = directory  
        
    if dir:
      files = mkFiles(dir,self.njobs_,self.jobId_,self.nf_)
      if fi_type!=0 and fi_type!=-99999 and map_c["tot"] == 0:
          allfiles = mkFiles(dir,-1,-1)
          for file_s in allfiles:
	      print "Getting N Processed Events for - ", file_s[0]
	      if self.sample_weights_file_==0 :
		nEventsInFile = getTreeEntry(file_s[0],"global_variables","processedEvents")
                map_c["tot"] = map_c["tot"] + nEventsInFile
		self.file_processed_events_[file_s[0]] = nEventsInFile
	      else:
		if file_s[0] in self.file_processed_events_: map_c["tot"] = map_c["tot"] + self.file_processed_events_[file_s[0]]
		else: (sys.exit("No Entry for %s found in %s, Please Delete %s and re-run with option --dryRun to regenerate it"%(file_s[0],self.sample_weights_file_,self.sample_weights_file_)))

      for file_s in files:
	  if file_s[1]: self.conf_.files.append((file_s[0],fi_type))
	  else:      self.conf_.files.append((None,fi_type))
          self.conf_.confs.append(map_c.copy())

          
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_input_branches(self,line):
    print line
    a, list = line.split(" ")
    branches = []
    if os.path.isfile(list):
      self.read_file(list,branches)
    else:
      branches = list.split(",")
    for b in branches:
      print b
      bname = b
      btype = 0
      if ":" in b:
        bname, sbtype = b.split(":")
        btype = int(sbtype)
      #if btype == 0 or not self.is_data_:
      print "Reading input branch %s " % bname   
      self.ut_.InputBranch(bname,btype)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  def read_output_branches(self,line):
    a, list = line.split(" ")
    branches = []
    if os.path.isfile(list):
      self.read_file(list,branches)
    else:
      branches = list.split(",")
    for b in branches:
      ## all outputs need to be read  first
      bname = b
      btype = 0
      if ":" in b:
        print b
        bname, sbtype = b.split(":")
        btype = int(sbtype)
     # if btype == 0 or not self.is_data_:
      self.ut_.InputBranch(bname,btype)
      self.ut_.OutputBranch(bname)

#--------------------------------------------------------//
#EOF 
