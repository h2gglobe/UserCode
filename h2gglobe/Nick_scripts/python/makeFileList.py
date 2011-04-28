import commands,sys

def makeDcFiles(dir):

   dcache_prepend = 'dcap://gfe02.grid.hep.ph.ic.ac.uk:22128//pnfs/hep.ph.ic.ac.uk/data/cms'
   dir = str(dir)
   return_files = []

   sc,flist = commands.getstatusoutput('srmls $DCACHE_SRM_ROOT/%s'%(dir))

   if not sc:
    files = flist.split('\n')
    for f in files:
        if len(f) < 1: continue
        f = (f.split()[-1]).split('/')[-1]
	if '.root' in f:
	  return_files.append(dcache_prepend+dir+'/'+f)

   else:
    sys.exit("No Such Directory: %s"%(dir))
   return return_files


