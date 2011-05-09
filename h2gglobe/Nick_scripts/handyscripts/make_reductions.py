import commands,sys,datetime
dir = str(sys.argv[1])
now = datetime.datetime.now()

sc,flist = commands.getstatusoutput('srmls $DCACHE_SRM_ROOT/%s'%(dir))

dirs = flist.split('\n')
dirs = dirs[4:]  # silly dCache

runfile = open("run_all_mc.sh","w")
runfile.write("#!/bin/bash")
added_files = []

if not sc:
  for f in dirs:
  
	if len(f) < 1: continue
	f = (f.split('/'))[-2]
	file = open("%s.dat"%f,"w")
	added_files.append("%s.dat"%f)
	file.write("output=/vols/cms02/h2g/ntuples/mc_V06/%s.root\n"%f)
	file.write("typ=1 DcDir=%s\n"%dir)
	file.write("->\n")
	file.write("Produced at : %s\n"%str(now))
	file.write("->\n")

else : sys.exit("No such directory as %s"%dir)

for add in added_files:
  runfile.write("python reduce_input.py %s\n"%add)
