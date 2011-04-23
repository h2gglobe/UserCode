from JSON import JSON

outfile=open("goodRuns_cpp.txt",'w')

print "Parsing JSON.py"
counts = 0
for p in JSON:
 if counts == 0: 
  outfile.write("if(   ( run==%d && (\n"%int(p))
 else:
  outfile.write("   || ( run==%d && (\n"%int(p))
 if len(JSON[p]) ==1:
   l = JSON[p][0]
   outfile.write( " 	(lumis > %d && lumis < %d)\n"%(l[0],l[1]))
 else:
   l = JSON[p][0]
   outfile.write( " 	(lumis > %d && lumis < %d)\n"%(l[0],l[1]))
   for l in (JSON[p])[1:]:
	 outfile.write( "  	|| (lumis > %d && lumis < %d)\n" %(l[0],l[1]))
 outfile.write( "   	)\n")
 outfile.write( "  	)\n")
 counts+=1

outfile.write( "  )\n")
print "%d Runs saved to goodRuns.txt"%counts
