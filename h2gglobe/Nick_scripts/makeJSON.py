from JSON import JSON
for p in JSON:
 print "	|| ( run==%d && ("%int(p)
 if len(JSON[p]) ==1:
  l = JSON[p][0]
  print " 	(lumis > %d && lumis < %d)"%(l[0],l[1])
 else:
  l = JSON[p][0]
  print " 	(lumis > %d && lumis < %d)"%(l[0],l[1])
  for l in (JSON[p])[1:]:
         
	 print "  	|| (lumis > %d && lumis < %d)" %(l[0],l[1])
 print "   	)"
 print "       )"
