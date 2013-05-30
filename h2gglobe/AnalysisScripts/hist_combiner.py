#!/usr/bin/env python

import os, sys, glob

taskdir=sys.argv.pop(1)
filestocomb=glob.glob("%s/filestocombine_*.dat"%taskdir)
print 'Getting files from ', filestocomb

opts = ""
yes = False
for a in sys.argv[1:]:
  if a == "-y":
    yes = True
  else:
    opts += "%s " % a
  
listofhists=" "

if len(filestocomb)==1:
  cfile = open( filestocomb[0],"r")
  in_comment = False
  for line in cfile.read().split('\n'):
    if line.startswith("->"):
      in_comment = not in_comment
    if line.startswith("#") or in_comment: 
      continue
    longl = line.split('Fil=')
    if len(longl)>1:
      path = os.path.dirname(longl[1])+"/histograms_"+os.path.basename(longl[1])
      if 'castor' in path:
        listofhists += "rfio:"+path+" "
      else:
        listofhists += path+" "

cmd = ("hadd %s -f "+taskdir+"/histograms_CMS-HGG.root"+listofhists) % opts
print "Will execute: \n"
print "\t\t %s" % cmd

if not yes:
  print "Happy to hadd?\n"
  raw_input()
  
os.system( cmd )





