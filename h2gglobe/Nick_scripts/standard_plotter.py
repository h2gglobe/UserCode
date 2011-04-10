import ROOT
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetLineWidth(2)
ROOT.gStyle.SetPalette(1)
import sys,os,getopt

# -- User Imput Options -------------------------------------------
optlist  = ['-signal','-norm','-logy','-nostack','-data']
tmp,opts = getopt.getopt(sys.argv[0:],'', \
			 longopts=optlist)
args = sys.argv[1:]

include_signal 		= False
scale_montecarlo	= False
set_logy		= False
dont_stack		= False	
plot_data		= False
plot_as_gif		= False

for o in opts[1:]:
  if   o=='-sig':
    include_signal = True
    args.remove(o)
  elif o=='-norm':
    scale_montecarlo = True
    args.remove(o)
  elif  o=='-logy':
    set_logy = True
    args.remove(o)
  elif o=='-nostack':
    dont_stack = True
    args.remove(o)
  elif o=='-data':
    plot_data = True
    args.remove(o)
  elif o=='-gif':
    plot_as_gif = True
    args.remove(o)
  else : sys.exit('No Option Available: %s'%o)

# -- Check Sanity -------------------------------------------------
if not plot_data and scale_montecarlo:
  scale_montecarlo = False
if plot_data and not dont_stack:
  plot_ratio = True
else: plot_ratio = False
# -- Parameters ---------------------------------------------------
random		 = ROOT.TRandom3()
current_lumi    = 36.1
data_name 	= 'hist/Data.root'
# ----------------------------------------------------------------

print "Higgs Analysis - Standard Plotter"
print "Plotting all histograms available in ./hist"

file_names	= {
		  'hist/QCD40_hist.root'	:[0,544.862	   	]
		  ,'hist/GJet20_hist.root'	:[1,2395.580	   	]
		  ,'hist/Born_10_hist.root'	:[2,2211.781	   	]
		  ,'hist/Born_25_hist.root'	:[2,24025.257	   	]
		  #,'hist/Born_250_hist.root'	:[2,67685208.127   	]
		  ,'hist/Box_10_hist.root'	:[3,2227.736	   	]
		  ,'hist/Box_25_hist.root'	:[3,62871.867	   	]
		  #,'hist/Box_250_hist.root'	:[3,3795528846.154 	]
		  ,'hist/DYEEM20_hist.root'	:[4,1636.621	   	]
		  ,'hist/DYEEM1020_hist.root'	:[4,726.965	   	]
		  ,'hist/WZTTH115_hist.root'	:[5,40475893.478    ]
		  ,'hist/VBF115_hist.root'	:[5,38717590.830    ]
		  ,'hist/GGH115_hist.root'	:[5,2848260.736     ]
		  ,'hist/WZTTH120_hist.root'	:[6,43907180.221    ]
		  ,'hist/VBF120_hist.root'	:[6,38470186.499    ]
		  ,'hist/GGH120_hist.root'	:[6,2939587.092     ]
		  ,'hist/WZTTH130_hist.root'	:[7,56914996.108    ]
		  ,'hist/VBF130_hist.root'	:[7,41722136.164    ]
		  ,'hist/GGH130_hist.root'	:[7,3427365.075     ]
		  }

#print "Using files/type/weight :"
#for k in file_names.keys(): print k, file_names[k] 

# should have one color and title per Channel
titles_colors = {
		'0':['QCD',ROOT.kRed+1				]
		,'1':['#gamma + Jet',ROOT.kOrange-3		]
	 	,'2':['Born',ROOT.kGreen+3			]
		,'3':['Box',ROOT.kMagenta+3			]
		,'4':['Z/#gamma* #rightarrow ee',ROOT.kOrange+2	]
		,'5':['Higgs M-115',ROOT.kBlue-2		]
		,'6':['Higgs M-120',ROOT.kTeal-2		]
		,'7':['Higgs M-130',ROOT.kCyan-2		]
	 	}

# -- Singal Index indicates which components are signal. -------------------
# -- These wont be included in normalisation		 -------------------
signal_index = [5,6,7]
# --------------------------------------------------------------------------

linked_files = [[filter(lambda x: (file_names[x])[0]==i,file_names),
		 i in signal_index] \
		for i in range(len(file_names))]

#remove signals from file list if no signal option
if not include_signal:
  remove_me=[]
  for si in signal_index: j = titles_colors.pop(str(si))
  for lnk in linked_files:
    if  lnk[1]: remove_me.append(lnk)
  for rm in remove_me:  linked_files.remove(rm)

#print linked_files 
#if not include_signal: file_list = filter(lambda x: not x[1], file_list)

file_list = [[[ROOT.TFile(f) for f in F[0]],F[1]] for F in linked_files]
file_list = filter(lambda x: len(x[0]) > 0, file_list)

# careful here
for f in file_names: 
  if not (os.path.isfile(f)):
   if include_signal and file_names[f][0] in signal_index:
    sys.exit("No file found %s"%f)
   elif not include_signal and file_names[f][0] not in signal_index:
    print file_names[f][0]
    sys.exit("No file found %s"%f)




if os.path.isfile(data_name):  	
  data_file = ROOT.TFile(data_name)
else:
  sys.exit('No file found names %s'%data_name)

if len(titles_colors) != len(file_list): 
  sys.exit("Incorrect number of colors to channels!") 

keys 	  = [[f.GetListOfKeys() for f in F[0]] 
	      for F in file_list]
data_k	  = data_file.GetListOfKeys()

#for i,key_list in enumerate(keys):
#  if len(data_k) != len(key_list):
#    print "In FILE", file_list[i][0].GetName()
#    print "number of data hists:", len(data_k) 
#    print "number of mc   hists:", len(key_list) 
#    sys.exit("Warning, Unmatched histograms!")

tmp=[0 for kk in file_list]

for i in range(len(keys[0][0])):


	leg=ROOT.TLegend(0.55,0.6,0.85,0.85)
	leg.SetFillColor(0)
	leg.SetBorderSize(0)
	
	stack = ROOT.THStack()	
        # ------------------------------------------------------
	c = ROOT.TCanvas()
	if plot_ratio:
          u_pad = ROOT.TPad(str(i)+'u',str(i)+'u'\
			 ,0.01,0.25,0.99,0.99)
          d_pad = ROOT.TPad(str(i)+'d',str(i)+'d'\
			 ,0.01,0.01,0.99,0.25)
	  u_pad.SetNumber(1)
	  d_pad.SetNumber(2)
          u_pad.Draw()
	  d_pad.Draw()

	  u_pad.SetGrid(True)	
	# ------------------------------------------------------

	to_be_stacked = []

	for l in range(len(file_list)):
	  for j in range(len(file_list[l][0])):

		f = file_list[l][0][j]
		h = keys[l][j][i].ReadObj()
		w = (file_names[f.GetName()])[1]

		h.Scale(float(current_lumi)/w)
                print f.GetName()
		if j == 0: 
			tmp[l] = h.Clone()
		else:  tmp[l].Add(h)

	  # -- Only include signal plots if specified -----------
	  print "DEBUG -- Making stacks"
          if   (include_signal and file_list[l][1] ) or\
	   not file_list[l][1]:
   
	    if dont_stack:	     
	      tmp[l].SetLineColor(titles_colors[str(l)][1])
	      tmp[l].SetLineWidth(2)	
	      if  tmp[l].Integral() != 0: tmp[l].Scale(1./tmp[l].Integral())
	      leg.AddEntry(tmp[l],titles_colors[str(l)][0],'L')
	    else:
	      tmp[l].SetFillColor(titles_colors[str(l)][1])
	      leg.AddEntry(tmp[l],titles_colors[str(l)][0],'F')

	    to_be_stacked.append(tmp[l])
	  print "DEBUG -- Made Stacks"
	  # -----------------------------------------------------
	

	data_hist = data_k[i].ReadObj()
	n_data	= data_hist.Integral()

	stacked_sum = 0

	for l,h in enumerate(to_be_stacked):
	  if not file_list[l][1]:
	    stacked_sum+=h.Integral()

	for h in to_be_stacked:
          if stacked_sum != 0 and scale_montecarlo:
	    h.Scale(n_data/stacked_sum)
	  stack.Add(h)  
        
	if set_logy and plot_ratio: u_pad.SetLogy()
        elif set_logy: c.SetLogy()
	else	   : data_hist.SetMinimum(0)
        # -- Data Styles ------------------------
	data_hist.SetMarkerStyle(21)
	data_hist.SetMarkerSize(0.7)
	data_hist.Sumw2()
	# ---------------------------------------	

        c.cd(1)
	if plot_data:
	  leg.AddEntry(data_hist, 'Data','PLE')
	  data_hist.Draw()

	
        if dont_stack:
	 to_be_stacked[0].SetMaximum(1.0)	 
	 to_be_stacked[0].Draw()
         for hist in to_be_stacked[1:]: hist.Draw("same")

	else: 
	  if plot_data: stack.Draw('same')
	  else:  stack.Draw()
	
	
	# Plot data again to ensure overlay
	if plot_data: data_hist.Draw('same')

	leg.Draw()

	if plot_ratio:
	  c.cd(2)
	  hist_rat = data_hist.Clone() 
	  hist_tot = data_hist.Clone() 
	  hist_1 = data_hist.Clone()

	  for b in range(hist_1.GetNbinsX()+1):
	    hist_1.SetBinContent(b,1) 
	    hist_tot.SetBinContent(b, \
		  sum([hh.GetBinContent(b)   \
		  for hh in to_be_stacked]))

          hist_rat.Divide(hist_tot)

	  hist_1.SetLineColor(4)
	  hist_1.SetLineWidth(3)
	  hist_1.SetMaximum(2)
	  hist_1.SetMinimum(0)
	  hist_1.SetTitle('')
	  hist_1.GetYaxis().SetLabelSize(0.12)
	  hist_1.GetYaxis().SetNdivisions(3)
	  hist_1.GetXaxis().SetLabelSize(0.12)
	  hist_1.GetXaxis().SetTitle('')
	  hist_1.GetYaxis().SetTitle('')

	  hist_1.Draw("hist")
	  hist_rat.Draw('same')
        # ----------------------------------------------------
	if plot_as_gif: c.SaveAs('plots/'+data_hist.GetName()+'.gif')
	else: c.SaveAs('plots/'+data_hist.GetName()+'.pdf')


#EOF



  

