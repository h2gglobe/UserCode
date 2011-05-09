import ROOT
# -- Parameters ---------------------------------------------------
current_lumi    = 94.4
file 		= 'Data.root'
plottingdir     = 'plots'
data_name	= 'Data'
file_names	= {
		  'QCD'		:[0,544.862	]
		  ,'GJet'	:[1,2395.580	]
		  ,'Born10'	:[2,2211.781	]
		  ,'Born25'	:[2,24025.257	]
		  ,'Box10'	:[3,2227.736	]
		  ,'Box25'	:[3,62871.867	]
		  ,'Dyee20'	:[4,1636.621	]
		  ,'Dyee10'	:[4,726.965	]
		  ,'WZTTH115'	:[5,40475893.478]
		  ,'VBF115'	:[5,38717590.830]
		  ,'GGH115'	:[5,2848260.736 ]
		  ,'WZTTH120'	:[6,43907180.221]
		  ,'VBF120'	:[6,38470186.499]
		  ,'GGH120'	:[6,2939587.092 ]
		  ,'WZTTH130'	:[7,56914996.108]
		  ,'VBF130'	:[7,41722136.164]
		  ,'GGH130'	:[7,3427365.075 ]
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

