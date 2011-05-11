#define PADEBUG 0

void LoopAll::TermRealPhotonAnalysis(Util* ut, int typerun) 
{
   if (typerun==3){	
      rooContainer->FitToData("cms-data_model","Data_mass");
      rooContainer->FitToData("bw"	      ,"Signal_mass");
      rooContainer->FitToData("exp"	      ,"Background_mass",100,115,125,200);
      rooContainer->FitToSystematicSet("bw","Signal_mass","e-scale");
      rooContainer->FitToSystematicSet("exp","Background_mass","e-scale");

      std::string outputfilename = (std::string) ut->histFileName;  // Not happy about doing this!!!

      rooContainer->WriteDataCard(outputfilename,"Data_mass","Signal_mass","Background_mass");
   }

}

void LoopAll::InitRealPhotonAnalysis(int typerun) {
  if(PADEBUG) 
    cout << "InitRealPhotonAnalysis START"<<endl;


  if (typerun == 2 || typerun == 1){

  }
  
  if (typerun == 3) {  
     // initialize the RooContainer with number of categories

     rooContainer = new RooContainer(4);

     //RooFitting type
     rooContainer->AddRealVar("Signal_mass",100.,200.);
     rooContainer->AddRealVar("Background_mass",100.,200.);
     rooContainer->AddRealVar("Data_mass",100.,200.);

     // Signal Parameters
     rooContainer->AddRealVar("width",1.,0.01,5.);
     rooContainer->AddRealVar("mean",130.,100.,150.);
     rooContainer->AddRealVar("sigma",1.,0.01,5.);

     // Background Parameters
     rooContainer->AddRealVar("c",+3.0,-5.,+5);
     rooContainer->AddRealVar("p1",-0.03,-1.,-0.001);
     rooContainer->AddRealVar("p2",+0.01,-1.,+1);

     rooContainer->AddRealVar("alpha",-0.03,-1.,-0.001);
     // All parameters for the Data
     rooContainer->AddRealVar("d_width",1.,0.1,3.);
     rooContainer->AddRealVar("d_mean",130.,100.,150.);
     rooContainer->AddRealVar("d_alpha",-0.04,-1.,-0.001);

     // Set Up The Signal Pdfs
     // -------------------------------------//
     std::vector<std::string> sig_pars(3,"p");	 
     sig_pars[0] = "Signal_mass";
     sig_pars[1] = "mean";
     sig_pars[2] = "width";
     rooContainer->AddGenericPdf("bw","bw",sig_pars,3,10.);
    
     std::vector<std::string> bkg_pars(2,"p");	 
     bkg_pars[0] = "Background_mass";
     bkg_pars[1] = "alpha";
     rooContainer->AddGenericPdf("exp"
	,"exp(@1*@0)",bkg_pars,0,100);
     // -------------------------------------//
     // Set Up The Background Pdfs

     // A Model For The Data
     // -------------------------------------//
     std::vector<std::string> data_sig(3,"p");
     data_sig[0] = "Data_mass";
     data_sig[1] = "d_mean";
     data_sig[2] = "d_width";

     std::vector<std::string> data_bkg(2,"p");
     data_bkg[0] = "Data_mass";
     data_bkg[1] = "d_alpha";

     rooContainer->AddGenericPdf("data_bw"
	,"1./((@0-@2)*(@0-@2)+@1*@1)"
	,data_sig,3,1.);
     rooContainer->AddGenericPdf("data_exp"
	,"exp((@0)*(@1))",data_bkg,1,10);
     // -------------------------------------//
     // Combine the Data pdfs
     std::vector<std::string> pdfs(2,"t");
     pdfs[0] = "data_exp";
     pdfs[1] = "data_bw";
     rooContainer->ComposePdf("cms-data_model"
			     ,"data_bkg+data_gaus"
			     ,pdfs);
     // -------------------------------------//
     // Create The DataSets
     rooContainer->CreateDataSet("Data_mass",50);
     rooContainer->CreateDataSet("Signal_mass",50);
     rooContainer->CreateDataSet("Background_mass",50);

     // Systematic Errors
     rooContainer->MakeSystematics("Signal_mass","e-scale");
     rooContainer->MakeSystematics("Background_mass","e-scale");
     rooContainer->MakeSystematics("Data_mass","e-scale");
  }

  if(PADEBUG) 
    cout << "InitRealPhotonAnalysis END"<<endl;

}

bool LoopAll::myPhotonAnalysisRunSelection(int cur_type){

  if (cur_type != 0)    return true;


  bool passes_run_selection = false; 

if (run < 163217) return true;  // This is the end of 2010 Runs + paolo
				 // -> JSON already applied

if(   ( run==163659 && (
 	(lumis > 1 && lumis < 374)
  	|| (lumis > 376 && lumis < 650)
  	|| (lumis > 652 && lumis < 705)
   	)
  	)
   || ( run==163658 && (
 	(lumis > 1 && lumis < 3)
   	)
  	)
   || ( run==163657 && (
 	(lumis > 1 && lumis < 140)
   	)
  	)
   || ( run==163630 && (
 	(lumis > 76 && lumis < 164)
  	|| (lumis > 176 && lumis < 185)
   	)
  	)
   || ( run==163655 && (
 	(lumis > 15 && lumis < 23)
   	)
  	)
   || ( run==163480 && (
 	(lumis > 1 && lumis < 92)
  	|| (lumis > 96 && lumis < 188)
  	|| (lumis > 190 && lumis < 191)
   	)
  	)
   || ( run==163478 && (
 	(lumis > 1 && lumis < 70)
   	)
  	)
   || ( run==163297 && (
 	(lumis > 1 && lumis < 219)
   	)
  	)
   || ( run==163296 && (
 	(lumis > 59 && lumis < 230)
  	|| (lumis > 232 && lumis < 585)
   	)
  	)
   || ( run==163481 && (
 	(lumis > 1 && lumis < 72)
  	|| (lumis > 74 && lumis < 77)
  	|| (lumis > 79 && lumis < 79)
   	)
  	)
   || ( run==163663 && (
 	(lumis > 1 && lumis < 106)
  	|| (lumis > 109 && lumis < 246)
   	)
  	)
   || ( run==163334 && (
 	(lumis > 1 && lumis < 35)
  	|| (lumis > 37 && lumis < 37)
  	|| (lumis > 156 && lumis < 556)
   	)
  	)
   || ( run==163252 && (
 	(lumis > 60 && lumis < 137)
   	)
  	)
   || ( run==163585 && (
 	(lumis > 1 && lumis < 32)
   	)
  	)
   || ( run==163337 && (
 	(lumis > 1 && lumis < 18)
  	|| (lumis > 27 && lumis < 201)
  	|| (lumis > 203 && lumis < 426)
  	|| (lumis > 434 && lumis < 461)
   	)
  	)
   || ( run==163583 && (
 	(lumis > 1 && lumis < 63)
  	|| (lumis > 65 && lumis < 92)
  	|| (lumis > 96 && lumis < 155)
  	|| (lumis > 157 && lumis < 173)
  	|| (lumis > 175 && lumis < 219)
   	)
  	)
   || ( run==163582 && (
 	(lumis > 1 && lumis < 22)
   	)
  	)
   || ( run==163332 && (
 	(lumis > 43 && lumis < 118)
  	|| (lumis > 224 && lumis < 264)
  	|| (lumis > 266 && lumis < 599)
  	|| (lumis > 601 && lumis < 639)
  	|| (lumis > 641 && lumis < 801)
   	)
  	)
   || ( run==163333 && (
 	(lumis > 1 && lumis < 106)
   	)
  	)
   || ( run==163661 && (
 	(lumis > 1 && lumis < 17)
   	)
  	)
   || ( run==163338 && (
 	(lumis > 1 && lumis < 164)
   	)
  	)
   || ( run==163339 && (
 	(lumis > 1 && lumis < 172)
   	)
  	)
   || ( run==163589 && (
 	(lumis > 1 && lumis < 49)
  	|| (lumis > 51 && lumis < 160)
   	)
  	)
   || ( run==163588 && (
 	(lumis > 1 && lumis < 8)
  	|| (lumis > 10 && lumis < 446)
   	)
  	)
   || ( run==163235 && (
 	(lumis > 1 && lumis < 461)
   	)
  	)
   || ( run==163234 && (
 	(lumis > 1 && lumis < 66)
   	)
  	)
   || ( run==163237 && (
 	(lumis > 1 && lumis < 213)
   	)
  	)
   || ( run==163374 && (
 	(lumis > 1 && lumis < 599)
  	|| (lumis > 603 && lumis < 863)
   	)
  	)
   || ( run==163375 && (
 	(lumis > 1 && lumis < 10)
   	)
  	)
   || ( run==163233 && (
 	(lumis > 1 && lumis < 283)
   	)
  	)
   || ( run==163232 && (
 	(lumis > 110 && lumis < 149)
   	)
  	)
   || ( run==163378 && (
 	(lumis > 1 && lumis < 81)
  	|| (lumis > 89 && lumis < 272)
  	|| (lumis > 306 && lumis < 615)
   	)
  	)
   || ( run==163270 && (
 	(lumis > 1 && lumis < 76)
  	|| (lumis > 79 && lumis < 96)
  	|| (lumis > 99 && lumis < 475)
  	|| (lumis > 479 && lumis < 527)
  	|| (lumis > 529 && lumis < 685)
  	|| (lumis > 695 && lumis < 928)
   	)
  	)
   || ( run==163482 && (
 	(lumis > 1 && lumis < 27)
  	|| (lumis > 48 && lumis < 48)
   	)
  	)
   || ( run==163238 && (
 	(lumis > 9 && lumis < 15)
   	)
  	)
   || ( run==163358 && (
 	(lumis > 39 && lumis < 63)
   	)
  	)
   || ( run==163476 && (
 	(lumis > 1 && lumis < 94)
  	|| (lumis > 98 && lumis < 212)
   	)
  	)
   || ( run==163587 && (
 	(lumis > 1 && lumis < 52)
   	)
  	)
   || ( run==163664 && (
 	(lumis > 1 && lumis < 119)
  	|| (lumis > 121 && lumis < 178)
   	)
  	)
   || ( run==163370 && (
 	(lumis > 1 && lumis < 147)
   	)
  	)
   || ( run==163402 && (
 	(lumis > 37 && lumis < 582)
  	|| (lumis > 586 && lumis < 801)
   	)
  	)
   || ( run==163376 && (
 	(lumis > 1 && lumis < 20)
  	|| (lumis > 22 && lumis < 246)
   	)
  	)
   || ( run==163660 && (
 	(lumis > 1 && lumis < 74)
   	)
  	)
   || ( run==163586 && (
 	(lumis > 1 && lumis < 75)
   	)
  	)
   || ( run==163371 && (
 	(lumis > 1 && lumis < 107)
  	|| (lumis > 148 && lumis < 363)
   	)
  	)
   || ( run==163584 && (
 	(lumis > 1 && lumis < 56)
   	)
  	)
   || ( run==163385 && (
 	(lumis > 52 && lumis < 240)
  	|| (lumis > 244 && lumis < 406)
   	)
  	)
   || ( run==163668 && (
 	(lumis > 1 && lumis < 53)
  	|| (lumis > 57 && lumis < 136)
  	|| (lumis > 140 && lumis < 213)
   	)
  	)
   || ( run==163387 && (
 	(lumis > 1 && lumis < 256)
   	)
  	)
   || ( run==163757 && (
 	(lumis > 1 && lumis < 40)
   	)
  	)
   || ( run==163289 && (
 	(lumis > 1 && lumis < 388)
   	)
  	)
   || ( run==163261 && (
 	(lumis > 1 && lumis < 3)
  	|| (lumis > 10 && lumis < 126)
   	)
  	)
   || ( run==163255 && (
 	(lumis > 1 && lumis < 359)
  	|| (lumis > 412 && lumis < 844)
  	|| (lumis > 846 && lumis < 846)
  	|| (lumis > 848 && lumis < 977)
   	)
  	)
   || ( run==163475 && (
 	(lumis > 30 && lumis < 295)
   	)
  	)
   || ( run==163662 && (
 	(lumis > 1 && lumis < 154)
   	)
  	)
   || ( run==163286 && (
 	(lumis > 112 && lumis < 401)
   	)
  	)
   || ( run==163483 && (
 	(lumis > 1 && lumis < 57)
   	)
  	)
   || ( run==163738 && (
 	(lumis > 34 && lumis < 311)
   	)
  	)
   || ( run==163340 && (
 	(lumis > 1 && lumis < 488)
   	)
  	)
   || ( run==163369 && (
 	(lumis > 1 && lumis < 94)
   	)
  	)
   || ( run==163301 && (
 	(lumis > 1 && lumis < 192)
   	)
  	)
   || ( run==163300 && (
 	(lumis > 1 && lumis < 616)
   	)
  	)
   || ( run==163479 && (
 	(lumis > 1 && lumis < 175)
   	)
  	)
   || ( run==163302 && (
 	(lumis > 1 && lumis < 190)
   	)
  	)
   || ( run==163596 && (
 	(lumis > 1 && lumis < 29)
   	)
  	)
   || ( run==163372 && (
 	(lumis > 1 && lumis < 52)
   	)
  	)
  )
  // Put here the Run Selection you want in the Data

  passes_run_selection = true; 

  // -----------------------------------------------

  return passes_run_selection;

}


void LoopAll::myGetEntryPhotonRedAnalysis(Util *ut, int jentry, int cur_type){

  if (cur_type ==0){
    b_event->GetEntry(jentry);
    b_lumis->GetEntry(jentry);
    b_run->GetEntry(jentry);
  }

  b_pho_n->GetEntry(jentry); 
  b_pho_p4->GetEntry(jentry); 
  b_pho_r9->GetEntry(jentry); 
  b_pho_calopos->GetEntry(jentry); 
  b_pho_hoe->GetEntry(jentry); 
  b_pho_sieie->GetEntry(jentry); 
  b_pho_ecalsumetconedr03->GetEntry(jentry); 
  b_pho_ecalsumetconedr04->GetEntry(jentry); 
  b_pho_hcalsumetconedr03->GetEntry(jentry); 
  b_pho_hcalsumetconedr04->GetEntry(jentry); 
  b_pho_trksumptsolidconedr03->GetEntry(jentry); 
  b_pho_trksumpthollowconedr04->GetEntry(jentry); 
  b_pho_haspixseed->GetEntry(jentry);

  if (cur_type!=0){

    b_gen_n->GetEntry(jentry);
    b_gen_p4->GetEntry(jentry);
    b_gen_pdgid->GetEntry(jentry);
    b_gen_status->GetEntry(jentry);
    b_gen_mother->GetEntry(jentry);
  } 
  

}

void LoopAll::myFillHistPhotonAnalysis(Util* ut, int jentry) {

  if(PADEBUG) 
    cout << "myFillHist START"<<endl;

  
  for (int i=0; i<pho_n; i++) {
    TLorentzVector *p4 = (TLorentzVector *) pho_p4->At(i);
    FillHist("pho_pt", p4->Pt());
  }
  
  Int_t in_endcap = 0;
  Float_t best_mass = 0;
  for (int i=0; i<pho_n-1; i++) {
    TLorentzVector *pg1= (TLorentzVector *) pho_p4->At(i);
    if (fabs(pg1->Eta()) > 1.479)
      in_endcap = 1;

    for (int j=i+1; j<pho_n; j++) {
      TLorentzVector *pg2= (TLorentzVector *) pho_p4->At(j);
      if (fabs(pg2->Eta()) > 1.479)
	in_endcap = 1;
      TLorentzVector higgs = (*pg1) + (*pg2);
      Float_t mass = higgs.M();
      if (mass > best_mass)
	best_mass = mass;
    }
  }     

  if (best_mass != 0) 
    FillHist("invmass", best_mass);
  
  if(PADEBUG) 
    cout<<"myFillHist END"<<endl;
}


void LoopAll::myFillHistPhotonAnalysisRed(Util * ut, int jentry) {

  if(PADEBUG) 
    cout << "myFillHistRed START"<<endl;

  counters[0]++;
  std::vector<PhotonCandidate> preselected_photons;  

  TVector3 *calopos;	
  TLorentzVector *p4;

  for (int i=0; i<pho_n; i++) {
    p4 = (TLorentzVector *) pho_p4->At(i);
    calopos  = (TVector3 *) pho_calopos->At(i);
    float pt  = p4->Pt(); 
    float eta = fabs(calopos->Eta());

    //Photon Selection
     if ( 
       (! pho_haspixseed[i])
       && pt > 30. 
       && pho_hoe[i] <  0.02
       && pho_trksumpthollowconedr04[i] < (1.5 + 0.001*pt)
       && pho_ecalsumetconedr04[i] < (2.0 + 0.006*pt)
       && pho_hcalsumetconedr04[i] < (2.0 + 0.0025*pt)
       && (   ((eta < 1.4442) && pho_sieie[i] < 0.01)
	   || ((eta > 1.566) && (eta < 2.5) && pho_sieie[i] < 0.028)
	  )
       ) {
         PhotonCandidate candidate;
         candidate.p4 		= p4;
	 candidate.calopos	= calopos;
         candidate.pixSeed 	= pho_haspixseed[i];
         candidate.trkIso 	= pho_trksumpthollowconedr04[i];
         candidate.ecalIso 	= pho_ecalsumetconedr04[i];
         candidate.hcalIso 	= pho_hcalsumetconedr04[i];
         candidate.sieie 	= pho_sieie[i];
         candidate.hoe 		= pho_hoe[i];
         candidate.r9 		= pho_r9[i];
         preselected_photons.push_back(candidate);
       }
  }

  //Event Selection
  int n_preselected_pho = preselected_photons.size();

  FillHist("n_sel",0,n_preselected_pho);

  // Sort Photons into Pt Order
  std::sort(preselected_photons.begin()
           ,preselected_photons.end()
           ,PhoP4greater); 
 
// Regular Event Selection begins here
  float best_mass = -1;
  float best_pt   = -1;
  int   category  = -1;
  float min_r9;
  float max_eta;

  if (n_preselected_pho > 1 ){
     PhotonCandidate leading, nleading;

     leading  = preselected_photons[0];
     nleading = preselected_photons[1];


     if (leading.p4->Pt() > 40.)  {
         min_r9  = min(leading.r9
		      ,nleading.r9);
	 max_eta = max(fabs(leading.calopos->Eta())
			   ,fabs(nleading.calopos->Eta()));

	 if (min_r9 < 0.93 && max_eta < 1.4442 ) category = 1;
	 if (min_r9 > 0.93 && max_eta < 1.4442 ) category = 2;
	 if (min_r9 < 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 3;
	 if (min_r9 > 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 4;

         // -------------------------------------------------------
         TLorentzVector Higgs = (*(leading.p4))
                                 +(*(nleading.p4));
           float mass = Higgs.M();
           float h_pt = Higgs.Pt();

             //if (mass > 100. && mass < 150.){

	       FillHist("pho_pt",0,leading.p4->Pt());
	       FillHist("pho_pt",0,nleading.p4->Pt());
	       FillHist("pho_pt",category,leading.p4->Pt());
	       FillHist("pho_pt",category,nleading.p4->Pt());

               best_mass = mass;
 	       best_pt   = h_pt;

            //}
     }
   }

  FillHist("mass",0, best_mass);
  FillHist("pt",0, best_pt);

  if (category > -1){
    FillHist("mass",category, best_mass);
    FillHist("pt",category, best_pt);
  }

  if(PADEBUG) 
    cout<<"myFillHistRed END"<<endl;
 
}

void LoopAll::myStatPhotonAnalysis(Util * ut, int jentry) {

  if(PADEBUG) 
    cout << "myStat START"<<endl;
  counters[0]++;

  float weight = sampleContainer[ut->current_type_index].event_weight;
  int cur_type = ut->current_type;

  std::vector<PhotonCandidate> preselected_photons;  

  TVector3 *calopos;	
  TLorentzVector *p4;

  for (int i=0; i<pho_n; i++) {
    p4 = (TLorentzVector *) pho_p4->At(i);
    calopos  = (TVector3 *) pho_calopos->At(i);
    float pt  = p4->Pt(); 
    float eta = fabs(calopos->Eta());

    //Photon Selection
     if ( 
       (! pho_haspixseed[i])
       && pt > 30. 
       && pho_hoe[i] <  0.02
       && pho_trksumpthollowconedr04[i] < (1.5 + 0.001*pt)
       && pho_ecalsumetconedr04[i] < (2.0 + 0.006*pt)
       && pho_hcalsumetconedr04[i] < (2.0 + 0.0025*pt)
       && (   ((eta < 1.4442) && pho_sieie[i] < 0.01)
	   || ((eta > 1.566) && (eta < 2.5) && pho_sieie[i] < 0.028)
	  )
       ) {
         PhotonCandidate candidate;
         candidate.p4 		= p4;
	 candidate.calopos	= calopos;
         candidate.pixSeed 	= pho_haspixseed[i];
         candidate.trkIso 	= pho_trksumpthollowconedr04[i];
         candidate.ecalIso 	= pho_ecalsumetconedr04[i];
         candidate.hcalIso 	= pho_hcalsumetconedr04[i];
         candidate.sieie 	= pho_sieie[i];
         candidate.hoe 		= pho_hoe[i];
         candidate.r9 		= pho_r9[i];
         preselected_photons.push_back(candidate);
       }
  }

  //Event Selection
  int n_preselected_pho = preselected_photons.size();
  // Sort Photons into Pt Order
  std::sort(preselected_photons.begin()
           ,preselected_photons.end()
           ,PhoP4greater); 
 
// Regular Event Selection begins here
  float best_mass = 0.;
  float best_pt   = -1;
  int   category  = -1;
  float min_r9;
  float max_eta;

  if (n_preselected_pho > 1 ){
     PhotonCandidate leading, nleading;

     leading  = preselected_photons[0];
     nleading = preselected_photons[1];


     if (leading.p4->Pt() > 40.)  {

         min_r9  = min(leading.r9
		      ,nleading.r9);
	 max_eta = max(fabs(leading.calopos->Eta())
			   ,fabs(nleading.calopos->Eta()));
	 if (min_r9 < 0.93 && max_eta < 1.4442 ) category = 0;
	 if (min_r9 > 0.93 && max_eta < 1.4442 ) category = 1;
	 if (min_r9 < 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 2;
	 if (min_r9 > 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 3;

         // -------------------------------------------------------
	
         TLorentzVector Higgs = (*(leading.p4))
                                 +(*(nleading.p4));
           float mass = Higgs.M();

	   // Want to Reacalculate the Mass to account for ES - uncertainty
	   // Energy Scale Error 1-sigma
	   float sys_error = 0.05;

	   std::vector<float> mass_errors;
           for (int sys=1;sys<31;sys++){
	        float incr = (float)sys/10;
		TLorentzVector lead_err = (1.-incr*sys_error)*(*leading.p4);
		TLorentzVector nlead_err = (1.-incr*sys_error)*(*nleading.p4);
                mass_errors.push_back((lead_err+nlead_err).M());
   	   }
           for (int sys=1;sys<31;sys++){
	        float incr = (float)sys/10;
		TLorentzVector lead_err = (1.+incr*sys_error)*(*leading.p4);
		TLorentzVector nlead_err = (1.+incr*sys_error)*(*nleading.p4);
                mass_errors.push_back((lead_err+nlead_err).M());
   	   }
            

  	   if (cur_type ==0){        // Data
	     rooContainer->InputDataPoint("Data_mass",category,mass);

           } else if (cur_type < 0){ // Signal
	     rooContainer->InputDataPoint("Signal_mass",category,mass,weight);
             rooContainer->InputSystematicSet("Signal_mass","e-scale",category,mass_errors);

	   } else if (cur_type > 0){ // Background
	     rooContainer->InputDataPoint("Background_mass",category,mass,weight);
             rooContainer->InputSystematicSet("Background_mass","e-scale",category,mass_errors);
	   }

     }
   }

}


void LoopAll::myReducePhotonAnalysis(Util * ut, int jentry) {

  if(PADEBUG) 
    cout<<"myReducePhotonAnalysis START"<<endl;

  //count all events
  countersred[0]++;

  if(outputFile) {
    if(makeOutputTree) {
      
      //first selection and fill output tree
      if(!myFillReducedVarPhotonAnalysis(ut, jentry)) 
	return;
      
      //additional selection
      if(!mySelectEventRedPhotonAnalysis(ut, jentry)) 
	return;

      countersred[1]++;

      outputEvents++;
      if(PADEBUG) 
	cout<<"before fill"<<endl;

      outputTree->Fill();
      if(PADEBUG) 
	cout<<"after fill"<<endl;

      if(outputEvents==100) {
	outputEvents=0;
	outputTree->Write(0,TObject::kWriteDelete);
      }
    }
  }

  if(PADEBUG) 
    cout<<"myReducePhotonAnalysis END"<<endl;
}


void LoopAll::myGetBranchPhotonAnalysis() {
  b_event = fChain->GetBranch("event");
  b_lumis = fChain->GetBranch("lumis");
  b_run = fChain->GetBranch("run");
  b_pho_n = fChain->GetBranch("pho_n");
  b_pho_p4 = fChain->GetBranch("pho_p4");
  b_pho_r9 = fChain->GetBranch("pho_r9");
  b_pho_calopos = fChain->GetBranch("pho_calopos");
  b_pho_hoe = fChain->GetBranch("pho_hoe");
  b_pho_sieie = fChain->GetBranch("pho_sieie");
  b_pho_ecalsumetconedr03 = fChain->GetBranch("pho_ecalsumetconedr03");
  b_pho_ecalsumetconedr04 = fChain->GetBranch("pho_ecalsumetconedr04");
  b_pho_hcalsumetconedr03 = fChain->GetBranch("pho_hcalsumetconedr03");
  b_pho_hcalsumetconedr04 = fChain->GetBranch("pho_hcalsumetconedr04");
  b_pho_trksumptsolidconedr03 = fChain->GetBranch("pho_trksumptsolidconedr03");
  b_pho_trksumpthollowconedr04 = fChain->GetBranch("pho_trksumpthollowconedr04");
  b_pho_isEB = fChain->GetBranch("pho_isEB");
  b_pho_isEE = fChain->GetBranch("pho_isEE");
  b_pho_haspixseed = fChain->GetBranch("pho_haspixseed");

  b_gen_n = fChain->GetBranch("gp_n");
  b_gen_p4 = fChain->GetBranch("gp_p4");
  b_gen_status = fChain->GetBranch("gp_status");
  b_gen_pdgid = fChain->GetBranch("gp_pdgid");
  b_gen_mother = fChain->GetBranch("gp_mother");
}


void LoopAll::mySetBranchAddressRedPhotonAnalysis() {
  
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
  fChain->SetBranchAddress("run", &run, &b_run);
  fChain->SetBranchAddress("pho_n", &pho_n, &b_pho_n);
  fChain->SetBranchAddress("pho_p4", &pho_p4, &b_pho_p4);
  fChain->SetBranchAddress("pho_r9", &pho_r9, &b_pho_r9);
  fChain->SetBranchAddress("pho_calopos", &pho_calopos, &b_pho_calopos);
  fChain->SetBranchAddress("pho_hoe", &pho_hoe, &b_pho_hoe);
  fChain->SetBranchAddress("pho_sieie", &pho_sieie, &b_pho_sieie);
  fChain->SetBranchAddress("pho_ecalsumetconedr03", &pho_ecalsumetconedr03, &b_pho_ecalsumetconedr03);
  fChain->SetBranchAddress("pho_ecalsumetconedr04", &pho_ecalsumetconedr04, &b_pho_ecalsumetconedr04);
  fChain->SetBranchAddress("pho_hcalsumetconedr03", &pho_hcalsumetconedr03, &b_pho_hcalsumetconedr03);
  fChain->SetBranchAddress("pho_hcalsumetconedr04", &pho_hcalsumetconedr04, &b_pho_hcalsumetconedr04);
  fChain->SetBranchAddress("pho_trksumptsolidconedr03", &pho_trksumptsolidconedr03, &b_pho_trksumptsolidconedr03);
  fChain->SetBranchAddress("pho_trksumpthollowconedr04", &pho_trksumpthollowconedr04, &b_pho_trksumpthollowconedr04);
  fChain->SetBranchAddress("pho_isEB", &pho_isEB, &b_pho_isEB);
  fChain->SetBranchAddress("pho_isEE", &pho_isEE, &b_pho_isEE);
  fChain->SetBranchAddress("pho_haspixseed", &pho_haspixseed, &b_pho_haspixseed);

  fChain->SetBranchAddress("gp_n", &gen_n, &b_gen_n);
  fChain->SetBranchAddress("gp_p4", &gen_p4, &b_gen_p4);
  fChain->SetBranchAddress("gp_status", &gen_status, &b_gen_status);
  fChain->SetBranchAddress("gp_pdgid", &gen_pdgid, &b_gen_pdgid);
  fChain->SetBranchAddress("gp_mother", &gen_mother, &b_gen_mother);
}


int LoopAll::myFillReducedVarPhotonAnalysis(Util * ut, int jentry) {
  if(PADEBUG) 
    cout<<"myFillReduceVar START"<<endl;
  
   return 1;

  if(PADEBUG) 
    cout<<"myFillReduceVar END"<<endl;

}

int LoopAll::mySelectEventRedPhotonAnalysis(Util * ut, int jentry) {
 
 
  // preselection at the end
  int selectevent=0;

  b_pho_n->GetEntry(jentry);

  if (pho_n > 1 )//&&  (run >= 163217)) 
    selectevent = 1;
  else
    selectevent = 0;

  return selectevent;

}

