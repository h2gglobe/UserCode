#define PADEBUG 0

void LoopAll::TermRealPhotonAnalysis(int typerun) 
{
   if (typerun==3){	
      //rooContainer->FitToData("cms-data_model","Data_mass");
//      rooContainer->FitToData("rbw"	      ,"Signal_mass");
      rooContainer->FitToData("exp"	      ,"Background_mass");
//      rooContainer->FitToSystematicSet("rbw","Signal_mass","e-scale");
      rooContainer->FitToSystematicSet("exp","Background_mass","e-scale");
   }

}

void LoopAll::InitRealPhotonAnalysis(int typerun) {
  if(PADEBUG) 
    cout << "InitRealPhotonAnalysis START"<<endl;


  if (typerun == 2 || typerun == 1){

  }
  
  if (typerun == 3) {  
     // initialize the RooContainer with number of categories

     rooContainer = new RooContainer(1);

     //RooFitting type
     rooContainer->AddRealVar("Signal_mass",90.,200.);
     rooContainer->AddRealVar("Background_mass",90.,200.);
     rooContainer->AddRealVar("Data_mass",90.,200.);

     // Signal Parameters
     rooContainer->AddRealVar("width",1.,0.01,5.);
     rooContainer->AddRealVar("mean",130.,100.,150.);
     rooContainer->AddRealVar("sigma",1.,0.01,5.);

     // Background Parameters
     rooContainer->AddRealVar("c",+3.0,-5.,+5);
     rooContainer->AddRealVar("p1",-0.03,-1.,-0.001);
     rooContainer->AddRealVar("p2",+0.01,-1.,+1);

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
     //sig_pars[3] = "sigma";
     rooContainer->AddGenericPdf("rbw"
	,"1./((@0*@0-@1*@1)*(@0*@0-@1*@1)+@2*@2*@1*@1)"
	,sig_pars,0,10.);
    
     //std::vector<std::string> tail_pars(2,"p");	 
     //tail_pars[0] = "Signal_mass";
     //tail_pars[1] = "alpha";
     //rooContainer->AddGenericPdf("s-exp"
//	,"expo"
//	,tail_pars,1,1.);

     //std::vector<std::string> s_pdfs(2,"t");
     //s_pdfs[0] = "bw";
     //s_pdfs[1] = "s-exp";
  
     //rooContainer->ComposePdf("signal-model","bw+s-exp",s_pdfs);
     // -------------------------------------//

     // Set Up The Background Pdfs
     std::vector<std::string> bkg_pars(4,"p");	 
     bkg_pars[0] = "Background_mass";
     bkg_pars[1] = "c";
     bkg_pars[2] = "p1";
     bkg_pars[3] = "p2";
     rooContainer->AddGenericPdf("exp"
	,"@1+@2*@0+@3*@0*@0",bkg_pars,0,10);
     // -------------------------------------//

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
	,data_sig,2,1.);
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
     rooContainer->CreateDataSet("Data_mass",25);
     rooContainer->CreateDataSet("Signal_mass",110);
     rooContainer->CreateDataSet("Background_mass",55);

     // Systematic Errors
     rooContainer->MakeSystematics("Signal_mass","e-scale");
     rooContainer->MakeSystematics("Background_mass","e-scale");
  }

  if(PADEBUG) 
    cout << "InitRealPhotonAnalysis END"<<endl;

}

bool LoopAll::myPhotonAnalysisRunSelection(int cur_type){
  
  if (cur_type != 0)    return true;


  bool passes_run_selection = false; 

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

  if (cur_type>0){

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
	 if (min_r9 < 0.93 && max_eta < 1.4442 ) category = 1;
	 if (min_r9 > 0.93 && max_eta < 1.4442 ) category = 2;
	 if (min_r9 < 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 3;
	 if (min_r9 > 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 4;

         // -------------------------------------------------------
         TLorentzVector Higgs = (*(leading.p4))
                                 +(*(nleading.p4));
           float mass = Higgs.M();
           float h_pt = Higgs.Pt();

             if (mass > 100. && mass < 150.){

	       FillHist("pho_pt",category,leading.p4->Pt());
	       FillHist("pho_pt",category,nleading.p4->Pt());
               best_mass = mass;
 	       best_pt   = h_pt;

            }
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
	 if (min_r9 < 0.93 && max_eta < 1.4442 ) category = 1;
	 if (min_r9 > 0.93 && max_eta < 1.4442 ) category = 2;
	 if (min_r9 < 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 3;
	 if (min_r9 > 0.93 && max_eta > 1.566 && max_eta < 2.5) category = 4;

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
	     rooContainer->InputDataPoint("Data_mass",0,mass);
	     //rooContainer->InputDataPoint("Data_mass",category,mass);

           } else if (cur_type < 0){ // Signal
	     rooContainer->InputDataPoint("Signal_mass",0,mass,weight);
	     //rooContainer->InputDataPoint("Signal_mass",category,mass,weight);

             rooContainer->InputSystematicSet("Signal_mass","e-scale",0,mass_errors);
             //rooContainer->InputSystematicSet("Signal_mass","e-scale",category,mass_errors);

	   } else if (cur_type > 0){ // Background
	     rooContainer->InputDataPoint("Background_mass",0,mass,weight);
	     //rooContainer->InputDataPoint("Background_mass",category,mass,weight);

             rooContainer->InputSystematicSet("Background_mass","e-scale",0,mass_errors);
             //rooContainer->InputSystematicSet("Background_mass","e-scale",category,mass_errors);
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

  if (pho_n > 1)
    selectevent = 1;
  else
    selectevent = 0;

  return selectevent;

}

