#include "Pythia8/Pythia.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <map>
#include <sstream>

using namespace Pythia8;

int main(int argc, char *argv[]) {

//	define a vector filled with command-line arguments
	std::vector<std::string> args(argv,argv+argc);

//	print command-line arguments that are passed to the program
	if (argc > 1) {
      		cout << "Argument #0 (filename) = " << argv[0] << endl;
      		cout << "Argument #1 (random number seed) = " << argv[1] << endl; 
      		cout << "Argument #2 (number of events to generate) = " << argv[2] << endl;
   	}

//	check that the first command-line argument is in fact of the expected variable type 
//	and within the range of values allowed by PYTHIA, i.e. between -1 and 900,000,000 
	istringstream run_counter(argv[1]); int X;
	if (!(run_counter >> X) || X<-1 || X>900000000) cerr << "Invalid integer: " << argv[1] << '\n';	
	if (X==0) cout << "Random number seed is set based on time... (the input value is " << argv[1] << ")." << '\n';	
	if (X==-1) cout << "Choosing default seed... (the input value is " << argv[1] << ")." << '\n';	

//	create a .cmnd file that contains all relevant run parameters; use the first command-line 
//	argument as random number seed and the second to specify the number of events to generate
	string ofilename = string("dijet_alpha1_#") + argv[1] + string(".cmnd");
	string RndmSeed = string("Random:seed = ") + argv[1];
	int nEvents = atoi(argv[2]);
	string Nevents = string("Main:numberOfEvents = ") + argv[2];
	ofstream ofileRunParameters(ofilename.c_str());
//	ofileRunParameters.open(ofilename.c_str());
	if (!ofileRunParameters) cerr << "Cannot open output file." << endl;
  	ofileRunParameters << "Beams:idA = 2212\n";
  	ofileRunParameters << "Beams:idB = 2212\n";
  	ofileRunParameters << "Beams:eCM = 13000\n";
  	ofileRunParameters << "SigmaProcess:alphaSvalue = 0.13\n";
  	ofileRunParameters << "SigmaProcess:alphaSorder = 1\n";
  	ofileRunParameters << "HardQCD:all = on\n";
  	ofileRunParameters << "Random:setSeed = on\n";
  	ofileRunParameters << RndmSeed.c_str() << endl;
  	ofileRunParameters << "PhaseSpace:pTHatMin = 380\n";
	ofileRunParameters << "PhaseSpace:mHatMin = 2400\n";
  	ofileRunParameters << Nevents.c_str() << endl;
  	ofileRunParameters << "Next:numberShowEvent = 10\n";
  	ofileRunParameters << "Next:numberCount = 1000\n";
  	ofileRunParameters.close();

	Pythia pythia; // declare Pythia object
	pythia.readFile(ofilename.c_str()); // import settings from the file created above
//	int nEvents = pythia.mode("Main:numberOfEvents");
	int nList = pythia.mode("Next:numberShowEvent");
	pythia.init();

//	define jet parameters 
	int Alg = -1;		// -1 = anti-kT, 0 = Cambridge/Aachen, 1 = kT
	double R = 0.4;		// radius of the jet cone
	double pTmin = 50.0;	// minimum jet pT
	double etaRange = 5.0;	// pseudorapidity range of the detector
	int nSel = 2;		// exclude neutrinos (and other invisibles) from study

//	set up SlowJet jet finder
	SlowJet slowjet(Alg,R,pTmin,etaRange,nSel);

//	create histogram for invariant mass distribution
	Hist invariantMass("invariant mass",1000,0.,6000.);	

//	create histograms for jet pT (transverse momentum)
	Hist leading_jet_pT("Leading jet pT",1000,0.,3000.);
	Hist subleading_jet_pT("Subleading jet pT",1000,0.,3000.);
	Hist subsubleading_jet_pT("Subsubleading jet pT",1000,0.,3000.);

//	create histograms for jet y (rapidity)
	Hist leading_jet_y("Leading jet y",1000,-5.,5.);
	Hist subleading_jet_y("Subleading jet y",1000,-5.,5.);
	Hist subsubleading_jet_y("Subsubleading jet y",1000,-5.,5.);

//	create chi histograms for specific mass intervals
	Hist chi_2500_2800("chi for 2.5 TeV < m < 2.8 TeV",1000,0.,30.);
	Hist chi_2800_3100("chi for 2.8 TeV < m < 3.1 TeV",1000,0.,30.);
	Hist chi_3100_3400("chi for 3.1 TeV < m < 3.4 TeV",1000,0.,30.);
	Hist chi_3400_4000("chi for 3.4 TeV < m < 4.0 TeV",1000,0.,30.);
	Hist chi_4000_4600("chi for 4.0 TeV < m < 4.6 TeV",1000,0.,30.);
	Hist chi_4600_5400("chi for 4.6 TeV < m < 5.4 TeV",1000,0.,30.);
	Hist chi_5400_13000("chi for 5.4 TeV < m",1000,0.,30.);

//	create histograms for yStar and yBoost
	Hist yStar("yStar",1000,0.,2.5);
	Hist yBoost("yBoost",1000,0.,2.5);

//	create histograms for theta, phi and eta
	Hist leading_jet_phi("leading jet phi",1000,-3.5,3.5);
	Hist subleading_jet_phi("subleading jet phi",1000,-3.5,3.5);
	Hist subsubleading_jet_phi("subsubleading jet phi",1000,-3.5,3.5);

// create variables to be defined in the event loop; with the
// exception of Chi and deltaR, these are used to apply selection cuts
	double invM,ystar,yboost,Chi;
	Vec4 p0,p1;
		
//	create counters to print additional run statistics
	bool at_least_two_jets;
	int if_loop_counter = 0;	
	int ystar_counter = 0;
	int yboost_counter = 0;
	int y_counter = 0;
	int chi_counter = 0;
	int invM_2500_counter = 0;
	int invM_2500_2800_counter = 0;
	int invM_2800_3100_counter = 0;
	int invM_3100_3400_counter = 0;
	int invM_3400_4000_counter = 0;
	int invM_4000_4600_counter = 0;
	int invM_4600_5400_counter = 0;
	int invM_5400_13000_counter = 0;

//	event loop
	for(int iEvent = 0; iEvent < nEvents; ++iEvent) {
		if(!pythia.next()) continue;
//	analyse properties of the SlowJet object
    slowjet.analyze(pythia.event);
  	if (iEvent < nList) slowjet.list();

		std::cout << "---------- Listing jets above threshold for event #" << iEvent+1 << ": "; 

		if (slowjet.sizeJet()>2 && slowjet.pT(0)>440 && slowjet.pT(1)>50) {		
			at_least_two_jets = true;
			if (at_least_two_jets == true) ++if_loop_counter;

			std::cout << slowjet.sizeJet() << " jets were found." << " ----------" << '\n';

			std::cout << "Leading jet pT = " << slowjet.pT(0) << '\n';
			std::cout << "Subleading jet pT = " << slowjet.pT(1) << '\n';
			std::cout << "Subsubleading jet pT = " << slowjet.pT(2) << '\n';

			std::cout << "Leading jet y = " << slowjet.y(0) << '\n';
			std::cout << "Subleading jet y = " << slowjet.y(1) << '\n';
			std::cout << "Subsubleading jet y = " << slowjet.y(2) << '\n';

			std::cout << "Leading jet phi = " << slowjet.phi(0) << '\n';
			std::cout << "Subleading jet phi = " << slowjet.phi(1) << '\n';
			std::cout << "Subsubleading jet phi = " << slowjet.phi(2) << '\n';

			ystar = abs(slowjet.y(0) - slowjet.y(1))/2;
			yboost = abs(slowjet.y(0) + slowjet.y(1))/2;
			Chi = exp(abs(slowjet.y(0) - slowjet.y(1)));
			p0 = slowjet.p(0);
			p1 = slowjet.p(1);
			invM = m(p0,p1);

			std::cout << "ystar = " << ystar << '\n';
			std::cout << "yboost = " << yboost << '\n';
			std::cout << "chi = " << Chi << '\n';
			std::cout << "invM = " << invM << '\n';
			
			if(Chi<30) ++chi_counter;
			if(ystar<1.7) ++ystar_counter;
			if(yboost<1.1) ++yboost_counter;
			if(ystar<1.7 && yboost<1.1) ++y_counter;
			if(ystar<1.7 && yboost<1.1 && 2500>invM) ++invM_2500_counter;
			if(ystar<1.7 && yboost<1.1 && 2500<invM && invM<2800) ++invM_2500_2800_counter;
			if(ystar<1.7 && yboost<1.1 && 2800<invM && invM<3100) ++invM_2800_3100_counter;
			if(ystar<1.7 && yboost<1.1 && 3100<invM && invM<3400) ++invM_3100_3400_counter;
			if(ystar<1.7 && yboost<1.1 && 3400<invM && invM<4000) ++invM_3400_4000_counter;
			if(ystar<1.7 && yboost<1.1 && 4000<invM && invM<4600) ++invM_4000_4600_counter;
			if(ystar<1.7 && yboost<1.1 && 4600<invM && invM<5400) ++invM_4600_5400_counter;
			if(ystar<1.7 && yboost<1.1 && 5400<invM) ++invM_5400_13000_counter;

//	fill histograms, apply angular selection cuts
			if(ystar<1.7 && yboost<1.1) {
				invariantMass.fill(invM);
				leading_jet_pT.fill(slowjet.pT(0));
				subleading_jet_pT.fill(slowjet.pT(1));
				subsubleading_jet_pT.fill(slowjet.pT(2));
    		leading_jet_y.fill(slowjet.y(0));
    		subleading_jet_y.fill(slowjet.y(1));
    		subsubleading_jet_y.fill(slowjet.y(2));
				yStar.fill(ystar);
				yBoost.fill(yboost);
				leading_jet_phi.fill(slowjet.phi(0));
				subleading_jet_phi.fill(slowjet.phi(1));
				subsubleading_jet_phi.fill(slowjet.phi(2));
			}

//	angular selection cuts + invariant mass cuts for chi histograms	
			if(ystar<1.7 && yboost<1.1 && 2500<invM && invM<2800) {
				chi_2500_2800.fill(Chi);
			}
			if(ystar<1.7 && yboost<1.1 && 2800<invM && invM<3100) {
				chi_2800_3100.fill(Chi);
			}
			if(ystar<1.7 && yboost<1.1 && 3100<invM && invM<3400) {
				chi_3100_3400.fill(Chi);
			}
			if(ystar<1.7 && yboost<1.1 && 3400<invM && invM<4000) {
				chi_3400_4000.fill(Chi);
			}
			if(ystar<1.7 && yboost<1.1 && 4000<invM && invM<4600) {
				chi_4000_4600.fill(Chi);
			}
			if(ystar<1.7 && yboost<1.1 && 4600<invM && invM<5400) {
				chi_4600_5400.fill(Chi);
			}
			if(ystar<1.7 && yboost<1.1 && 5400<invM) {
				chi_5400_13000.fill(Chi);
			}
		}

		else std::cout << "Threshold conditions for leading and/or subleading jet pT not met. Not starting algorithm." << " ----------" << '\n';
/*
		else if (slowjet.sizeJet() == 1)
		std::cout << "Only " << slowjet.sizeJet() << " jet was found. Not starting algorithm." << " ----------" << '\n';
		else if (slowjet.sizeJet() == 0)
		std::cout << "No jets were found about the specified threshold (" << pTmin << " GeV)." << " ----------" << '\n';
*/
	}
//	end of event loop

//	define ratios
	double multiple_jets_ratio = 1.0*if_loop_counter/nEvents;
	double chi_ratio = 1.0*chi_counter/if_loop_counter;
	double ystar_ratio = 1.0*ystar_counter/if_loop_counter;
	double yboost_ratio = 1.0*yboost_counter/if_loop_counter;
	double y_ratio = 1.0*y_counter/if_loop_counter;
	double chi_2500_ratio = 1.0*invM_2500_counter/y_counter;
	double chi_2500_2800_ratio = 1.0*invM_2500_2800_counter/y_counter;
	double chi_2800_3100_ratio = 1.0*invM_2800_3100_counter/y_counter;
	double chi_3100_3400_ratio = 1.0*invM_3100_3400_counter/y_counter;
	double chi_3400_4000_ratio = 1.0*invM_3400_4000_counter/y_counter;
	double chi_4000_4600_ratio = 1.0*invM_4000_4600_counter/y_counter;
	double chi_4600_5400_ratio = 1.0*invM_4600_5400_counter/y_counter;
	double chi_5400_13000_ratio = 1.0*invM_5400_13000_counter/y_counter;
	
//	print additional run statistics
	std::cout << "\n";
	std::cout << "---------- ADDITIONAL RUN STATISTICS ----------" << "\n\n";
	std::cout << if_loop_counter << " out of " << nEvents << " events " << "(" << 100*multiple_jets_ratio << " %)" << " contained at least two jets above the threshold (" << pTmin << " GeV).\n"; 
	std::cout << chi_counter << " out of the " << if_loop_counter << "/" << nEvents << " events " << "(" << 100*chi_ratio << " %)" << " with at least two jets yield a chi value < 30.\n";
	std::cout << ystar_counter << " out of the " << if_loop_counter << "/" << nEvents << " events " << "(" << 100*ystar_ratio << " %)" << " with at least two jets yield a y* value < 1.7.\n";
	std::cout << yboost_counter << " out of the " << if_loop_counter << "/" << nEvents << " events " << "(" << 100*yboost_ratio << " %)" << " with at least two jets yield a yboost value < 1.1.\n";
	std::cout << y_counter << " out of the " << if_loop_counter << "/" << nEvents << " events " << "(" << 100*y_ratio << " %)" << " with at least two jets yield both a y* value < 1.7 and yboost value < 1.1.\n";
	std::cout << invM_2500_counter << " out of the " << y_counter << "/" << nEvents << " events " << "(" << 100*chi_2500_ratio << " %)" << " with at least two jets yield a invM value < 2500 after angular selection cuts have been applied.\n";
	std::cout << invM_2500_2800_counter << " out of the " << y_counter << "/" << nEvents << " events " << "(" << 100*chi_2500_2800_ratio << " %)" << " with at least two jets yield a invM value > 2500 and < 2800 after angular selection cuts have been applied.\n";
	std::cout << invM_2800_3100_counter << " out of the " << y_counter << "/" << nEvents << " events " << "(" << 100*chi_2800_3100_ratio << " %)" << " with at least two jets yield a invM value > 2800 and < 3100 after angular selection cuts have been applied.\n";
	std::cout << invM_3100_3400_counter << " out of the " << y_counter << "/" << nEvents << " events " << "(" << 100*chi_3100_3400_ratio << " %)" << " with at least two jets yield a invM value > 3100 and < 3400 after angular selection cuts have been applied.\n";
	std::cout << invM_3400_4000_counter << " out of the " << y_counter << "/" << nEvents << " events " << "(" << 100*chi_3400_4000_ratio << " %)" << " with at least two jets yield a invM value > 3400 and < 4000 after angular selection cuts have been applied.\n";
	std::cout << invM_4000_4600_counter << " out of the " << y_counter << "/" << nEvents << " events " << "(" << 100*chi_4000_4600_ratio << " %)" << " with at least two jets yield a invM value > 4000 and < 4600 after angular selection cuts have been applied.\n";
	std::cout << invM_4600_5400_counter << " out of the " << y_counter << "/" << nEvents << " events " << "(" << 100*chi_4600_5400_ratio << " %)" << " with at least two jets yield a invM value > 4600 and < 5400 after angular selection cuts have been applied.\n";
	std::cout << invM_5400_13000_counter << " out of the " << y_counter << "/" << nEvents << " events " << "(" << 100*chi_5400_13000_ratio << " %)" << " with at least two jets yield a invM value > 5400 after angular selection cuts have been applied.\n";

//	print statistics
	pythia.stat();

//	print histogram for the invariant mass distribution
	cout << invariantMass;
//	print histograms for leading/subleading jet pT and y
	cout << leading_jet_pT << subleading_jet_pT << subsubleading_jet_pT << leading_jet_y << subleading_jet_y << subsubleading_jet_y;
//	print histograms for yStar and yBoost
	cout << yStar << yBoost; 
//	print chi histograms for sepcific mass intervals
	cout << chi_2500_2800 << chi_2800_3100 << chi_3100_3400 << chi_3400_4000 << chi_4000_4600 << chi_4600_5400 << chi_5400_13000;
// print histograms for phi
	cout << leading_jet_phi << subleading_jet_phi << subsubleading_jet_phi;

//	create a map linking each PYTHIA histogram with a filename
	map<string,Hist> histograms;
	histograms["dijet_invMass_alpha1"] = invariantMass;
	histograms["dijet_leading_pT_alpha1"] = leading_jet_pT;
	histograms["dijet_subleading_pT_alpha1"] = subleading_jet_pT;
	histograms["dijet_subsubleading_pT_alpha1"] = subsubleading_jet_pT;
	histograms["dijet_leading_y_alpha1"] = leading_jet_y;
	histograms["dijet_subleading_y_alpha1"] = subleading_jet_y;
	histograms["dijet_subsubleading_y_alpha1"] = subsubleading_jet_y;
	histograms["dijet_yStar_alpha1"] = yStar;
	histograms["dijet_yBoost_alpha1"] = yBoost;
	histograms["dijet_chi_2500_2800_alpha1"] = chi_2500_2800;
	histograms["dijet_chi_2800_3100_alpha1"] = chi_2800_3100;
	histograms["dijet_chi_3100_3400_alpha1"] = chi_3100_3400;
	histograms["dijet_chi_3400_4000_alpha1"] = chi_3400_4000;
	histograms["dijet_chi_4000_4600_alpha1"] = chi_4000_4600;
	histograms["dijet_chi_4600_5400_alpha1"] = chi_4600_5400;
	histograms["dijet_chi_5400_13000_alpha1"] = chi_5400_13000;
	histograms["dijet_leading_phi_alpha1"] = leading_jet_phi;
	histograms["dijet_subleading_phi_alpha1"] = subleading_jet_phi;
	histograms["dijet_subsubleading_phi_alpha1"] = subsubleading_jet_phi;

//	iterate through the map elements (by default in alphabetical order) and write histogram data to .txt files
	for (map<string,Hist>::const_iterator itMap=histograms.begin(); itMap!= histograms.end(); ++itMap) {
		const char* filepath = "MCdata/Dijet/";
		string filename = filepath + itMap->first + string("_#") + argv[1] + string("_") + argv[2] + string(".txt");
//	const char *FileName = FileName.c_str();
		itMap->second.table(filename.c_str(),false,false);
	}

	return 0;
}
