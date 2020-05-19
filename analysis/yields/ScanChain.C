#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH2F.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "Math/VectorUtil.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TMVA/Reader.h"
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include "../misc/class_files/v8.02/SS.h"
#include "../../common/CORE/Tools/dorky/dorky.h"
#include "../../common/CORE/Tools/utils.h"
#include "../misc/common_utils.h"
#include "../misc/tqdm.h"

#define MZ 91.2

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>> Vec4;
// using namespace tas;
// float lumiAG = getLumi();
bool STOP_REQUESTED = false;
// float lumiAG = 36.3;
struct HistCol2D {
    map<string, TH2D> in;
    HistCol2D(vector<string> regions, const string& name, int nbinsx, float lowx, float highx, int nbinsy, float lowy, float highy, vector<HistCol2D*>* registry=nullptr) {
        for (string region : regions) {
            string base_name = region + "_" + name;
            string base_title = region + " " + name;
            in.emplace(region, TH2D((base_name + "_in").c_str(), (base_title + " in").c_str(), nbinsx, lowx, highx, nbinsy, lowy, highy));
        }
        if (registry != nullptr) registry->push_back(this);
    }
    void Fill(const string& region, int id1, int id2, float valx, float valy, float weight) { in[region].Fill(valx, valy, weight); }
    void Write() { for (auto p : in) p.second.Write(); }
};

struct HistCol {
    map<string, TH1D> in;
    map<string, TH1D> ee;
    map<string, TH1D> em;
    map<string, TH1D> mm;

    HistCol(vector<string> regions, const string& name, int nbins, const float* bins, vector<HistCol*>* registry=nullptr) {
        for (string region : regions) {
            string base_name = region + "_" + name;
            string base_title = region + " " + name;
            in.emplace(region, TH1D((base_name + "_in").c_str(), (base_title + " in").c_str(), nbins, bins));
            ee.emplace(region, TH1D((base_name + "_ee").c_str(), (base_title + " ee").c_str(), nbins, bins));
            em.emplace(region, TH1D((base_name + "_em").c_str(), (base_title + " em").c_str(), nbins, bins));
            mm.emplace(region, TH1D((base_name + "_mm").c_str(), (base_title + " mm").c_str(), nbins, bins));

        }
        if (registry != nullptr)
            registry->push_back(this);
    }

    HistCol(vector<string> regions, const string& name, int nbins, float low, float high, vector<HistCol*>* registry=nullptr) {
        for (string region : regions) {
            string base_name = region + "_" + name;
            string base_title = region + " " + name;
            in.emplace(region, TH1D((base_name + "_in").c_str(), (base_title + " in").c_str(), nbins, low, high));
            ee.emplace(region, TH1D((base_name + "_ee").c_str(), (base_title + " ee").c_str(), nbins, low, high));
            em.emplace(region, TH1D((base_name + "_em").c_str(), (base_title + " em").c_str(), nbins, low, high));
            mm.emplace(region, TH1D((base_name + "_mm").c_str(), (base_title + " mm").c_str(), nbins, low, high));
        }
        if (registry != nullptr)
            registry->push_back(this);
    }

    void Fill(const string& region, int id1, int id2, float val, float weight) {
        in[region].Fill(val, weight);

        if (abs(id1) == 11 and abs(id2) == 11) {
            ee[region].Fill(val, weight);
        } else if (abs(id1) == 13 and abs(id2) == 13) {
            mm[region].Fill(val, weight);
        } else if ((abs(id1) == 11 and abs(id2) == 13) or
                (abs(id1) == 13 and abs(id2) == 11)) {
            em[region].Fill(val, weight);
        } else {
            cout << "These ids are garbage: (" << id1 << ", " << id2 << ")\n";
        }    

    }

    void Write() {
        for (auto p : in) p.second.Write();
        for (auto p : ee) p.second.Write();
        for (auto p : em) p.second.Write();
        for (auto p : mm) p.second.Write();

    }
};

float calcDeltaPhi(float phi1, float phi2){
    float dPhi = phi1 - phi2;
    while (dPhi  >  TMath::Pi()) dPhi -= 2*TMath::Pi();
    while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
    return fabs(dPhi);
}

float calcDeltaR(float eta1, float phi1, float eta2, float phi2){
    return TMath::Sqrt(pow(calcDeltaPhi(phi1, phi2),2)+pow((eta1-eta2),2));

}

float calcMT(float pt1, float phi1, float pt2, float phi2){
    return sqrt( 2 * pt1 * pt2 * ( 1 - cos( phi1 - phi2 ) ) );
}


float minDR(const Vec4 & lep , const vector<Vec4> & jet)
{
    int size = (int)jet.size();
    //cout<<"jet size "<<size<<endl;
    float mindr = 999;  
    if(size)for(int i=0;i<size;i++){
        float dr = calcDeltaR(lep.eta(), lep.phi(), jet[i].eta(), jet[i].phi());
        if(dr<mindr) mindr = dr;
    }
    return mindr;
}

float PtMaxEta(const vector<Vec4> & jet)
{
    int size = (int)jet.size();
    //cout<<"jet size "<<size<<endl;
    float maxeta=0, eta=0;
    int index=0;
    for(int i=0;i<size;i++){
        eta = jet[i].eta();
        if(eta>maxeta){
            maxeta = eta;
            index = i;
        }
    }
    return jet[index].pt();
}

struct lepton
{
    Vec4 v;
    int id;
    bool isgood;
    float miniiso;
    float dxy;
    float dz;   

}; 

float MOSSF(vector<lepton> lep)
{
    int size = (int)lep.size();
    float mossf = 999999;
    float diff = 999999;
    if(size>=2)for(int i=0;i<size;i++){
        for(int j=i+1;j<size;j++){	
            if(lep[i].id*lep[j].id==-169||lep[i].id*lep[j].id==-121){
                float mass = (lep[i].v+lep[j].v).M();	  
                if(TMath::Abs(mass-MZ)<diff)
                {
                    mossf = mass;
                    diff = TMath::Abs(mass-MZ);
                }
            }	  
        }
    }

    return mossf;
}


float isr_reweight(bool useIsrWeight, int year, int nisrmatch, int sample) {
    if (!useIsrWeight) return 1;
    if (ss::is_real_data()) return 1;
    return isrWeight(year, nisrmatch, sample); // 10 is ttbar
}

float nb_reweight(int nbtags) {
    if (ss::is_real_data()) return 1;
    std::vector<float> weights = { 1.06, 0.96, 0.99, 1.29, 1.31 };
    // std::vector<float> weights = { 1.00, 1.00, 1.00, 1.29, 1.00 }; // FIXME only reweighting nb==3
    return weights[min(nbtags,4)];
}


int ScanChain(TChain *ch, TString options="", TString outputdir="outputs"){

    signal(SIGINT, [](int){
            cout << "SIGINT Caught, stopping after current event" << endl;
            STOP_REQUESTED=true;
            });
    STOP_REQUESTED=false;
    bool doFakes = options.Contains("doFakes");
    bool doFlips = options.Contains("doFlips");
    bool useNonIsoTriggers = options.Contains("useNonIsoTriggers");
    bool quiet = options.Contains("quiet");
    bool minPtFake18 = options.Contains("minPtFake18");
    bool new2016FRBins = options.Contains("new2016FRBins");
    bool noISRWeights = options.Contains("noISRWeights");
    bool isData = options.Contains("doData");

    //FG what is this for???
    //ana_t analysis = FTANA;
    //if (doSS) {
    ana_t analysis = SSANA;
    //}

    TString proc(ch->GetTitle());

    //cout<<"proc "<<proc<<endl;
    // bool useIsrWeight = proc.Contains("tt_");
    // FG: Don't think these are needed
    bool useIsrWeight = proc.Contains("tt_") 
        or proc.Contains("ttw") 
        or proc.Contains("ttz") 
        or proc.Contains("tth");

    bool doScaleUnc = proc.Contains("tt_");
    if (doScaleUnc && !quiet)
        cout << "Calculating scale uncertainty for process " << proc << endl;

    // this isn't needed
    bool useTTBB = proc.Contains("tt_") 
        or proc.Contains("ttw") 
        or proc.Contains("ttz") 
        or proc.Contains("tth");

    useTTBB = false;
    if (useTTBB and !quiet)
        cout << "Applying TTBB scale factors for process " << proc << endl;

    //cout<<"using isr "<<useIsrWeight<<" TTBB "<<useTTBB<<endl;  
    if (noISRWeights) 
        useIsrWeight = false;
    if (useIsrWeight and !quiet)
        cout << "Applying ISR weights to process " << proc << endl;

    // We may have derived the fake rate map throwing away leptons with pT<18 (e.g., 2017), so
    // we need to apply this cut here to be consistent
    //float min_pt_fake = minPtFake18 ? 18. : -1;
    float min_pt_fake = minPtFake18 ? 18. : -1;

    int year;
    float lumiAG = 0.;
    bool is2016(false), is2017(false), is2018(false);
    if (options.Contains("Data2016")) {
        lumiAG = getLumi(2016);
        year = 2016;
        is2016 = true;
    } else if (options.Contains("Data2017")) {
        lumiAG = useNonIsoTriggers ? 36.529: getLumi(2017);
        year = 2017;
        is2017 = true;
    } else if (options.Contains("Data2018")) {
        lumiAG = getLumi(2018);
        year = 2018;
        is2018 = true;
    } else {
        cout << "Need to specify year!\n";
        return -1;
    }
    // Clear already-seen list
    duplicate_removal::clear_list();

    // Used to determine which "era" a MC event is in
    TRandom *tr1 = new TRandom();

    if (!quiet) cout << "Working on " << proc << endl;

    vector<string> regions =
    {			    
        "ssbr"
    };

    vector<HistCol*> registry;
    vector<HistCol2D*> registry2D;
    HistCol h_met           (regions, "met"             , 20, 0   , 500 ,  &registry);
    HistCol h_nleps         (regions, "nleps"           , 5, -0.5 , 4.5 ,  &registry);
    HistCol h_type          (regions, "type"            , 8 , -0.5, 7.5 ,  &registry);

    // 2D hists
    // Declare a bunch of event variables to be filled below in the loop
    int nleps;
    int type;
    float met;
    float weight = 1.;

    //add list of trees
    float tree_weight = -1;
    float tree_met = -1;
    int tree_type = -1;
    int tree_nleps = -1;

    TFile *  f1 = new TFile(Form("%s/histos_%s.root", outputdir.Data(), ch->GetTitle()), "RECREATE");
    f1->cd();
    TTree* out_tree = new TTree("t","fortraining");

    out_tree->Branch("weight", &tree_weight);
    out_tree->Branch("nleps", &tree_nleps);
    out_tree->Branch("type", &tree_type);
    out_tree->Branch("met", &tree_met);

    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);

    tqdm bar;
    //  bar.set_theme_braille();

    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();
    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        if (STOP_REQUESTED) break;
        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("t");
        samesign.Init(tree);

        TString filename(currentFile->GetTitle());

        auto tokens = filename.Tokenize("/");
        auto basename = ((TObjString*)(tokens->At(tokens->GetEntries()-1)))->String().Data();
        bar.set_label(basename);

        for(unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {
            if (STOP_REQUESTED) break;
            samesign.GetEntry(event);
            nEventsTotal++;            
            if (!quiet) bar.progress(nEventsTotal, nEventsChain);

            weight = ss::is_real_data() ? 1 : ss::scale1fb()*lumiAG;

            // Simple cuts first to speed things up
            int lep1id = ss::lep1_id();
            int lep2id = ss::lep2_id();

            bool pass_trig =  ss::fired_trigger_ss();
            if (!pass_trig) continue;
            if (!ss::passes_met_filters()) continue;
            //
            // Reject duplicates
            if (ss::is_real_data()) {
                duplicate_removal::DorkyEventIdentifier id(ss::run(), ss::event(), ss::lumi());
                if (duplicate_removal::is_duplicate(id)) continue;
            }

            // Save a bunch of event info for quick reference later
            met = ss::met();

            int lep1good = ss::lep1_passes_id();
            int lep2good = ss::lep2_passes_id();
            int lep3good = ss::lep3_passes_id();
            nleps = (lep3good) ? ((ss::lep4_passes_id() and (ss::lep4_p4().pt() > (abs(ss::lep4_id())==11 ? 15 : 10))) ? 4 : 3) : 2;      
            vector<lepton> lep;
            lepton temp;
            temp.v = ss::lep1_p4(); temp.id = ss::lep1_id(); 
            lep.push_back(temp);
            temp.v = ss::lep2_p4(); temp.id = ss::lep2_id();
            lep.push_back(temp);

            /* hyp_class
             * 1: SS, loose-loose
             * 2: SS, tight-loose (or loose-tight)
             * 3: SS, tight-tight
             * 4: OS, tight-tight
             * 5: SS, inSituFR
             * 6: SS, tight-tight and fails Z-veto (lies! hyp_class==6 != lep1good and lep2good)
             */
            int hyp_class = ss::hyp_class();      
            type = hyp_class;
            float ht = ss::ht();
            int nisrmatch = ss::nisrMatch();      	    
            float lep1ccpt = ss::lep1_coneCorrPt();
            float lep2ccpt = ss::lep2_coneCorrPt();
            float lep3ccpt = ss::lep3_coneCorrPt();
            float lep1eta = ss::lep1_p4().eta();
            float lep2eta = ss::lep2_p4().eta();
            float lep3eta = ss::lep3_p4().eta();
            float lep1pt = ss::lep1_p4().pt();
            float lep2pt = ss::lep2_p4().pt();
            float lep3pt = ss::lep3_p4().pt();
            int lep3id = ss::lep3_id();
            if (!ss::is_real_data()) {
                weight *= getTruePUw(year, ss::trueNumInt()[0]);
                if (lep1good) weight *= leptonScaleFactor(year, lep1id, lep1ccpt, lep1eta, ht, FTANA);
                if (lep2good) weight *= leptonScaleFactor(year, lep2id, lep2ccpt, lep2eta, ht, FTANA);
                if (lep3good) weight *= leptonScaleFactor(year, lep3id, lep3ccpt, lep3eta, ht, FTANA);

                if (!lep3good) {
                    weight *= triggerScaleFactor(year, lep1id, lep2id, lep1pt, lep2pt, lep1eta, lep2eta, ht, analysis, 0);
                }
                weight *= ss::weight_btagsf();
                if (proc.Contains("ttw"))  weight *= isrWeight(year, nisrmatch, 1);
                if (proc.Contains("ttz"))  weight *= isrWeight(year, nisrmatch, 2);
                if (proc.Contains("tt_"))  weight *= isrWeight(year, nisrmatch, 10);

                if (year == 2016) weight *= ss::prefire2016_sf();
                if (year == 2017) weight *= ss::prefire2017_sf();		
                if (useTTBB and (ss::extragenb() >= 2)) weight *= 1.7;	
                weight *=ss::decayWSF();	
            }

            bool class6Fake = false;
            if (doFakes) {
                if (hyp_class == 6) {
                    bool lep1_lowpt_veto = lep1pt < (abs(lep1id) == 11 ? 15 : 10);
                    bool lep2_lowpt_veto = lep2pt < (abs(lep2id) == 11 ? 15 : 10);
                    bool lep3_lowpt_veto = lep3pt < (abs(lep3id) == 11 ? 15 : 10);
                    int nfakes = 0;
                    if (ss::lep3_fo() and !ss::lep3_tight() and !lep3_lowpt_veto and lep1good and lep2good && lep3pt>min_pt_fake) {  // lep3 fake
                        float fr = fakeRate(year, lep3id, lep3ccpt, lep3eta, ht, analysis, new2016FRBins, !minPtFake18);
                        class6Fake = true;
                        nfakes++;
                        weight *= fr / (1-fr);
                    }
                    if (ss::lep2_fo() and !ss::lep2_tight() and !lep2_lowpt_veto and lep1good and lep3good && lep2pt>min_pt_fake) {  // lep2 fake
                        float fr = fakeRate(year, lep2id, lep2ccpt, lep2eta, ht, analysis, new2016FRBins, !minPtFake18);
                        class6Fake = true;
                        nfakes++;
                        weight *= fr / (1-fr);
                    }
                    if (ss::lep1_fo() and !ss::lep1_tight() and !lep1_lowpt_veto and lep2good and lep3good && lep1pt>min_pt_fake) {  // lep1 fake
                        float fr = fakeRate(year, lep1id, lep1ccpt, lep1eta, ht, analysis, new2016FRBins, !minPtFake18);
                        class6Fake = true;
                        nfakes++;
                        weight *= fr / (1-fr);
                    }
                    if (!class6Fake) {
                        continue; // No fakes!
                    }
                    if (nfakes == 2) weight *= -1;
                } else if (hyp_class == 1 or hyp_class == 2) {
                    bool foundGoodLoose = false;
                    if (ss::lep1_passes_id()==0 && lep1pt>min_pt_fake) {
                        float fr = fakeRate(year, lep1id, lep1ccpt, lep1eta, ht, analysis, new2016FRBins, !minPtFake18);
                        weight *= fr/(1.-fr);
                        foundGoodLoose = true;
                    }
                    if (ss::lep2_passes_id()==0 && lep2pt>min_pt_fake) {
                        float fr = fakeRate(year, lep2id, lep2ccpt, lep2eta, ht, analysis, new2016FRBins, !minPtFake18);
                        weight *= fr/(1.-fr);
                        foundGoodLoose = true;
                    }
                    if (!foundGoodLoose)
                        continue;
                    // subtract double FO (why is this?)
                    if (hyp_class == 1 && lep1pt>min_pt_fake && lep2pt>min_pt_fake) weight *= -1.;
                    // subtract prompt MC
                    if (isData and !ss::is_real_data()) weight *= -1.;
                    hyp_class = 3; // we've faked a SS Tight-Tight with a SS LL or SS TL
                    // Basically just update this so it gets put in the SR	  	  
                } else {
                    continue; // Not a fakeing hyp_class
                }
            }      

            //flips
            if (doFlips) {
                if (hyp_class == 4) hyp_class = 3; // we've flipped an OS to a SS
                // else if (hyp_class == 6) class6Fake = true;
                else continue;
                float flipFact = 0.;
                if (abs(lep1id) == 11) {
                    float flr = flipRate(year, lep1pt, lep1eta, FTANA);
                    flipFact += (flr/(1-flr));
                }
                if (abs(lep2id) == 11) {
                    float flr = flipRate(year, lep2pt, lep2eta, FTANA);
                    flipFact += (flr/(1-flr));
                }
                weight *= flipFact;
                if (weight == 0.0) continue; // just quit if there are no flips.
            }

            auto getmll = [](const Vec4& p1, const Vec4& p2, float ccpt1=-1, float ccpt2=-1) {
                /* Calculate dilepton mass with optional rescaling based on cone-corrected lepton pt */
                if (ccpt1 == -1) return (p1 + p2).M();
                else             return (p1*ccpt1/p1.pt() + p2*ccpt2/p2.pt()).M();
            };
            float m12 = getmll(ss::lep1_p4(), ss::lep2_p4(), lep1ccpt, lep2ccpt);
            float m13 = getmll(ss::lep1_p4(), ss::lep3_p4(), lep1ccpt, lep3ccpt);
            float m23 = getmll(ss::lep2_p4(), ss::lep3_p4(), lep2ccpt, lep3ccpt);
            float m3l = (ss::lep1_p4()+ss::lep2_p4()+ss::lep3_p4()).M();

            auto z_cand = [](int id1, int id2, float mll) {
                return abs(id1) == abs(id2) and  // Same flavor
                    id1*id2<0 and             // Opposite sign
                    fabs(mll - MZ) < 15;     // Z-mass window
            };
            bool zcand12 = z_cand(lep1id, lep2id, m12);
            bool zcand13 = z_cand(lep1id, lep3id, m13);
            bool zcand23 = z_cand(lep2id, lep3id, m23);
            float mllos = fabs(m13 - MZ) < fabs(m23 - MZ) ? m13 : m23;	    

            auto fill_region = [&](const string& region, float weight) {
                if (std::find(regions.begin(), regions.end(), region) == regions.end()) return;	      
                // Fill all observables for a region
                auto do_fill = [region, lep1id, lep2id, weight](HistCol& h, float val, float extraweight=1.) {
                    h.Fill(region, lep1id, lep2id, val, weight*extraweight);
                };
                auto do_fill2D = [region, lep1id, lep2id, weight](HistCol2D& h, float valx, float valy) {
                    h.Fill(region, lep1id, lep2id, valx, valy, weight);
                };
                do_fill(h_met, met);
                do_fill(h_nleps, nleps);

                int looseleg = -1;
                if (hyp_class == 2) {
                    looseleg = (lep1good ? 2 : 1);
                }
                int type = ss::hyp_type();
                do_fill(h_type,   type>1 ? type-1 : type);
            };

            bool MossfOnZ = fabs(MOSSF(lep)-MZ)<15.;
            if  ( nleps > 1 ){ 
                fill_region("ssbr", weight); 
                tree_met = met;
                tree_weight = weight;	
                tree_nleps=nleps;	
                tree_type=type;
                out_tree->Fill();	
            }
        }//event loop

        //delete tree;
        delete file;
    }//file loop

    f1->cd();
    for (HistCol* coll : registry) coll->Write();
    for (HistCol2D* coll : registry2D) coll->Write();
    out_tree->Write();
    f1->Close();
    if (!quiet) cout << "\n Done!" << endl;
    return 0;
}
