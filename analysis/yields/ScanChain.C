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
#include <fstream>

#define MZ 91.2

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>> Vec4;
// using namespace tas;
// float lumiAG = getLumi();
bool STOP_REQUESTED = false;
// float lumiAG = 36.3;
//
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

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

struct Jet
{
    Jet () : v(LorentzVector(0,0,0,0)),islep(false),pt(0),eta(0),phi(0) {}
    Jet (LorentzVector v_) : v(v_),islep(false),pt(v.pt()),eta(v.eta()),phi(v.phi()) {}
    LorentzVector v;
    bool islep;
    float pt;
    float eta;
    float phi;
};

std::pair<int,int> getMaxMjjJets (std::vector<Jet> vj) {
    float mass = 0.;
    std::pair<int,int> ret = std::make_pair(-1,-1);
    for (unsigned int i = 0; i < vj.size(); i++) {
        for (unsigned int j =i+1; j < vj.size(); j++) {
            float tmp_mass = (vj.at(i).v+vj.at(j).v).mass();
            if (tmp_mass > mass) 
                ret = std::make_pair(i,j);
        }
    }
    return ret;
}

float max_mjj(std::vector<Jet> vj) {
    std::pair<int,int> indices = getMaxMjjJets(vj);    
    return (vj.at(indices.first).v+vj.at(indices.second).v).mass();
}

float max_mjj_deta(std::vector<Jet> vj) {
    std::pair<int,int> indices = getMaxMjjJets(vj);    
    return (vj.at(indices.first).eta-vj.at(indices.second).eta);
}

float max_mjj_dphi(std::vector<Jet> vj) {
    std::pair<int,int> indices = getMaxMjjJets(vj);    
    return calcDeltaPhi(vj.at(indices.first).phi,vj.at(indices.second).phi);
}

float max_mjj_dr(std::vector<Jet> vj) {
    std::pair<int,int> indices = getMaxMjjJets(vj);    
    float dphi = calcDeltaPhi(vj.at(indices.first).phi,vj.at(indices.second).phi);
    float deta = vj.at(indices.first).eta-vj.at(indices.second).eta;
    return TMath::Sqrt(dphi*dphi + deta*deta);
}

struct Lepton
{
    Lepton () : name(""),v(LorentzVector(0,0,0,0)),id(0),isgood(false),isfo(false),ccpt(0.),miniiso(999.),dxy(999.),dz(999.),pt(0),eta(0),phi(0) {}
    Lepton (string name_, LorentzVector v_, int id_, bool isgood_, bool isfo_, float ccpt_, float miniiso_, float dxy_, float dz_) : 
            name(name_),v(v_),id(id_),isgood(isgood_),isfo(isfo_),ccpt(ccpt_),miniiso(miniiso_),dxy(dxy_),dz(dz_) {pt = v.pt(); eta = v.eta(); phi = v.phi();}
    Lepton (const Lepton &lep) {name=lep.name; v=lep.v; id=lep.id; isgood=lep.isgood; isfo=lep.isfo; ccpt=lep.ccpt; miniiso=lep.miniiso; dxy=lep.dxy; dz=lep.dz; pt=v.pt(); eta=v.eta(); phi=v.phi();}
    string name;
    LorentzVector v;
    int id;
    bool isgood;
    bool isfo;
    float ccpt; 
    float miniiso;
    float dxy;
    float dz;   
    float pt;
    float eta;
    float phi;

    void print() {std::cout << this->name << ": pt = " << this->v.pt() << ", eta = " << this->v.eta() << ", is_good = " << this->isgood << ", isfo = " << this->isfo << std::endl;}
}; 

bool lep_comparator(const Lepton& lhs, const Lepton& rhs) {
    if (abs(lhs.id) > abs(rhs.id)) return true;
    else if (abs(lhs.id) < abs(rhs.id)) return false;
    else return (lhs.pt > rhs.pt);
}

std::vector<Jet> remove_overlaps(std::vector<Jet> vj, std::vector<Lepton> vl) {
    std::vector<Jet> ret;
    for (auto j : vj) {
        bool islep = false;
        for (auto l : vl) {
            if (ROOT::Math::VectorUtil::DeltaR(j.v,l.v) < 0.4) islep = true;
        }
        if (islep) continue;
        ret.push_back(j);
    }

    return ret;
}

struct Dilepton
{
    Dilepton (Lepton l1, Lepton l2);
    Lepton lep1;
    Lepton lep2;
    int type; 
    LorentzVector v;
    bool isgood;
    bool passpt;
    bool isss; 
};

Dilepton::Dilepton(Lepton l1, Lepton l2) : v(lep1.v+lep2.v)
{
    lep1 = (l1.v.pt() > l2.v.pt()) ? l1 : l2; 
    lep2 = (l1.v.pt() < l2.v.pt()) ? l1 : l2; 
    int q2 = abs(lep1.id*lep2.id);
    //
    // dilepton types
    // 0: ee
    // 1: em
    // 2: me
    // 3: mm
    if (q2 == 121) type = 0;
    else if (q2 == 169) type = 3;
    else if (abs(lep1.id) == 11) type = 1;
    else type = 2;

    isgood = (l1.isgood and l2.isgood);
    passpt = lep1.pt > 25 and lep2.pt > 20;
    isss = (lep1.id*lep2.id > 0);
}

vector<Dilepton> make_hyps(std::vector<Lepton> vl) 
{
    std::vector<Dilepton> ret;

    for (unsigned int i = 0; i < vl.size(); i++) {
        for (unsigned int j=i+1; j < vl.size(); j++) {
            ret.push_back(Dilepton(vl.at(i),vl.at(j)));
        }
    }

    return ret;
}

bool dilep_comparator(const Dilepton& lhs, const Dilepton& rhs) {
    if (lhs.type > rhs.type) return true;
    else if (lhs.type < rhs.type) return false;
    else return (lhs.v.pt() > rhs.v.pt());
}

std::vector<Dilepton> sort_dileps (std::vector<Dilepton> dl) {
    if (dl.size() < 2) return dl;
    sort(dl.begin(),dl.end(),&dilep_comparator);
    return dl; 
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

bool pass_cut(bool cut_val, int& cut_index,TH1D *h, float weight) {
    if (cut_val) cut_index++;
    h->Fill(cut_index,weight); 
    return cut_val;
}

void write_debug(std::ofstream& inf, int r, int l, int e, int c) {
    inf << Form("%d,%d,%d,%d",r,l,e,c) << endl;
    return;
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
    bool debug = options.Contains("debug");


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

    // open file to write debug
    std::ofstream dfile;
    std::cout << "debug: " << debug << ", outputdir  = " << outputdir.Data() << ", title = " << ch->GetTitle() << endl;
    if (debug and not dfile.is_open()) dfile.open(Form("%s/debug_%s.dat", outputdir.Data(), ch->GetTitle()), ios::out);

    vector<string> regions =
    {			    
        "br",
        "ssbr",
        "mlbr",
        "ssbr_mjj500",
        "mlbr_mjj500",
        "ssbr_mjj500_njgt2",
        "mlbr_mjj500_njgt2",
        "br_deta",
        "ssbr_deta",
        "mlbr_deta",
        "ssbr_mjj500_deta",
        "mlbr_mjj500_deta",
        "ssbr_mjj500_deta_njgt2",
        "mlbr_mjj500_deta_njgt2",
    };

    vector<HistCol*> registry;
    vector<HistCol2D*> registry2D;
    HistCol h_met           (regions, "met"             , 20, 0   , 400 ,  &registry);
    HistCol h_nleps         (regions, "nleps"           , 5, -0.5 , 4.5 ,  &registry);
    HistCol h_ss_type       (regions, "ss_type"         , 4 , -0.5, 3.5 ,  &registry);
    HistCol h_ml_type       (regions, "ml_type"         , 4 , -0.5, 3.5 ,  &registry);
    HistCol h_njets         (regions, "njets"           , 8, -0.5, 7.5  ,  &registry);
    HistCol h_ptl1          (regions, "ptl1"            , 40, 0, 400,  &registry);
    HistCol h_ptl2          (regions, "ptl2"            , 40, 0, 400,  &registry);
    HistCol h_ss_mll        (regions, "ss_mll"          , 40, 0, 400,  &registry);
    HistCol h_ml_mll        (regions, "ml_mll"          , 40, 0, 400,  &registry);
    HistCol h_mjj           (regions, "mjj"             , 20, 0, 1000,  &registry);
    HistCol h_mjj_deta      (regions, "mjj_deta"        , 100, -5, 5,  &registry);
    HistCol h_mjj_dphi      (regions, "mjj_dphi"        , 32, 0, 3.2,  &registry);
    HistCol h_mjj_dr        (regions, "mjj_dr"          , 60, 0, 6,  &registry);
    HistCol h_ptj1          (regions, "ptj1"            , 50, 0, 500,  &registry);
    HistCol h_ptj2          (regions, "ptj2"            , 50, 0, 500,  &registry);
    
    TH1D *h_cuts = new TH1D("cuts", "cuts", 12, -0.5, 11.5);

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

            // count which cut we're on
            int cut_index = 0;

            bool pass_trig =  ss::fired_trigger_ss();
            if (!pass_trig) continue;
            if (!ss::passes_met_filters()) continue;
            //
            // Reject duplicates
            if (ss::is_real_data()) {
                duplicate_removal::DorkyEventIdentifier id(ss::run(), ss::event(), ss::lumi());
                if (duplicate_removal::is_duplicate(id)) continue;
            }

            pass_cut(true,cut_index,h_cuts,1); 
            // cut_index incremented to 1

            // Save a bunch of event info for quick reference later
            met = ss::met();

            Lepton tmp_lep1("lep1", ss::lep1_p4(), ss::lep1_id(), ss::lep1_passes_id(), ss::lep1_fo(), ss::lep1_coneCorrPt(), ss::lep1_miniIso(), 99, 99);
            Lepton tmp_lep2("lep2", ss::lep2_p4(), ss::lep2_id(), ss::lep2_passes_id(), ss::lep2_fo(), ss::lep2_coneCorrPt(), ss::lep2_miniIso(), 99, 99);
            Lepton tmp_lep3("lep3", ss::lep3_p4(), ss::lep3_id(), ss::lep3_passes_id(), ss::lep3_fo(), ss::lep3_coneCorrPt(), 99, 99, 99);
            Lepton tmp_lep4("lep4", ss::lep4_p4(), ss::lep4_id(), ss::lep4_passes_id(), ss::lep4_fo(), ss::lep4_coneCorrPt(), 99, 99, 99);
            std::vector<Lepton> vlep = {tmp_lep1,tmp_lep2,tmp_lep3,tmp_lep4};
            std::vector<Dilepton> vdilep = make_hyps(vlep);
            std::vector<Dilepton> vdilep_good_ss;
            for (auto d : vdilep) {
                if (d.isgood and d.passpt and d.isss) vdilep_good_ss.push_back(d);
            }

            if (vdilep_good_ss.size() == 0) continue;

            nleps = 0;
            for (auto l : vlep) {
                if (l.isgood and l.pt > 20) nleps++;
            }

            std::vector<Dilepton> vdilep_good_ss_sorted = sort_dileps(vdilep_good_ss);
            Lepton lep1 = vdilep_good_ss_sorted.at(0).lep1;
            Lepton lep2 = vdilep_good_ss_sorted.at(0).lep2;

            std::vector<Lepton> vleps_good_not_hyp_sorted;
            for (auto l : vlep) {
                if (l.name == lep1.name or l.name == lep2.name) continue;
                if (not l.isgood) continue;
                if (l.pt < 20) continue;
                vleps_good_not_hyp_sorted.push_back(l); 
            }
            sort(vleps_good_not_hyp_sorted.begin(),vleps_good_not_hyp_sorted.end(),&lep_comparator);

            Lepton lep3; 
            Lepton lep4; 
            if (vleps_good_not_hyp_sorted.size() > 0) 
                lep3 = vleps_good_not_hyp_sorted.at(0);
            if (vleps_good_not_hyp_sorted.size() > 1) 
                lep4 = vleps_good_not_hyp_sorted.at(1);
            
            int hyp_class = ss::hyp_class();      
            type = hyp_class;

            if (!ss::is_real_data()) {
                weight *= getTruePUw(year, ss::trueNumInt()[0]);
                if (lep1.isgood) weight *= leptonScaleFactor(year, lep1.id, lep1.ccpt, lep1.eta, ss::ht(), FTANA);
                if (lep2.isgood) weight *= leptonScaleFactor(year, lep2.id, lep2.ccpt, lep2.eta, ss::ht(), FTANA);
                if (lep3.isgood) weight *= leptonScaleFactor(year, lep3.id, lep3.ccpt, lep3.eta, ss::ht(), FTANA);

                if (!lep3.isgood) {
                    weight *= triggerScaleFactor(year, lep1.id, lep2.id, lep1.pt, lep2.pt, lep1.eta, lep2.eta, ss::ht(), analysis, 0);
                }
                weight *= ss::weight_btagsf();
                if (proc.Contains("ttw"))  weight *= isrWeight(year, ss::nisrMatch(), 1);
                if (proc.Contains("ttz"))  weight *= isrWeight(year, ss::nisrMatch(), 2);
                if (proc.Contains("tt_"))  weight *= isrWeight(year, ss::nisrMatch(), 10);

                if (year == 2016) weight *= ss::prefire2016_sf();
                if (year == 2017) weight *= ss::prefire2017_sf();		
                //if (useTTBB and (ss::extragenb() >= 2)) weight *= 1.7;	
                weight *=ss::decayWSF();	
            }

            if (not pass_cut(vdilep_good_ss.size(),cut_index,h_cuts,weight)) {
               if (debug and dfile.is_open()) write_debug(dfile,ss::run(),ss::lumi(),ss::event(),cut_index); 
               continue;
            }
            // cut_index incremented to 2

            auto getmll = [](const Vec4& p1, const Vec4& p2, float ccpt1=-1, float ccpt2=-1) {
                /* Calculate dilepton mass with optional rescaling based on cone-corrected lepton pt */
                if (ccpt1 == -1) return (p1 + p2).M();
                else             return (p1*ccpt1/p1.pt() + p2*ccpt2/p2.pt()).M();
            };

            if (pass_cut(ss::madeExtraZ(),cut_index,h_cuts,weight)) {
               if (debug and dfile.is_open()) write_debug(dfile,ss::run(),ss::lumi(),ss::event(),cut_index);               
               continue;
            }
            // cut_index incremented to 3

            // store jets
            std::vector<Lepton> vlep_good;
            for (auto l : vlep) {
                if (l.pt > 20 and l.isgood) vlep_good.push_back(l);
            }

            std::vector<Jet> jets;
            for (auto j : ss::jets()) {
                if (j.pt() < 40.) continue;
                jets.push_back(Jet(j)); 
            }
            std::vector<Jet> jets_ro = remove_overlaps(jets,vlep_good);
            unsigned int njets = jets_ro.size();
            sort(jets_ro.begin(),jets_ro.end(),[](const Jet &lhs, const Jet &rhs) { return lhs.v.pt() > rhs.v.pt(); } );

            if (not pass_cut(njets>1,cut_index,h_cuts,weight)) {
               if (debug and dfile.is_open()) write_debug(dfile,ss::run(),ss::lumi(),ss::event(),cut_index); 
               continue; // cut_index incremented to 4 
            }

            if (pass_cut(ss::nbtags(),cut_index,h_cuts,weight)) {
               if (debug and dfile.is_open()) write_debug(dfile,ss::run(),ss::lumi(),ss::event(),cut_index); 
               continue; // cut_index incremented to 5 
            }

            float mjj_max = max_mjj(jets_ro);
            float mjj_max_deta = max_mjj_deta(jets_ro);
            float mjj_max_dphi = max_mjj_dphi(jets_ro);
            float mjj_max_dr   = max_mjj_dr(jets_ro);

            if (not pass_cut(mjj_max>120, cut_index,h_cuts,weight)) {
               if (debug and dfile.is_open()) write_debug(dfile,ss::run(),ss::lumi(),ss::event(),cut_index); 
               continue; // cut_index incremented to 4 
            }

            pass_cut(mjj_max>120, cut_index,h_cuts,weight);
            pass_cut(fabs(mjj_max_deta)>2.5, cut_index,h_cuts,weight);
            pass_cut(mjj_max>500, cut_index,h_cuts,weight);
            pass_cut(njets>2, cut_index,h_cuts,weight);
            pass_cut(njets>3, cut_index,h_cuts,weight);
            pass_cut(njets>4, cut_index,h_cuts,weight);
            pass_cut(njets>5, cut_index,h_cuts,weight);

            auto fill_region = [&](const string& region, float weight) {
                if (std::find(regions.begin(), regions.end(), region) == regions.end()) return;	      
                // Fill all observables for a region
                auto do_fill = [region, lep1, lep2, weight](HistCol& h, float val, float extraweight=1.) {
                    h.Fill(region, lep1.id, lep2.id, val, weight*extraweight);
                };
                auto do_fill2D = [region, lep1, lep2, weight](HistCol2D& h, float valx, float valy) {
                    h.Fill(region, lep1.id, lep2.id, valx, valy, weight);
                };
                do_fill(h_met, met);
                do_fill(h_nleps, nleps);
                if (nleps == 2)
                    do_fill(h_ss_type,vdilep_good_ss_sorted.at(0).type);
                else
                    do_fill(h_ml_type,int((sgn(lep1.id)+sgn(lep2.id)+sgn(lep3.id)+3)/2));
                do_fill(h_njets, njets);
                do_fill(h_ptl1, lep1.pt);
                do_fill(h_ptl2, lep2.pt);
                do_fill(h_ptj1, jets_ro.at(0).pt);
                do_fill(h_ptj2, jets_ro.at(1).pt);
                do_fill(h_mjj, mjj_max);
                do_fill(h_mjj_deta, mjj_max_deta);
                do_fill(h_mjj_dphi, mjj_max_dphi);
                do_fill(h_mjj_dr, mjj_max_dr);
                if (nleps == 2)
                    do_fill(h_ss_mll, (lep1.v+lep2.v).mass());
                else {           
                    std::vector<Lepton> good_vls = {lep1,lep2};
                    if (nleps>2) good_vls.push_back(lep3);   
                    if (nleps>3) good_vls.push_back(lep4);   
                    for (unsigned int i = 0; i < good_vls.size(); i++) {
                        Lepton l1 = good_vls.at(i);
                        for (unsigned int j = i+1; j < good_vls.size(); j++) {
                            Lepton l2 = good_vls.at(j);
                            if (l1.id*l2.id < 0 and abs(l1.id) == abs(l2.id)) do_fill(h_ml_mll, (l1.v+l2.v).mass());
                        }
                    }
                }   
            };

            //bool MossfOnZ = fabs(MOSSF(lep)-MZ)<15.;
            if  ( nleps > 1 ){ 
                tree_met = met;
                tree_weight = weight;	
                tree_nleps=nleps;	
                tree_type=type;
                out_tree->Fill();	

                fill_region("br",weight);
                if (fabs(mjj_max_deta) > 2.5) fill_region("br_deta", weight);
                if (nleps == 2) {
                    fill_region("ssbr",weight);
                    if (fabs(mjj_max_deta) > 2.5) fill_region("ssbr_deta",weight);
                    if (mjj_max>500) {
                        fill_region("ssbr_mjj500",weight);
                        if (fabs(mjj_max_deta) > 2.5) fill_region("ssbr_mjj500_deta",weight);
                        if (fabs(mjj_max_deta) > 2.5 and njets>2) fill_region("ssbr_mjj500_deta_njgt2",weight);
                        if (njets > 2) fill_region("ssbr_mjj500_njgt2",weight);
                    }
                } else {
                    fill_region("mlbr",weight);
                    if (mjj_max>500) {
                        fill_region("mlbr_mjj500",weight);
                        if (fabs(mjj_max_deta) > 2.5) fill_region("mlbr_mjj500_deta",weight);
                        if (fabs(mjj_max_deta) > 2.5 and njets>2) fill_region("mlbr_mjj500_deta_njgt2",weight);
                        if (njets > 2) fill_region("mlbr_mjj500_njgt2",weight);
                    }
                }
            }
        }
    }//file loop

    f1->cd();
    for (HistCol* coll : registry) coll->Write();
    for (HistCol2D* coll : registry2D) coll->Write();
    h_cuts->Write();
    out_tree->Write();
    f1->Close();
    delete h_cuts;
    if (dfile.is_open()) dfile.close();
    if (!quiet) cout << "\n Done!" << endl;
    return 0;
}
