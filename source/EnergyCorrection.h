// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ENERGYCORRECTION_H
#define ENERGYCORRECTION_H

#include <fun4all/SubsysReco.h>
#include <TH1.h>

#include <string>
#include <iostream>

class PHCompositeNode;

class EnergyCorrection : public SubsysReco
{
public:
    EnergyCorrection(const std::string &name = "EnergyCorrection");

    ~EnergyCorrection() override;

    int Init(PHCompositeNode *topNode) override;

    /** Called for first event when run number is known.
        Typically this is where you may want to fetch data from
        database, because you know the run number. A place
        to book histograms which have to know the run number.
     */

    /** Called for each event.
        This is where you do the real work.
     */

    int process_event(PHCompositeNode *topNode) override;

    /// Called at the end of all processing.
    int End(PHCompositeNode *topNode) override;

    void SetHitNodeName(const std::string &name) { m_HitNodeName = name; }
    void SetUpweightTruth(bool upweight) { m_upweighttruth = upweight; }
    void SetRapidityDep(bool rapidity = true) { rapiditydep = rapidity; }

    void SetGeneratorType(const std::string &type) { m_generatortype = type; }

    void SetMinEta(float min) { mineta = min; }
    void SetMaxEta(float max) { maxeta = max; }

    void SetReweightHeavierBaryons(bool reweight) { reweightheavierbaryons = reweight; }

private:
    std::string m_HitNodeName {"G4HIT_CEMC"};
    std::string m_generatortype {"HIJING"};
    
    bool m_upweighttruth = false;
    bool rapiditydep = false;

    int m_npart =  -1;

    float mineta = -2.5;
    float maxeta = 2.5;

    bool reweightheavierbaryons = true;

    static const int ncentbins = 5;
    TH1F *h_pimi[ncentbins] = {nullptr};
    TH1F *h_pip[ncentbins] = {nullptr};
    TH1F *h_p[ncentbins] = {nullptr};
    TH1F *h_pbar[ncentbins] = {nullptr};
    TH1F *h_kp[ncentbins] = {nullptr};
    TH1F *h_kmi[ncentbins] = {nullptr};

    TH1F *h_lambda[ncentbins] = {nullptr};
    TH1F *h_lambdabar[ncentbins] = {nullptr};   

    //pi
    std::vector<std::vector<float>> rapIntervalsPi = {{-0.1,0.},{0.,0.1},{0.4,0.6},{0.6,0.8},{0.8,1.0}, {1.0,1.2},{1.2,1.4},{2.1,2.3},{2.4,2.6},{3.0,3.1},{3.1,3.2},{3.2,3.3}, {3.3,3.4},{3.4,3.66}};
    //K
    std::vector<std::vector<float>> rapIntervalsK = {{-0.1,0.},{0.,0.1},{0.4,0.6},{0.6,0.8},{0.8,1.0}, {1.0,1.2},{2.0,2.2},{2.3,2.5},{2.9,3.0},{3.0,3.1},{3.1,3.2},{3.2,3.4}};
    //p
    std::vector<std::vector<float>> rapIntervalsP = {{-0.1,0.1},{0.75,0.95},{1.7,2.4},{2.7,3.1}};
    //pbar
    std::vector<std::vector<float>> rapIntervalsPbar = {{-0.1,0.1},{0.75,0.95},{1.7,2.4},{2.7,3.1}};

    //rapidity dependent correction
    static const int Pihistosize = 14;
    static const int Khistosize = 12;
    static const int Phistosize = 4;
    static const int Pbarhistosize = 4;

    TH1F *hPiplus[Pihistosize] = {nullptr};
    TH1F *hPiminus[Pihistosize] = {nullptr};
    TH1F *hKplus[Khistosize] = {nullptr};
    TH1F *hKminus[Khistosize] = {nullptr};
    TH1F *hP[Phistosize] = {nullptr};
    TH1F *hPbar[Pbarhistosize] = {nullptr};
    TH1F *hPratio = nullptr;
    TH1F *hPbarratio = nullptr;

    const std::vector<std::vector<float>>& getRapidityIntervals(int pid) {
    switch(pid) {
        case 211: return rapIntervalsPi;
        case -211: return rapIntervalsPi;
        case 321: return rapIntervalsK;
        case -321: return rapIntervalsK;
        case 2212: return rapIntervalsP;
        case -2212: return rapIntervalsPbar;
        default: throw std::invalid_argument("Invalid particle ID");
    }
}   

    TH1F* getHist(int pid, int i) {
    switch(pid) {
        case 211: return hPiplus[i];
        case -211: return hPiminus[i];
        case 321: return hKplus[i];
        case -321: return hKminus[i];
        case 2212: return hP[i];
        case -2212: return hPbar[i];

        default: throw std::invalid_argument("Invalid particle ID");
    }

    }

    float findrapscale(int pid, float pt, float y){
        float scale = 1;
        y = std::abs(y);
        const auto& rapIntervals = getRapidityIntervals(pid);
        int Lowerbin = -1;
        int Upperbin = 1000;

        std::vector<float> meanrapidity;
        for (int i = 0; i < (int)rapIntervals.size(); i++)
        {
            meanrapidity.push_back(abs(rapIntervals[i][0] + rapIntervals[i][1]) / 2);
        }
        // find which two bins the y falls in between, and interpolate
        for (int i = 0; i < (int) meanrapidity.size() - 1; i++)
        {

            if (y >= meanrapidity[i] && y <= meanrapidity[i + 1])
            {
                Lowerbin = i;
                Upperbin = i + 1;
                break;
            }
        }
        if(y < meanrapidity[0]){
            Lowerbin = -1;
            Upperbin = 0;
        }
        if(y > meanrapidity[meanrapidity.size()-1]){
            Lowerbin = meanrapidity.size()-1;
            Upperbin = 1000;
        }
        //interpolation
        if (Lowerbin == -1 && Upperbin == 0)
        {
            TH1F *h = getHist(pid, Upperbin);
            scale = h->Interpolate(pt);
        }
        else if (Upperbin == 1000)
        {
            TH1F *h = getHist(pid, Lowerbin);
            scale = h->Interpolate(pt);
        }
        else
        {
            TH1F *h1 = getHist(pid, Lowerbin);
            TH1F *h2 = getHist(pid, Upperbin);
            scale = (h1->Interpolate(pt) * (meanrapidity[Upperbin] - y) + h2->Interpolate(pt) * (y - meanrapidity[Lowerbin])) / (meanrapidity[Upperbin] - meanrapidity[Lowerbin]);
        }
        return scale;
    }

    float findrapcorrection(int pid, float pt, float y, int npart)
    {
        float scale = 1;
        if (pid == 211 || pid == -211 || pid == 321 || pid == -321 || pid == 2212 || pid == -2212)
        {
            scale = findrapscale(pid, pt, y);
        }
        else if(pid == 111){
            scale = (findrapscale(211, pt, y) + findrapscale(-211, pt, y)) / 2.;
        }
        else if(pid == 130 || pid == 310 || pid == 311){
            scale = (findrapscale(321, pt, y) + findrapscale(-321, pt, y)) / 2.;
        }
        //if proton or neutron
        else if (pid == 2212 || pid == 2112)
        {
            scale = findcorrection(npart,2212, pt) * hPratio->Interpolate(std::abs(y)) / hPratio->Interpolate(0);
        }
        //if antiproton or antineutron
        else if (pid == -2212 || pid == -2112)
        {
            scale = findcorrection(npart,-2212, pt) * hPbarratio->Interpolate(std::abs(y)) / hPratio->Interpolate(0);
        }
        //other baryons just use PHENIX and STAR data
        else if(pid > 2000 && pid < 4000){
            scale = findcorrection(npart,2212, pt);
        }
        else if(pid < -2000 && pid > -4000){
            scale = findcorrection(npart,-2212, pt);
        }
        return scale;
    }





    
    float avgcent[ncentbins] = {325.8, 236.1, 141.5, 61.6, 14.7};
    float avgcentlambda[ncentbins] = {356.192, 238.602, 144.272, 64.9728, 23.566};
    

    float findcorrection(int npart, int pid, float pt)
    {
        float weight[ncentbins] = {0};
        float weightlambda[ncentbins] = {0};
        float scale = 0;
        if (npart > avgcentlambda[0] || npart < avgcentlambda[ncentbins - 1])
        {
            if (npart > avgcentlambda[0])
                weightlambda[0] = 1;
            if (npart < avgcentlambda[ncentbins - 1])
                weightlambda[ncentbins - 1] = 1;           
        }
        else{

             // use interpolation here
            // first find which two bins the npart falls in between

            int lowerBinlambda = -1;
            int upperBinlambda = -1;


            for (int i = 0; i < ncentbins - 1; i++)
            {
                if (npart <= avgcentlambda[i] && npart >= avgcentlambda[i + 1])
                {
                    lowerBinlambda = i;
                    upperBinlambda = i + 1;
                    break;
                }
            }

            weightlambda[upperBinlambda] = (avgcentlambda[lowerBinlambda] - npart) / (avgcentlambda[lowerBinlambda] - avgcentlambda[upperBinlambda]);
            weightlambda[lowerBinlambda] = (npart - avgcentlambda[upperBinlambda]) / (avgcentlambda[lowerBinlambda] - avgcentlambda[upperBinlambda]);
            
        }


        if (npart > avgcent[0] || npart < avgcent[ncentbins - 1])
        {
            if (npart > avgcent[0])
                weight[0] = 1;
            if (npart < avgcent[ncentbins - 1])
                weight[ncentbins - 1] = 1;
        }
        else
        {
            // use interpolation here
            // first find which two bins the npart falls in between
            int lowerBin = -1;
            int upperBin = -1;

        

            for (int i = 0; i < ncentbins - 1; i++)
            {
                if (npart <= avgcent[i] && npart >= avgcent[i + 1])
                {
                    lowerBin = i;
                    upperBin = i + 1;
                    break;
                }
            }

         
            // interpolate
            weight[upperBin] = (avgcent[lowerBin] - npart) / (avgcent[lowerBin] - avgcent[upperBin]);
            weight[lowerBin] = (npart - avgcent[upperBin]) / (avgcent[lowerBin] - avgcent[upperBin]);

           // print all weights and the sum
        }
        // if pi0 then use the average of pi+ and pi-
        if (pid == 111)
        {
            for (int i = 0; i < ncentbins; i++)
            {
                scale += weight[i] * (h_pip[i]->Interpolate(pt) + h_pimi[i]->Interpolate(pt)) / 2.;
            }
        }
        else if (pid == 130 || pid == 310 || pid == 311)
        {
            for (int i = 0; i < ncentbins; i++)
            {
                scale += weight[i] * (h_kmi[i]->Interpolate(pt) + h_kp[i]->Interpolate(pt)) / 2.;
            }
        }
        else if (pid == 211)
        {
            // loop over cent bins
            for (int i = 0; i < ncentbins; i++)
            {

                scale += weight[i] * h_pip[i]->Interpolate(pt);
            }
        }
        else if (pid == -211)
        {
            // loop over cent bins
            for (int i = 0; i < ncentbins; i++)
            {
                scale += weight[i] * h_pimi[i]->Interpolate(pt);
            }
        }
      
        // k+

        else if (pid == 321)
        {
            // loop over cent bins
            for (int i = 0; i < ncentbins; i++)
            {
                scale += weight[i] * h_kp[i]->Interpolate(pt);
            }
        }
        // k-
        else if (pid == -321)
        {
            // loop over cent bins
            for (int i = 0; i < ncentbins; i++)
            {

                scale += weight[i] * h_kmi[i] -> Interpolate(pt);
            }
        }
        //lambda and sigmas
        else if (pid == 3122 || pid == 3222 || pid == 3212 || pid == 3112)
        {
            // loop over cent bins
            for (int i = 0; i < ncentbins; i++)
            {
                scale += weightlambda[i] * h_lambda[i]->Interpolate(pt);
            }
            std::cout<<"pid = "<<pid<<" pt = "<<pt<<" scale = "<<scale<<std::endl;
            
        }
        //antilambda and antisigmas
        else if (pid == -3122 || pid == -3222 || pid == -3212 || pid == -3112)
        {
            // loop over cent bins
            for (int i = 0; i < ncentbins; i++)
            {
                scale += weightlambda[i] * h_lambdabar[i]->Interpolate(pt);
            }
        }

        else if (pid > 2000 && pid < 4000)
        {

            // loop over cent bins
            for (int i = 0; i < ncentbins; i++)
            {
                scale += weight[i] * h_p[i]->Interpolate(pt);
            }
            if (!reweightheavierbaryons)
            {
                scale = 1;
            }
        }
        else if (pid < -2000 && pid > -4000)
        {

            // loop over cent bins
            for (int i = 0; i < ncentbins; i++)
            {
                scale += weight[i] * h_pbar[i]->Interpolate(pt);
            }
            if (!reweightheavierbaryons)
            {
                scale = 1;
            }
        }
        else
        {
            scale = 1;
        }
        return scale;
    }
};

#endif // ENERGYCORRECTION_H
