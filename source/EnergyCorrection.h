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

    void SetGeneratorType(const std::string &type) { m_generatortype = type; }

    void SetMinEta(float min) { mineta = min; }
    void SetMaxEta(float max) { maxeta = max; }

    void SetReweightHeavierBaryons(bool reweight) { reweightheavierbaryons = reweight; }

private:
    std::string m_HitNodeName {"G4HIT_CEMC"};
    std::string m_generatortype {"HIJING"};
    
    bool m_upweighttruth = false;

    int m_npart =  -1;

    float mineta = -2.5;
    float maxeta = 2.5;

    bool reweightheavierbaryons = true;

    static const int ncentbins = 5;
    TH1F *h_pimi[5] = {nullptr};
    TH1F *h_pip[5] = {nullptr};
    TH1F *h_p[5] = {nullptr};
    TH1F *h_pbar[5] = {nullptr};
    TH1F *h_kp[5] = {nullptr};
    TH1F *h_kmi[5] = {nullptr};

    TH1F *h_lambda[5] = {nullptr};
    TH1F *h_lambdabar[5] = {nullptr}; 

    
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
