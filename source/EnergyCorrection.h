// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ENERGYCORRECTION_H
#define ENERGYCORRECTION_H

#include <fun4all/SubsysReco.h>
#include <TH1.h>

#include <string>

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
    void SetUpWeightTruth(bool upweight) { m_upweighttruth = upweight; }

private:
    std::string m_HitNodeName = "G4HIT_CEMC";
    bool m_upweighttruth = false;

    int m_npart =  -1;

    static const int ncentbins = 5;
    TH1F *h_pimi[5] = {nullptr};
    TH1F *h_pip[5] = {nullptr};
    TH1F *h_p[5] = {nullptr};
    TH1F *h_pbar[5] = {nullptr};
    TH1F *h_kp[5] = {nullptr};
    TH1F *h_kmi[5] = {nullptr};

    
    const float avgcent[ncentbins] = {325.8, 236.1, 141.5, 61.6, 14.7};
    

    float findcorrection(int npart, int pid, float pt)
    {
        float weight[ncentbins] = {0};
        float scale = 0;

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
        // if pi0 then 50% chance use the correction for pi+ and 50% chance use the correction for pi-
        if (pid == 111)
        {
            pid = 211;
            if (rand() % 2 == 0)
            {
                pid = -211;
            }
        }
        if (pid == 130 || pid == 310 || pid == 311)
        {
            pid = 211;
            if (rand() % 2 == 0)
            {
                pid = -211;
            }
        }
        if (pid == 211)
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
        else if (pid > 2000 && pid < 4000)
        {

            // loop over cent bins
            for (int i = 0; i < ncentbins; i++)
            {
                scale += weight[i] * h_p[i]->Interpolate(pt);
            }
        }
        else if (pid < -2000 && pid > -4000)
        {

            // loop over cent bins
            for (int i = 0; i < ncentbins; i++)
            {
                scale += weight[i] * h_pbar[i]->Interpolate(pt);
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
        else
        {
            scale = 1;
        }
       
        return scale;
    }
};

#endif // ENERGYCORRECTION_H
