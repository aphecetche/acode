/// macro to get an idea of the number of clusters / cm2 
/// from the cluster maps in the QAresults.root

#include <vector>
#include "TFile.h"
#include "TGrid.h"
#include "TH2.h"
#include "TString.h"
#include "TCanvas.h"
#include "Riostream.h"
#include "TPolyLine.h"

struct IntegrationRange
{
    double xmin; 
    double xmax;
    double ymin;
    double ymax;

    double Surface() { return (xmax-xmin)*(ymax-ymin); }
};

void ClusterDensityFromQA(const char* qaFile="")
{
    if ( TString(qaFile).BeginsWith("alien"))
    {
        TGrid::Connect("alien://");
    }

    TFile* f = TFile::Open(qaFile);


    std::vector<IntegrationRange> rangesLow;
    std::vector<IntegrationRange> rangesHigh;

    rangesLow.push_back(IntegrationRange{-10,10, 80, 85});
    rangesLow.push_back(IntegrationRange{-10,10, 80, 85});
    rangesLow.push_back(IntegrationRange{-10,10, 95,100});
    rangesLow.push_back(IntegrationRange{-10,10, 95,100});

    for ( int i = 0; i < 4; ++i )
    {
        IntegrationRange ref = rangesLow[i];

        TObjArray* a = static_cast<TObjArray*>(f->Get("MUON_QA/expert"));

        TH2* h = static_cast<TH2*>(a->FindObject(Form("hClusterHitMapInCh%d",i+1)));

        // derive other symetric ranges from that one
        
        std::vector<IntegrationRange> ranges;

        double ysize = ref.ymax - ref.ymin;
        double xsize = (ref.xmax - ref.xmin)/2.0;

        ranges.push_back(ref);
        ranges.push_back(IntegrationRange{ref.xmin,ref.xmax,-ref.ymax,-ref.ymin});
        ranges.push_back(IntegrationRange{ref.ymin+ysize/2.0,ref.ymax,ref.xmin,ref.xmax});

        TCanvas* c = new TCanvas(Form("Chamber%d",i+1),Form("Chamber%d",i+1));

        h->Draw("colz");

        std::cout << "CHAMBER " << i+1 << " LOW = ";

        for ( auto r : ranges )
        {
            double count = h->Integral(
                h->GetXaxis()->FindBin(r.xmin),
                h->GetXaxis()->FindBin(r.xmax),
                h->GetYaxis()->FindBin(r.ymin),
                h->GetYaxis()->FindBin(r.ymax)
                );

            std::cout << " " << count << "(" << r.Surface() << " cm^2)";
            std::vector<double> x = { r.xmin,r.xmax,r.xmax,r.xmin,r.xmin };
            std::vector<double> y = { r.ymax,r.ymax,r.ymin,r.ymin,r.ymax };

            TPolyLine* l = new TPolyLine(5,&x[0],&y[0]);
            l->SetLineColor(1);
            l->SetLineStyle(9);
            l->Draw();
        }

        std::cout << std::endl;

    }
}

