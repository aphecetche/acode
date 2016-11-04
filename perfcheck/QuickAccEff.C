#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBManager.h"
#include "AliESDEvent.h"
#include "AliESDMuonCluster.h"
#include "AliESDMuonTrack.h"
#include "AliGeomManager.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONCDB.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONPadStatusMaker.h"
#include "AliMUONRecoParam.h"
#include "AliMUONTrackerData.h"
#include "AliMUONVCluster.h"
#include "AliMUONVStore.h"
#include "AliMpBusPatch.h"
#include "AliMpCDB.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "AliMpDEManager.h"
#include "AliMpDetElement.h"
#include "AliMpManuIterator.h"
#include "AliMpPad.h"
#include "AliMpSegmentation.h"
#include "AliMpVSegmentation.h"
#include "Riostream.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGraph.h"
#include "TH1.h"
#include "TTree.h"
#include <cassert>
#include <map>
#include <set>
#include <vector>

#include "TObjectTable.h"


void ReadManuStatus(const char* inputfile,
    std::map<int,std::vector<UInt_t> >& manuStatusForRuns);

typedef std::pair<int,int> ManuPair;

UInt_t ENCODE(Int_t a16, Int_t b16)
{
    return ( ( a16 & 0x0000FFFF ) << 16 ) | ( b16 & 0x0000FFFF );
}

Int_t DECODELOW(UInt_t e)
{
    return e & 0x0000FFFF;
}

Int_t DECODEHIGH(UInt_t e)
{
    return (e & 0xFFFF0000) >> 16;
}

void DECODE(UInt_t e, Int_t& a, Int_t& b)
{
    a = DECODEHIGH(e);
    b = DECODELOW(e);
}

void ReadIntegers(const char* filename,
        std::vector<int>& integers,
        Bool_t resetVector=kTRUE)
{
    /// Read integers from filename, where integers are either
    /// separated by "," or by return carriage

    if ( gSystem->AccessPathName(gSystem->ExpandPathName(filename))==kTRUE )
    {
        return;
    }
    std::ifstream in(gSystem->ExpandPathName(filename));
    int i;

    std::set<int> runset;

    if (!resetVector)
    {
        for ( std::vector<int>::size_type j = 0; j < integers.size(); ++ j )
        {
            runset.insert(integers[j]);
        }
    }

    char line[10000];

    in.getline(line,10000,'\n');

    TString sline(line);

    if (sline.Contains(","))
    {
        TObjArray* a = sline.Tokenize(",");
        TIter next(a);
        TObjString* s;
        while ( ( s = static_cast<TObjString*>(next()) ) )
        {
            runset.insert(s->String().Atoi());
        }
        delete a;
    }
    else
    {
        runset.insert(sline.Atoi());

        while ( in >> i )
        {
            runset.insert(i);
        }
    }

    integers.clear();

    for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it )
    {
        integers.push_back((*it));
    }

    std::sort(integers.begin(),integers.end());
}

void GetRunList(const char* runlist, std::vector<int>& vrunlist)
{
    // Read the runlist from an ASCII file or a comma separated list
    // or a space separated list

    vrunlist.clear();

    if ( TString(runlist).Contains(",") || TString(runlist).Contains(" ") )
    {
        TObjArray* runs = 0x0;
        if ( TString(runlist).Contains(",") )
        {
            runs = TString(runlist).Tokenize(",");
        }
        else
        {
            runs = TString(runlist).Tokenize(" ");
        }
        TIter next(runs);
        TObjString* s;
        std::set<int> runset;

        while ( ( s = static_cast<TObjString*>(next()) ) )
        {
            runset.insert(s->String().Atoi());
        }

        for ( std::set<int>::const_iterator it = runset.begin(); it != runset.end(); ++it )
        {
            vrunlist.push_back((*it));
        }

        std::sort(vrunlist.begin(),vrunlist.end());

        delete runs;
    }
    else
    {
        ReadIntegers(runlist,vrunlist);
    }
}

struct ClusterLocation
{
    ClusterLocation(int b=0, int nb=0)
        : mBendingManuIx(b), mNonBendingManuIx(nb) 
    {}

    int DetElemId() const;

    int BendingManuIndex() const { return mBendingManuIx; }
    int NonBendingManuIndex() const { return mNonBendingManuIx; }

    friend std::ostream& operator<<(std::ostream& out,const ClusterLocation& cl);

private:

    int mBendingManuIx;
    int mNonBendingManuIx;
};

struct CompactTrack
{
    CompactTrack(Double_t x=0.0, Double_t y=0.0, Double_t z=0.0) 
        : mPx(x), mPy(y), mPz(z), mClusters() {}
    Double_t Px() const { return mPx; }
    Double_t Py() const { return mPy; }
    Double_t Pz() const { return mPz; }

    Double_t mPx,mPy,mPz;
    std::vector<ClusterLocation> mClusters;

    friend std::ostream& operator<<(std::ostream& out,
            const CompactTrack& ct);

};

struct CompactEvent
{
    CompactEvent() : mTracks() {}
    std::vector<CompactTrack> mTracks;
};

struct CompactMapping
{
    CompactMapping() : mManuIds(), mManuMap() {}

    // array containing the 32bits encoded
    // pair (detElemId,manuId) for each
    // absolute manu index
    // mManuIds[abs]=(de << 16) & local
    std::vector<UInt_t> mManuIds;

    // associate each encoded pair (detElemId,manuId)
    // to an index in mManuIds
    // (i.e. reverse structure of mManuIds
    std::map<int,int> mManuMap;

    friend std::ostream& operator<<(std::ostream& os, const CompactMapping& cm);

    Int_t GetDetElemIdFromAbsManuIndex(Int_t index) const
    {
        return GetDetElemIdFromAbsManuId(AbsManuId(index));
    }
        
    Int_t AbsManuId(UInt_t index) const 
    {
        if ( index>mManuIds.size())
        {
            std::cout << index << " > " << mManuIds.size()
                << std::endl;
        }
        return mManuIds[index];
    }
    
    Int_t GetManuIdFromAbsManuId(UInt_t absManuId) const
    {
        return DECODELOW(absManuId);
    }

    Int_t GetDetElemIdFromAbsManuId(UInt_t absManuId) const
    {
        return DECODEHIGH(absManuId);
    }
};

UInt_t MANUBADPEDMASK = ( 1 << 0 );
UInt_t MANUBADHVMASK = ( 1 << 1 );
UInt_t MANUBADLVMASK =  (1 << 2 );
UInt_t MANUBADOCCMASK = ( 1 << 3 );
UInt_t MANUOUTOFCONFIGMASK = ( 1 << 4 );

std::string CauseAsString(UInt_t cause)
{
    std::string rv = "";

    if ( cause & MANUBADPEDMASK ) rv += "_ped";
    if ( cause & MANUBADHVMASK ) rv += "_hv";
    if ( cause & MANUBADLVMASK ) rv += "_lv";
    if ( cause & MANUBADOCCMASK ) rv += "_occ";
    if ( cause & MANUOUTOFCONFIGMASK ) rv += "_config";

    return rv;
}

CompactMapping* GetCompactMapping()
{
    static CompactMapping* cm(0x0);

    if (!cm)
    {
        cm = new CompactMapping;

        AliCDBManager* man = AliCDBManager::Instance();
        if (!man->IsDefaultStorageSet())
        {
            man->SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2016/OCDB");
            man->SetRun(264000);
        }
        AliMpCDB::LoadAll();

        // first get a sorted list of the detElemId

        AliMpDEIterator deit;
        deit.First();
        std::vector<int> deids;
        while (!deit.IsDone())
        {
            Int_t detElemId = deit.CurrentDEId();

            if ( AliMpDEManager::GetStationType(detElemId) != AliMp::kStationTrigger )
            {
                deids.push_back(detElemId);
            }
            deit.Next();
        }

        std::sort(deids.begin(),deids.end());

        // then for each DE get one sorted list of manuIds
        int totalNofManus(0);

        for ( std::vector<int>::size_type i = 0; i < deids.size(); ++i )
        {
            AliMpManuIterator it;
            Int_t detElemId, manuId;
            std::vector<int> bendingManuids;
            std::vector<int> nonBendingManuids;

            // get the list of manu ids on both planes
            // for this detection element
            while ( it.Next(detElemId,manuId) )
            {
                if ( detElemId == deids[i] )
                {
                    if ( manuId >= 1024 )
                    {
                        nonBendingManuids.push_back(manuId);
                    }
                    if ( manuId < 1024 )
                    {
                        bendingManuids.push_back(manuId);
                    }
                }
            }

            detElemId = deids[i];

            // insure manuids are sorted (should be the case
            // already, though)
            std::sort(bendingManuids.begin(),bendingManuids.end());
            std::sort(nonBendingManuids.begin(),nonBendingManuids.end());

            // add those manus to the global array of AbsManuIds
            // and update the relevant "links"
            

            std::vector<int>& allManuOfThisDE = bendingManuids;
            allManuOfThisDE.insert(allManuOfThisDE.end(),
                    nonBendingManuids.begin(),
                    nonBendingManuids.end());
          
            UInt_t ix = cm->mManuIds.size();

            for ( std::vector<int>::size_type i = 0; i < allManuOfThisDE.size(); ++i )
            {
                Int_t manuId = allManuOfThisDE[i];
                UInt_t encodedManu = ENCODE(detElemId,manuId);
                cm->mManuMap[encodedManu]=ix;
                cm->mManuIds.push_back(encodedManu);
                ++ix;
            }

            totalNofManus += allManuOfThisDE.size();

            // std::cout << Form("DE %04d (%3d) ",deids[i],
            //         allManuOfThisDE.size());
            //
            // for ( std::vector<int>::size_type j = 0; 
            //         j < allManuOfThisDE.size(); ++j )
            // {
            //     std::cout << Form("%04d ",allManuOfThisDE[j]);
            // }
            // std::cout << std::endl;
        }

        // std::cout << "Total number of manus : " << totalNofManus << std::endl;
        assert(totalNofManus==16828);
    }
    return cm;
}

std::ostream& operator<<(std::ostream& os,
        const CompactMapping& cm)
{
    // consistency check
    for ( unsigned int i = 0; i < cm.mManuIds.size(); ++i )
    {
        Int_t absManuId = cm.AbsManuId(i);

        std::cout << Form("ENCODEDMANUID[%6d]=%04d : DE %04d MANU %04d",
                i,absManuId,
                cm.GetDetElemIdFromAbsManuId(absManuId),
                cm.GetManuIdFromAbsManuId(absManuId))
            << std::endl;
    }

    return os;
}

std::ostream& operator<<(std::ostream& out, const ClusterLocation& cl)
{
    out << Form(" DE %04d B %6d NB %6d ",
            cl.DetElemId(),cl.BendingManuIndex(),cl.NonBendingManuIndex());
    return out;
}

Int_t ClusterLocation::DetElemId() const
{
    // FIXME: not ideal to have this method 
    // depends on the compact mapping object.
    // but still usefull for debug
    // or is it just fine ?
    
    CompactMapping* cm = GetCompactMapping();
    Int_t absManuIndex = BendingManuIndex();
    if ( absManuIndex < 0 ) 
    {
        absManuIndex = NonBendingManuIndex();
    }
    
    assert(absManuIndex>=0);
    assert(absManuIndex<=static_cast<int>(cm->mManuIds.size()));

    return cm->GetDetElemIdFromAbsManuIndex(absManuIndex);
}

std::ostream& operator<<(std::ostream& out, const CompactTrack& track)
{
    out << "CompactTrack ";
    for ( std::vector<ClusterLocation>::size_type i = 0;
            i < track.mClusters.size(); ++i )
    {
        const ClusterLocation& cl = track.mClusters[i];
        out << cl << " ";
    }
    return out;
}

AliMUONGeometryTransformer* Transformer()
{
    static AliMUONGeometryTransformer* t = 0x0;
    if (!t)
    {
        t = new AliMUONGeometryTransformer();
        t->LoadGeometryData();
    }
    return t;
}

Int_t FindManuAbsIndex(Int_t detElemId, Int_t manuId)
{
    UInt_t encoded = ENCODE(detElemId,manuId);
    CompactMapping* cm = GetCompactMapping();

    return cm->mManuMap[encoded];
}

void GetClusterLocation(Int_t detElemId,
        Double_t xg, 
        Double_t yg, 
        Double_t zg,
        Int_t& bendingManuAbsIndex,
        Int_t& nonBendingManuAbsIndex)
{
    /// Get the pad under the center of the cluster
    Double_t x,y,z;

    Transformer()->Global2Local(detElemId,
            xg,yg,zg,x,y,z);

    AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

    const AliMpVSegmentation* segB = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,de->GetCathodType(AliMp::kBendingPlane));
    const AliMpVSegmentation* segNB = AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,de->GetCathodType(AliMp::kNonBendingPlane));

    AliMpPad padB = segB->PadByPosition(x,y);
    AliMpPad padNB = segNB->PadByPosition(x,y);

    Int_t manuBending = padB.GetManuId();
    Int_t manuNonBending = padNB.GetManuId();

    bendingManuAbsIndex=-1;
    nonBendingManuAbsIndex=-1;

    if ( padB.IsValid() )
    { 
        bendingManuAbsIndex = FindManuAbsIndex(de->GetId(),manuBending);
    }

    if ( padNB.IsValid() )
    {
        nonBendingManuAbsIndex = FindManuAbsIndex(de->GetId(),manuNonBending);
    }
}

void ConvertEvent(AliESDEvent& esd, CompactEvent& compactEvent)
{
    compactEvent.mTracks.clear();

    for ( Int_t i = 0; i < esd.GetNumberOfMuonTracks(); ++i )
    {
        AliESDMuonTrack* track = esd.GetMuonTrack(i);

        if  (!track->ContainTrackerData()) continue;

        CompactTrack compactTrack(track->Px(),track->Py(),track->Pz());

        for ( Int_t j = 0; j < track->GetNClusters(); ++j )
        {
            UInt_t id = track->GetClusterId(j);
            AliESDMuonCluster* cluster = esd.FindMuonCluster(id);
            Int_t b,nb;
            GetClusterLocation(cluster->GetDetElemId(),
                    cluster->GetX(),
                    cluster->GetY(),
                    cluster->GetZ(),
                    b,
                    nb);
            ClusterLocation cl(b,nb);
            compactTrack.mClusters.push_back(cl);
        }

        compactEvent.mTracks.push_back(compactTrack);
    }
}

Bool_t ValidateCluster(const ClusterLocation& cl,
        const std::vector<UInt_t>& manuStatus,
        UInt_t causeMask)
{
    UInt_t bendingMask = 0;
    UInt_t nonBendingMask = 0;

    if ( cl.BendingManuIndex() >= 0 ) 
    {
        assert(cl.BendingManuIndex()<=(int)manuStatus.size());
        bendingMask = manuStatus[cl.BendingManuIndex()];
    }
    if ( cl.NonBendingManuIndex() >= 0 )
    {
        assert(cl.NonBendingManuIndex()<=(int)manuStatus.size());
        nonBendingMask=manuStatus[cl.NonBendingManuIndex()];
    }

    if ( ( bendingMask & causeMask ) ||
            ( nonBendingMask & causeMask ) )
    {
        return kFALSE;
    }
    return kTRUE;
}

Bool_t ValidateTrack(const CompactTrack& track,
        const std::vector<UInt_t>& manuStatus,
        UInt_t causeMask)
{
    /// We remove from the track all the clusters 
    /// located on a bad manu.
    /// Then we consider the track survived if we get
    /// at least one cluster per station

    if ( manuStatus.empty() || causeMask == 0 ) return kTRUE;

    Int_t currentCh;
    Int_t currentSt;
    Int_t previousCh = -1;
    Int_t nChHitInSt4 = 0;
    Int_t nChHitInSt5 = 0;
    UInt_t presentStationMask = 0;
    const UInt_t requestedStationMask = 0x1F;
    const Bool_t request2ChInSameSt45 = kTRUE;

    CompactMapping* cm = GetCompactMapping();

    for ( std::vector<ClusterLocation>::size_type i = 0;
            i < track.mClusters.size(); ++i )
    {
        const ClusterLocation& cl = track.mClusters[i];

        if (!ValidateCluster(cl,manuStatus,causeMask))
        {
            continue;
        }

        Int_t detElemId = cl.DetElemId(); 

        currentCh = AliMpDEManager::GetChamberId(detElemId);
        currentSt = currentCh/2;

        // build present station mask
        presentStationMask |= ( 1 << currentSt );

        // count the number of chambers hit in station 4 that contain cluster(s)
        if (currentSt == 3 && currentCh != previousCh) {
            nChHitInSt4++;
            previousCh = currentCh;
        }

        // count the number of chambers hit in station 5 that contain cluster(s)
        if (currentSt == 4 && currentCh != previousCh) {
            nChHitInSt5++;
            previousCh = currentCh;
        }

    }

    // at least one cluster per requested station
    if ((requestedStationMask & presentStationMask) != requestedStationMask) 
    {
        return kFALSE;
    }

    if (request2ChInSameSt45) 
    {
        // 2 chambers hit in the same station (4 or 5)
        return (nChHitInSt4 == 2 || nChHitInSt5 == 2);
    }
    else 
    {
        // or 2 chambers hit in station 4 & 5 together
        return (nChHitInSt4+nChHitInSt5 >= 2);
    }

    return kTRUE;
}

// void FilterTracks(const std::vector<CompactEvent>& idealEvents,
//         std::vector<CompactEvent>& filteredEvents,
//         const std::map<int,std::set<int> >& badManuList)
// {
//     /// Loop over events and remove the tracks
//     /// that are affected by the holes in the badManuList
//
//     filteredEvents.clear();
//
//     for ( std::vector<CompactEvent>::size_type i = 0;
//             i < idealEvents.size(); ++i )
//     {
//         const CompactEvent& e = idealEvents[i];
//         CompactEvent fe;
//
//         for ( std::vector<CompactTrack>::size_type j = 0;
//                 j < e.mTracks.size(); ++j ) 
//         {
//             const CompactTrack& track = e.mTracks[j];
//             if (ValidateTrack(track,badManuList))
//             {
//                 fe.mTracks.push_back(track);
//             }
//         }
//         filteredEvents.push_back(fe);
//     }
// }

TH1* ComputeMinv(const std::vector<CompactEvent>& events,
        std::vector<UInt_t> manustatus,
        UInt_t causeMask)
{
    TH1* h = new TH1F("hminv","hminv",300,0.0,15.0);

    const double m2 = 0.1056584*0.1056584;

    Int_t nTracks=0;
    Int_t nValidatedTracks = 0;

    for ( std::vector<CompactEvent>::size_type i = 0;
            i < events.size(); ++i )
    {
        const CompactEvent& e = events[i];

        for ( std::vector<CompactTrack>::size_type j = 0;
                j < e.mTracks.size(); ++j ) 
        {
            const CompactTrack& t1 = e.mTracks[j];

            ++nTracks;
            if (!ValidateTrack(t1,manustatus,causeMask)) continue;
            ++nValidatedTracks;

            for ( std::vector<CompactTrack>::size_type k = j+1;
                    k < e.mTracks.size(); ++k )
            {
                const CompactTrack& t2 = e.mTracks[k];

                if (!ValidateTrack(t2,manustatus,causeMask)) continue;

                double p1square = t1.mPx*t1.mPx +
                    t1.mPy*t1.mPy +
                    t1.mPz*t1.mPz;

                double p2square = t2.mPx*t2.mPx +
                    t2.mPy*t2.mPy +
                    t2.mPz*t2.mPz;

                double minv = TMath::Sqrt(2.0*( 
                            m2 
                            + TMath::Sqrt(m2+p1square)*
                            TMath::Sqrt(m2+p2square)
                            - (t1.mPx*t2.mPx+t1.mPy*t2.mPy+
                                t1.mPz*t2.mPz)));
                h->Fill(minv);

            }
        }
    }

    std::cout << Form("nTracks %d nValidated %d",nTracks,
            nValidatedTracks) << std::endl;

    return h;
}

Int_t ConvertESD(const char* inputfile,
        const char* outputfile)
{
    AliCDBManager::Instance()->SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2016/OCDB");

    AliCDBManager::Instance()->SetRun(0);
    AliMpCDB::LoadAll();

    AliGeomManager::LoadGeometry(Form("%s/geometry.root",
                gSystem->DirName(inputfile)));
//    if (!AliGeomManager::ApplyAlignObjsFromCDB("MUON")) return -1;

    TFile* f = TFile::Open(inputfile);
    if (!f->IsOpen()) return -1;

    TTree* tree = static_cast<TTree*>(f->Get("esdTree"));

    AliESDEvent esd;

    esd.ReadFromTree(tree);

    TFile* fout = TFile::Open(outputfile,"recreate");
    TTree* out = new TTree("compactevents","a tree with compacted tracks");
    CompactEvent compactEvent;
    out->Branch("event",&compactEvent);

    for ( Long64_t i = 0; i < tree->GetEntries(); ++i )
    {
        tree->GetEntry(i);

        if ( esd.GetNumberOfMuonTracks() >= 2 )
        {
            ConvertEvent(esd,compactEvent);
            out->Fill();
        }
    }

    out->Write();
    delete fout;

    return 0;
}

UInt_t GetEvents(TTree* tree,std::vector<CompactEvent>& events)
{
    /// Read events from the tree
    events.clear();
    CompactEvent* compactEvent=0x0;
    tree->SetBranchAddress("event",&compactEvent);

    for ( Long64_t i = 0; i < tree->GetEntries(); ++i )
    {
        tree->GetEntry(i);
        events.push_back(*compactEvent);
    }
    return events.size();
}

UInt_t GetEvents(const char* treeFile, std::vector<CompactEvent>& events)
{
    TFile* f = TFile::Open(treeFile);
    if (!f->IsOpen()) return 0;

    TTree* tree = static_cast<TTree*>(f->Get("compactevents"));
    if (!tree) return 0;

    UInt_t rv = GetEvents(tree,events);

    delete f;

    return rv;
}

void GetManuStatus(Int_t runNumber, std::vector<UInt_t>& manustatus, Bool_t print=kFALSE)
{
    AliCDBManager* man = AliCDBManager::Instance();
    if (!man->IsDefaultStorageSet())
    {
        man->SetDefaultStorage("raw://");
    }

    man->SetRun(runNumber);

    if (!AliMpDDLStore::Instance())
    {
        AliMpCDB::LoadAll();
    }
    
    AliMUONCalibrationData cd(runNumber,true);

    AliMUONPadStatusMaker statusMaker(cd);

    AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();

    statusMaker.SetLimits(*recoParam);

    manustatus.resize(16828,0);

    AliMpManuIterator it;
    Int_t detElemId, manuId;

    Int_t pedCheck = (
            AliMUONPadStatusMaker::kPedMeanZero |
            AliMUONPadStatusMaker::kPedMeanTooLow |
            AliMUONPadStatusMaker::kPedMeanTooHigh |
            AliMUONPadStatusMaker::kPedSigmaTooLow |
            AliMUONPadStatusMaker::kPedSigmaTooHigh );

    Int_t hvCheck = (
            AliMUONPadStatusMaker::kHVError |
            AliMUONPadStatusMaker::kHVTooLow |
            AliMUONPadStatusMaker::kHVTooHigh |
            AliMUONPadStatusMaker::kHVChannelOFF |
            AliMUONPadStatusMaker::kHVSwitchOFF );

    Int_t occCheck = (
            AliMUONPadStatusMaker::kManuOccupancyTooHigh
            );

    Int_t lvCheck = ( AliMUONPadStatusMaker::kLVTooLow );

    Int_t ntotal(0);
    Int_t nbad(0);
    Int_t nbadped=0;
    Int_t nbadocc=0;
    Int_t nbadhv=0;
    Int_t nbadlv=0;
    Int_t nmissing=0;
    Int_t nreco=0;

    while ( it.Next(detElemId,manuId) )
    {
        AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);

        UInt_t manuStatus = 0;

        Int_t manubadped=0;
        Int_t manubadocc=0;
        Int_t manubadhv=0;
        Int_t manubadlv=0;
        Int_t manumissing=0;

        for ( Int_t manuChannel = 0; manuChannel < AliMpConstants::ManuNofChannels(); ++manuChannel )
        {
            if ( de->IsConnectedChannel(manuId,manuChannel) )
            {
                ++ntotal;

                UInt_t status = statusMaker.PadStatus(detElemId, manuId, manuChannel);

                // if ( runNumber == 263529 )
                // {
                //     // FIXME: manual hack for testing
                //     if (detElemId >= 1007 && detElemId<=1013)
                //     {
                //         manuStatus |= MANUBADLVMASK;
                //         status = AliMUONPadStatusMaker::BuildStatus(0,0,lvCheck,0);
                //     }
                //     if (detElemId >= 907 && detElemId<=913)
                //     {
                //         manuStatus |= MANUBADLVMASK;
                //         status = AliMUONPadStatusMaker::BuildStatus(0,0,lvCheck,0);
                //     }
                // }
                
                if (!status) continue;

                if ( status & AliMUONPadStatusMaker::BuildStatus(pedCheck,0,0,0) )
                {
                    ++manubadped;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(0,hvCheck,0,0) )
                {
                    ++manubadhv;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(0,0,lvCheck,0) )
                {
                    ++manubadlv;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(0,0,0,occCheck) )
                {
                    ++manubadocc;
                }

                if ( status & AliMUONPadStatusMaker::BuildStatus(AliMUONPadStatusMaker::kMissing,0,0,AliMUONPadStatusMaker::kMissing) )
                {
                    ++manumissing;
                }

            }
            if (manubadped>=0.7*de->NofChannels())
            {
                manuStatus |= MANUBADPEDMASK;
            }

            if ( manubadhv )
            {
                manuStatus |= MANUBADHVMASK;
            }

            if ( manubadlv )
            {
                manuStatus |= MANUBADLVMASK;
            }
            if ( manubadocc ) 
            {
                manuStatus |= MANUBADOCCMASK;
            }
            if ( manumissing) 
            {
                manuStatus |= MANUOUTOFCONFIGMASK;
            }

            Int_t manuAbsIndex = FindManuAbsIndex(detElemId,manuId);
            manustatus[manuAbsIndex] = manuStatus;
        }
    }

    assert(ntotal==1064008);

    if ( print )
    {
        CompactMapping* cm = GetCompactMapping();

        for ( std::vector<UInt_t>::size_type i = 0; i < manustatus.size(); ++i )
        {
            Int_t absManuId = cm->AbsManuId(i);
            Int_t detElemId = cm->GetDetElemIdFromAbsManuId(absManuId);
            Int_t manuId = cm->GetManuIdFromAbsManuId(absManuId);
            if ( manustatus[i])
            {
                std::cout << Form("status[%04d]=%6x (DE %04d MANU %04d) %s",i,manustatus[i],detElemId,manuId,CauseAsString(manustatus[i]).c_str()) << std::endl;
            }
        }
    }

    man->ClearCache();
}

void GetBadManuListFromBPOccupancy(const char* ocdbpath,
        Int_t runNumber,
        std::map<int, std::set<int> >& badManuList)
{
    /// We do rough selection here. 
    /// Assume that all the manus belonging to a buspatch
    /// that has an occupancy above 30% are bad.

    badManuList.clear();

    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage(ocdbpath);
    man->SetRun(runNumber);

    AliMpCDB::LoadAll();

    AliCDBEntry* e = static_cast<AliCDBEntry*>(man->Get("MUON/Calib/OccupancyMap",runNumber));
    AliMUONVStore* occupancyStore = static_cast<AliMUONVStore*>(e->GetObject());

    AliMUONTrackerData td("occ","occ",*occupancyStore);

    TIter nextBP(AliMpDDLStore::Instance()->CreateBusPatchIterator());
    AliMpBusPatch* bp;

    while ( ( bp = static_cast<AliMpBusPatch*>(nextBP())))
    {
        double occ = td.BusPatch(bp->GetId(),2);
        if ( occ>0.3 )
            // || 
            //     ( bp->GetId() >= 1601 &&
            //       bp->GetId() <= 1614 ))
        {
            std::cout << Form("BP %04d OCC %7.2f %%",
                    bp->GetId(),occ*100.0) << std::endl;
            int de = AliMpDDLStore::Instance()->GetDEfromBus(bp->GetId());
            std::set<int> manus;
            for ( int i = 0; i < bp->GetNofManus(); ++i )
            {
                manus.insert(bp->GetManuId(i));
            }
            badManuList[de].insert(manus.begin(),manus.end());
        }
    }
}

void ComputeEvolution(const std::vector<CompactEvent>& events, 
        std::vector<int>& vrunlist,
        const std::map<int,std::vector<UInt_t> >& manuStatusForRuns,
        const char* outputfile)
{
    TH1* h = ComputeMinv(events,std::vector<UInt_t>(),0);
    Int_t b1 = h->GetXaxis()->FindBin(2.8);
    Int_t b2 = h->GetXaxis()->FindBin(3.4);
    Double_t referenceNofJpsi = h->Integral(b1,b2);

    std::vector<TGraph*> gdrop;
    // one graph for each "bad" cause (but on 
    // manu level only)
    // - ped is a bit ill-defined for manu, we consider a manu
    // "bad for ped" if 70% of its channels are bad
    // - manu occupancy
    // - hv
    // - lv
    // - missing (i.e. buspatch removed from configuration)

    std::vector<UInt_t> causes;

    causes.push_back(MANUOUTOFCONFIGMASK);
    causes.push_back(MANUOUTOFCONFIGMASK |
                     MANUBADPEDMASK);
    causes.push_back(MANUOUTOFCONFIGMASK | 
                     MANUBADPEDMASK  | 
                     MANUBADOCCMASK);
    causes.push_back(MANUOUTOFCONFIGMASK | 
                     MANUBADPEDMASK  | 
                     MANUBADOCCMASK  |
                     MANUBADHVMASK);
    causes.push_back(MANUOUTOFCONFIGMASK | 
                     MANUBADPEDMASK  | 
                     MANUBADOCCMASK  |
                     MANUBADHVMASK   |
                     MANUBADLVMASK);
    for ( std::vector<UInt_t>::size_type i = 0; i < causes.size(); ++i )
    {
        TGraph* g = new TGraph(vrunlist.size());
        gdrop.push_back(g);
        g->SetName(Form("acceffdrop%s",CauseAsString(causes[i]).c_str()));
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.5);
    }

    for ( std::vector<int>::size_type i = 0; i < vrunlist.size(); ++i )
    {
        Int_t runNumber = vrunlist[i];

        std::map<int, std::vector<UInt_t> >::const_iterator it = manuStatusForRuns.find(runNumber);

        const std::vector<UInt_t>& manustatus = it->second; 

        for ( std::vector<UInt_t>::size_type icause = 0; icause < causes.size(); ++icause )
        {
            auto nbad =std::count_if(manustatus.begin(),
                    manustatus.end(),
                    [&](int n) { return (n & causes[icause]); }); 
            std::cout << Form("RUN %6d %30s rejected manus = %6d",
                runNumber,
                CauseAsString(causes[icause]).c_str(),
                nbad
                ) << std::endl;
            TH1* h = ComputeMinv(events,manustatus,causes[icause]);
            h->SetName(Form("hminv%6d",runNumber));
            Double_t drop = 100.0*(1.0-h->Integral(b1,b2)/referenceNofJpsi);
            std::cout << Form("RUN %6d %30s AccxEff drop %7.2f %%",
                    runNumber," ",drop) << std::endl;
            gdrop[icause]->SetPoint(i,runNumber,drop);
            delete h;
        }
    }


    TFile* fout = TFile::Open(outputfile,"recreate");
    for ( std::vector<UInt_t>::size_type icause = 0; icause < causes.size(); ++icause )
    {
        gdrop[icause]->Write();
    }
    delete fout;
}

void ComputeEvolutionFromManuStatus(const char* treeFile,
        const char* runList,
        const char* outputfile,
        const char* manustatusfile) 
{
    std::vector<CompactEvent> events;

    if (!GetEvents(treeFile,events))
    {
        return ;
    }

    std::vector<int> vrunlist;
    GetRunList(runList,vrunlist);

    std::map<int,std::vector<UInt_t> > manuStatusForRuns;
   
    ReadManuStatus(manustatusfile,manuStatusForRuns);

    ComputeEvolution(events,vrunlist,manuStatusForRuns,outputfile);
}

void ComputeEvolution(const char* treeFile,
        const char* runList,
        const char* outputfile,
        const char* ocdbpath="local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2016/OCDB")
{
    std::vector<CompactEvent> events;

    if (!GetEvents(treeFile,events))
    {
        return ;
    }

    std::vector<int> vrunlist;
    GetRunList(runList,vrunlist);

    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage(ocdbpath);

    man->SetRun(vrunlist[0]);
    AliMpCDB::LoadAll();

    std::map<int,std::vector<UInt_t> > manuStatusForRuns;

    for ( std::vector<int>::size_type i = 0; i < vrunlist.size(); ++i )
    {
        Int_t runNumber = vrunlist[i];

        std::vector<UInt_t> manustatus;

        GetManuStatus(runNumber,manustatus);

        manuStatusForRuns[runNumber]=manustatus;

        AliCDBManager::Instance()->ClearCache();
    }

    ComputeEvolution(events,vrunlist,manuStatusForRuns,outputfile);
}


void WriteCompactMappingForO2(const char* outputfile)
{
    CompactMapping* cm = GetCompactMapping();

    if (!cm) return;

    std::ofstream out(outputfile,std::ios::binary);

    std::vector<UInt_t> manus;

    assert(cm->mManuIds.size()==16828);

    manus.resize(cm->mManuIds.size());

    // some devices might need the number of channels
    // for each manu (e.g. the occupancy device)
    std::vector<uint8_t> npad;
    npad.resize(cm->mManuIds.size());

    for ( std::vector<UInt_t>::size_type i = 0; i < 
            cm->mManuIds.size(); ++i )
    {
        // "remap" to the VDigit UniqueID convention
        // == detElemId | manuID << 12 | manuChannel << 24 |
        // cathode << 30
        UInt_t id = cm->mManuIds[i];
        Int_t detElemId;
        Int_t manuId;
        DECODE(id,detElemId,manuId);
        manus[i] = detElemId | (manuId << 12);
        AliMpDetElement* de = AliMpDDLStore::Instance()->GetDetElement(detElemId);
        npad[i] = de->NofChannelsInManu(manuId);
    }

    std::vector<UInt_t>::size_type nmanus = manus.size();

    out.write((char*)&nmanus,sizeof(int));
    out.write((char*)&manus[0],nmanus*sizeof(UInt_t));
    out.write((char*)&npad[0],nmanus*sizeof(uint8_t));

    out.close();
}

void WriteManuStatus(const char* runlist, const char* outputfile, Bool_t print=kFALSE)
{
    std::vector<int> vrunlist;
    GetRunList(runlist,vrunlist);

    std::ofstream out(outputfile,std::ios::binary);

    std::vector<int>::size_type nruns = vrunlist.size();

    out.write((char*)&nruns,sizeof(int));
    out.write((char*)&vrunlist[0],nruns*sizeof(int));
    
    for ( std::vector<int>::size_type i = 0; i < vrunlist.size(); ++i )
    {
        Int_t runNumber = vrunlist[i];
        std::vector<UInt_t> manuStatus;
        GetManuStatus(runNumber,manuStatus,print);
        out.write((char*)&manuStatus[0],
                manuStatus.size()*sizeof(int));
        assert(manuStatus.size()==16828);

    std::cout << Form("RUN %6d",runNumber) << std::endl;
    gObjectTable->Print();

    }
    out.close();
}

void ReadManuStatus(const char* inputfile,
    std::map<int,std::vector<UInt_t> >& manuStatusForRuns)
{
    std::ifstream in(inputfile,std::ios::binary);

    int nruns;

    in.read((char*)&nruns,sizeof(int));

    std::cout << "nruns=" << nruns << std::endl;

    std::vector<int> vrunlist;

    vrunlist.resize(nruns,0);

    std::vector<int> manuStatus;

    in.read((char*)&vrunlist[0],sizeof(int)*nruns);

    for ( std::vector<int>::size_type i = 0; i < vrunlist.size(); ++i ) 
    {
        Int_t runNumber = vrunlist[i];
        std::cout << runNumber << " ";
        manuStatus.resize(16828,0);
        in.read((char*)&manuStatus[0],sizeof(int)*manuStatus.size());
        for ( std::vector<int>::size_type j = 0; j < manuStatus.size(); ++j )
        {
            manuStatusForRuns[runNumber].push_back(manuStatus[j]);
        }
    }
    std::cout << std::endl;
}

void CompareWithFullAccEff(const char* fullacceff="/data/EfficiencyJPsiRun_from_astrid.root", const char* quickacceff="lhc15.root")
{
    TFile* f = TFile::Open(fullacceff);
    TH1* hfullacceff = static_cast<TH1*>(f->Get("Eff"));
    hfullacceff->SetDirectory(0);
    delete f;

    f = TFile::Open(quickacceff);

    TGraph* g = static_cast<TGraph*>(f->Get("acceffdrop_ped_hv_lv_occ_config"));

    delete f;

    TAxis* x = hfullacceff->GetXaxis();

    TH1* hquick = static_cast<TH1*>(hfullacceff->Clone("hquick"));

    for ( Int_t i = 1; i <= x->GetNbins(); ++i )
    {
        TString srn = x->GetBinLabel(i);
        Int_t rn = static_cast<Int_t>(g->GetX()[i]);
        assert(rn=srn.Atoi());
        double quick = 0.2*(1 - g->GetY()[i-1]/100.0);
        // 0.2 is arbitrary scaling due to error
        // in original simulation that used 
        // a cut on child in AliGenParam...
        std::cout << x->GetBinLabel(i) << " " 
            << hfullacceff->GetBinContent(i)
            << " " 
            << quick << std::endl;
        hquick->SetBinContent(i,quick);
    }

    hfullacceff->Draw("histe");
    hquick->SetLineColor(2);
    hquick->Draw("histsame");
}
