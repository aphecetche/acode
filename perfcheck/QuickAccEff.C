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
#include <map>
#include <set>
#include <vector>

typedef std::pair<int,int> ManuPair;

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
    ClusterLocation(int detelemid=0, int b=0, int nb=0)
        : mDetElemId(detelemid), mBendingManu(b), mNonBendingManu(nb) {}

    int mDetElemId;
    int mBendingManu;
    int mNonBendingManu;

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
};

struct CompactEvent
{
    CompactEvent() : mTracks() {}
    std::vector<CompactTrack> mTracks;
};

struct CompactMapping
{
    CompactMapping() : mBendingDeIx(),
    mNonBendingDeIx(), mManuIds(), mDeToManuMap() {}

    // associate each detElemId to an index
    // in the mManuIds array (indicating the
    // first bending manuId of that DE)
    std::map<int,int> mBendingDeIx;

    // associate each detElemId to an index
    // in the mManuIds array (indicating the
    // first non bending manuId of that DE)
    std::map<int,int> mNonBendingDeIx;

    // array containing the (local) manuId
    // mManuIds[abs]=local
    std::vector<int> mManuIds;

    // associate each detElemId to a map of
    // manuId -> absmanuid
    std::map<int, std::map<int,int> > mDeToManuMap;

};

UInt_t MANUBADPEDMASK = ( 1 << 0 );
UInt_t MANUBADHVMASK = ( 1 << 1 );
UInt_t MANUBADLVMASK =  (1 << 2 );
UInt_t MANUBADOCCMASK = ( 1 << 3 );
UInt_t MANUOUTOFCONFIG = ( 1 << 4 );

std::string CauseAsString(UInt_t cause)
{
    std::string rv = "";

    if ( cause & MANUBADPEDMASK ) rv += " ped ";
    if ( cause & MANUBADHVMASK ) rv += " hv ";
    if ( cause & MANUBADLVMASK ) rv += " lv ";
    if ( cause & MANUBADOCCMASK ) rv += " occ ";
    if ( cause & MANUOUTOFCONFIG ) rv += " config ";

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

        AliMpDEIterator it;
        it.First();
        std::vector<int> deids;
        while (!it.IsDone())
        {
            Int_t detElemId = it.CurrentDEId();

            if ( AliMpDEManager::GetStationType(detElemId) != AliMp::kStationTrigger )
            {
                deids.push_back(detElemId);
            }
            it.Next();
        }

        std::sort(deids.begin(),deids.end());

        // then for each DE get one sorted list of manuIds
        int totalNofManus(0);

        for ( std::vector<int>::size_type i = 0; i < deids.size(); ++i )
        {
            AliMpManuIterator it;
            Int_t detElemId, manuId;
            std::vector<int> manuids;
            Int_t bix=-1, nbix=-1;

            while ( it.Next(detElemId,manuId) )
            {
                if ( detElemId == deids[i] )
                {
                    cm->mDeToManuMap[detElemId][manuId]=manuids.size();
                    manuids.push_back(manuId);
                    cm->mManuIds.push_back(manuId);
                    if ( manuId >= 1024 && bix < 0 ) 
                    {
                        bix = cm->mManuIds.size()-1;
                        cm->mBendingDeIx[detElemId] = bix;
                    }
                    if ( manuId < 1024 && nbix < 0 ) 
                    {
                        nbix = cm->mManuIds.size()-1;
                        cm->mNonBendingDeIx[detElemId] = nbix;
                    }
                }
                std::sort(manuids.begin(),manuids.end());
                totalNofManus += manuids.size();

                std::cout << Form("DE %04d (%3d) ",deids[i],
                        manuids.size());

                for ( std::vector<int>::size_type j = 0; 
                        j < manuids.size(); ++j )
                {
                    std::cout << Form("%04d ",manuids[j]);
                }
                std::cout << std::endl;
            }
            std::cout << "Total number of manus : " << totalNofManus << std::endl;


            // consistency check
            for ( unsigned int i = 0; i < cm->mManuIds.size(); ++i )
            {
                std::cout << Form("mManuIds[%4d]=%4d",i,cm->mManuIds[i])
                    << std::endl;
            }
        }
        return cm;
    }
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

Int_t FindManuAbsId(Int_t detElemId, Int_t manuId)
{
    CompactMapping* cm = GetCompactMapping();

    int bendingIx = cm->mBendingDeIx[detElemId];
    int nonBendingIx = cm->mNonBendingDeIx[detElemId];

    std::cout << Form("Bending %4d NonBendng %4d",bendingIx,nonBendingIx) << std::endl;
}

void GetClusterLocation(Int_t detElemId,
        Double_t xg, 
        Double_t yg, 
        Double_t zg,
        Int_t& bendingManuAbsId,
        Int_t& nonBendingManuAbsId)
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

    bendingManuAbsId = FindManuAbsId(de->GetId(),manuBending);
    nonBendingManuAbsId = FindManuAbsId(de->GetId(),manuNonBending);
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
            ClusterLocation cl(cluster->GetDetElemId(),b,nb);
            compactTrack.mClusters.push_back(cl);

            // std::cout << cluster->GetDetElemId() << " ("
            //     << b << "," << nb << ")";
        }

        compactEvent.mTracks.push_back(compactTrack);
    }
}

UInt_t GetTrackStationMask(const CompactTrack& track)
{
    UInt_t mask = 0;

    for ( std::vector<ClusterLocation>::size_type i = 0;
            i < track.mClusters.size(); ++i )
    {
        const ClusterLocation& cl = track.mClusters[i];
        Int_t stationId = AliMpDEManager::GetChamberId(cl.mDetElemId)/2;
        mask |= ( 1 << stationId);
    }
    return mask;
}

// Bool_t ValidateTrack(const CompactTrack& track,
//         const std::map<int,std::set<int> >& badManuList)
// {
//     /// We remove from the track all the clusters 
//     /// located on a bad manu.
//     /// Then we consider the track survived if we get
//     /// at least one cluster per station
//
//     UInt_t originalMask = 0; 
//     UInt_t filteredMask = 0;
//
//     for ( std::vector<ClusterLocation>::size_type i = 0;
//             i < track.mClusters.size(); ++i )
//     {
//         const ClusterLocation& cl = track.mClusters[i];
//         std::map<int,std::set<int> >::const_iterator it = badManuList.find(cl.mDetElemId);
//
//         Int_t stationId = AliMpDEManager::GetChamberId(cl.mDetElemId)/2;
//         originalMask |= ( 1 << stationId);
//         if  ( it != badManuList.end() ) 
//         {
//             const std::set<int>& badones = it->second;
//             if ( std::find(badones.begin(),badones.end(),cl.mBendingManu) != badones.end() || std::find(badones.begin(),badones.end(),cl.mNonBendingManu) != badones.end())
//             {
//                 continue; 
//             }
//         }
//         filteredMask |= ( 1 << stationId );
//     }
//
//     if ( filteredMask != originalMask )
//     {
//         return kFALSE;
//     }
//     return kTRUE;
// }

Bool_t ValidateTrack(const CompactTrack& track,
        const std::vector<UInt_t>& manuStatus,
        UInt_t causeMask)
{

    /// We remove from the track all the clusters 
    /// located on a bad manu.
    /// Then we consider the track survived if we get
    /// at least one cluster per station

    if ( manuStatus.empty() || causeMask == 0 ) return kTRUE;

    UInt_t originalMask = 0; 
    UInt_t filteredMask = 0;

    for ( std::vector<ClusterLocation>::size_type i = 0;
            i < track.mClusters.size(); ++i )
    {
        const ClusterLocation& cl = track.mClusters[i];

        UInt_t bendingMask = manuStatus[cl.mBendingManu];
        UInt_t nonBendingMask = manuStatus[cl.mNonBendingManu];

        Int_t stationId = AliMpDEManager::GetChamberId(cl.mDetElemId)/2;
        originalMask |= ( 1 << stationId);

        if ( ( bendingMask & causeMask == causeMask ) ||
                ( nonBendingMask & causeMask == causeMask ) )
        {
            continue; 
        }
        filteredMask |= ( 1 << stationId );
    }

    if ( filteredMask != originalMask )
    {
        return kFALSE;
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

    for ( std::vector<CompactEvent>::size_type i = 0;
            i < events.size(); ++i )
    {
        const CompactEvent& e = events[i];


        for ( std::vector<CompactTrack>::size_type j = 0;
                j < e.mTracks.size(); ++j ) 
        {
            const CompactTrack& t1 = e.mTracks[j];

            if (!ValidateTrack(t1,manustatus,causeMask)) continue;

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

    return h;
}

Int_t ConvertESD(const char* inputfile,
        const char* outputfile)
{
    AliCDBManager::Instance()->SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2016/OCDB");

    AliCDBManager::Instance()->SetRun(0);
    AliMpCDB::LoadAll();

    AliGeomManager::LoadGeometry("geometry.root");

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

void GetManuStatus(AliMUONPadStatusMaker& statusMaker,
        std::vector<UInt_t>& manustatus)
{
    manustatus.clear();

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
                manuStatus |= MANUOUTOFCONFIG;
            }

            Int_t manuAbsId = FindManuAbsId(detElemId,manuId);
            manustatus[manuAbsId] = manuStatus;
        }
    }

    assert(ntotal==1064008);
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

    causes.push_back(MANUBADPEDMASK);
    causes.push_back(MANUBADPEDMASK | MANUBADOCCMASK);

    for ( std::vector<UInt_t>::size_type i = 0; i < causes.size(); ++i )
    {
        TGraph* g = new TGraph(vrunlist.size());
        gdrop.push_back(g);
        g->SetName(Form("acceffdrop_%s",CauseAsString(causes[i]).c_str()));
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.5);
    }

    AliCDBManager* man = AliCDBManager::Instance();
    man->SetDefaultStorage(ocdbpath);

    for ( std::vector<int>::size_type i = 0; i < vrunlist.size(); ++i )
    {
        Int_t runNumber = vrunlist[i];
        man->SetRun(runNumber);

        AliMpCDB::LoadAll();
        AliMUONCalibrationData cd(runNumber,true);

        AliMUONPadStatusMaker statusMaker(cd);

        AliMUONRecoParam* recoParam = AliMUONCDB::LoadRecoParam();

        statusMaker.SetLimits(*recoParam);

        std::vector<UInt_t> manustatus;

        GetManuStatus(statusMaker,manustatus);

        for ( std::vector<UInt_t>::size_type icause = 0; icause < causes.size(); ++icause )
        {
            TH1* h = ComputeMinv(events,manustatus,causes[icause]);
            h->SetName(Form("hminv%6d",runNumber));
            Double_t drop = h->Integral(b1,b2)/referenceNofJpsi;
            std::cout << Form("RUN %6d %20s AccxEff drop %7.2f",
                    runNumber,CauseAsString( causes[icause] ).c_str(),drop) << std::endl;
            gdrop[icause]->SetPoint(i,runNumber,drop);
            delete h;
        }

        AliCDBManager::Instance()->ClearCache();
    }


    TFile* fout = TFile::Open(outputfile,"recreate");
    for ( std::vector<UInt_t>::size_type icause = 0; icause < causes.size(); ++icause )
    {
        gdrop[icause]->Write();
    }
    delete fout;
}
