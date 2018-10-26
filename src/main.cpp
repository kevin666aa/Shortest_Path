#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <exception>
#include <vector>
#include <cmath>
#include <iomanip>
#include "MapRouter.h"


int main(int argc, char *argv[]){
    int NumRoutes = 1000;
    
    if(2 > argc){
        std::cerr<<"Syntax Error: proj4 file [numpoints]"<<std::endl;
        return EXIT_FAILURE;
    }
    if(3 <= argc){
        NumRoutes = std::stoul(argv[2]);
    }
    std::cerr<<"Loading"<<std::endl;
    
    auto LoadStart = std::chrono::steady_clock::now();
    CMapRouter Router;

    std::ifstream InFile(argv[1]);
    Router.LoadMap(InFile);
    auto LoadDuration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-LoadStart);
    double MinLat, MinLon, MaxLat, MaxLon;
    
    std::vector< std::tuple<double, double, double, double> > RandomSourceDest;
    std::vector< std::vector< CMapRouter::TNodeID > > ShortestPaths;
    std::vector< double > ShortestDistance;
    std::vector< std::vector< CMapRouter::TNodeID > > FastestPaths;
    std::vector< double > FastestTime;

    std::cerr<<"Loaded"<<std::endl;
    std::cerr<<"Generating End Points"<<std::endl;
    Router.GetMapExtents(MinLat, MinLon, MaxLat, MaxLon);
    //std::cout << MinLat << " "<< MinLon << " - " << MaxLat << " "<< MaxLon << std::endl;
   //  std::cout << std::endl;
    srand(1);
    for(int Index = 0; Index < NumRoutes; Index++){
        double SLat = MinLat + (MaxLat - MinLat) * rand() / RAND_MAX;
        double SLon = MinLon + (MaxLon - MinLon) * rand() / RAND_MAX;
        double DLat = MinLat + (MaxLat - MinLat) * rand() / RAND_MAX;
        double DLon = MinLon + (MaxLon - MinLon) * rand() / RAND_MAX;
        RandomSourceDest.push_back(std::make_tuple(SLat, SLon, DLat, DLon));
       // std::cout << SLat << " " << SLon << " - " << DLat << " " << DLon << std::endl;
    }
   
    ShortestPaths.resize(NumRoutes);
    ShortestDistance.resize(NumRoutes);
    FastestPaths.resize(NumRoutes);
    FastestTime.resize(NumRoutes);
    //std::cout << std::endl;
    std::cerr<<"Finding Paths"<<std::endl;    
    auto ProcessingStart = std::chrono::steady_clock::now();
    for(int Index = 0; Index < NumRoutes; Index++){
        auto CurSourceDest = RandomSourceDest[Index];
        auto SourceID = Router.FindClosestNode(std::get<0>(CurSourceDest),std::get<1>(CurSourceDest));
        auto DestID = Router.FindClosestNode(std::get<2>(CurSourceDest),std::get<3>(CurSourceDest));
       // std::cout << " "<<SourceID << " " << DestID <<std::endl;
        ShortestDistance[Index] = Router.FindShortestPath(SourceID, DestID, ShortestPaths[Index]);
        FastestTime[Index] = Router.FindFastestPath(SourceID, DestID, FastestPaths[Index]);
    }
    
    auto ProcessingDuration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-ProcessingStart);
    std::cerr<<"Outputting Results"<<std::endl;    
    
    for(int Index = 0; Index < NumRoutes; Index++){
        std::vector< std::string > StreetNames;
        Router.GetPathStreetNames(ShortestPaths[Index], StreetNames);
        std::cout<<Index<<" SN:";
        bool First = true;
        for(auto Street : StreetNames){
            if(!First){
                std::cout<<",";   
            }
            std::cout<<" "<<Street;
            First = false;
        }
        std::cout<<std::endl;
        std::cout<<Index<<" SP:";
        if((std::numeric_limits<double>::max() == ShortestDistance[Index]) or ShortestPaths[Index].empty()){
            std::cout<<" (N/A)";   
        }
        else{
            std::cout<<" ("<<std::setprecision(2)<<ShortestDistance[Index]<<" mi)";    
        }
        for(auto Node : ShortestPaths[Index]){
            std::cout<<" "<<Node;
        }
        std::cout<<std::endl;
        std::cout<<Index<<" FP:";
        if((std::numeric_limits<double>::max() == FastestTime[Index]) or FastestPaths[Index].empty()){
            std::cout<<" (N/A)";   
        }
        else{
            int IntegralTime;
            double Remainder = FastestTime[Index];
            bool PrecededOutput = false;
            std::cout<<" (";
            IntegralTime = Remainder;
            if(IntegralTime){
                std::cout<<IntegralTime<<":";
                PrecededOutput = true;
            }
            Remainder -= IntegralTime;
            Remainder *= 60;
            IntegralTime = Remainder;
            if(IntegralTime or PrecededOutput){
                if(PrecededOutput){
                    std::cout<<std::setw(2)<<std::setfill('0');
                }
                std::cout<<IntegralTime<<":";
                PrecededOutput = true;
            }
            Remainder -= IntegralTime;
            Remainder *= 60;
            IntegralTime = Remainder;
            if(PrecededOutput){
                std::cout<<std::setw(2)<<std::setfill('0');
            }
            std::cout<<IntegralTime;
            if(PrecededOutput){
                std::cout<<")";
            }
            else{
                Remainder -= IntegralTime;
                Remainder *= 100;
                IntegralTime = Remainder;
                std::cout<<"."<<std::setw(2)<<std::setfill('0')<<IntegralTime<<" s)";    
            }
        }
        for(auto Node : FastestPaths[Index]){
            std::cout<<" "<<Node;
        }
        std::cout<<std::endl;
    }
    
    std::cout<<"Duration (load): "<<LoadDuration.count()<<std::endl;
    std::cerr<<"Duration (load): "<<LoadDuration.count()<<std::endl;
    std::cout<<"Duration (proc): "<<ProcessingDuration.count()<<std::endl;
    std::cerr<<"Duration (proc): "<<ProcessingDuration.count()<<std::endl;
    auto SamplesPerDay = long((86400000 - LoadDuration.count()) / (double(ProcessingDuration.count()) / NumRoutes));
    auto MarginOfError = long(double(SamplesPerDay) / sqrt(NumRoutes));
    std::cout<<"Queries per day: "<<SamplesPerDay<<" (+-"<<MarginOfError<<"), "<<(SamplesPerDay - MarginOfError)<<" min"<<std::endl;
    std::cerr<<"Queries per day: "<<SamplesPerDay<<" (+-"<<MarginOfError<<"), "<<(SamplesPerDay - MarginOfError)<<" min"<<std::endl;
    
    
    return EXIT_SUCCESS;
}
