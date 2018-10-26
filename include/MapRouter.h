#ifndef MAPROUTER_H
#define MAPROUTER_H

#include <memory>
#include <vector>
#include <istream>

class CMapRouter{
    public:
        using TNodeID = unsigned long;
    private:
        class CImplementation;
        std::unique_ptr< CImplementation > DData;

    public:
        CMapRouter();
        ~CMapRouter();
        
        static double HaversineDistance(double lat1, double lon1, double lat2, double lon2);        
        
        void LoadMap(std::istream &is);
        void GetMapExtents(double &minlat, double &minlon, double &maxlat, double &maxlon) const;
        TNodeID FindClosestNode(double lat, double lon);
        double FindShortestPath(TNodeID src, TNodeID dest, std::vector< TNodeID > &path);
        double FindFastestPath(TNodeID src, TNodeID dest, std::vector< TNodeID > &path);
        bool GetPathStreetNames(const std::vector< TNodeID > &path, std::vector< std::string > &streetnames) const;
};

#endif
