#include "MapRouter.h"
#include "XMLParser.h"
#include <cmath>
#include <limits>
#include <unordered_map>
#include <stdlib.h>
#include <iostream>
#include <tuple>
#include <iomanip>
#include <queue>
#define DEGREES_TO_RADIANS(angle)   (M_PI * (angle) / 180.0)
//retrived from https://stackoverflow.com/questions/20834838/using-tuple-in-unordered-map
    namespace std{
        namespace{
            template <class T>
            inline void hash_combine(std::size_t& seed, T const& v)
            {
                seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            }
            
            // Recursive template code derived from Matthieu M.
            template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
            struct HashValueImpl
            {
                static void apply(size_t& seed, Tuple const& tuple)
                {
                    HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
                    hash_combine(seed, get<Index>(tuple));
                }
            };
            
            template <class Tuple>
            struct HashValueImpl<Tuple,0>
            {
                static void apply(size_t& seed, Tuple const& tuple)
                {
                    hash_combine(seed, get<0>(tuple));
                }
            };
        }
        
        template <typename ... TT>
        struct hash<std::tuple<TT...>>
        {
            size_t
            operator()(std::tuple<TT...> const& tt) const
            {
                size_t seed = 0;
                HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
                return seed;
            }
        };
    };
    
    
    class CMapRouter::CImplementation : public CXMLParser{
    public:
        static double MaxDeltaLongitude(double dist, double lat){
            const double EarthCircumferenceMiles = 3959.88 * M_PI * 2.0;
            return 360.0 * dist / (cos(DEGREES_TO_RADIANS(lat)) * EarthCircumferenceMiles);
        };
        
        static double MaxDeltaLatitude(double dist){
            const double EarthCircumferenceMiles = 3959.88 * M_PI * 2.0;
            return 360.0 * dist / EarthCircumferenceMiles;
        };
        
        virtual void StartElement(const std::string &name, const std::vector< TAttribute > &attrs){
            
            /* if another map, 1000000*/
            if (name == "node") {
                the_name_before = name;
                TNodeID temp_id = atof(attrs[0].DValue.c_str());
                double la = atof(attrs[1].DValue.c_str());
                vertices.push_back(Vertex(la, atof(attrs[2].DValue.c_str()), temp_id));
                original.push_back( std::make_tuple(count_v, (la - 38.5) * 10000000) );
                IDToIndex.insert({temp_id, count_v});
                count_v++;
                
            } else if(name == "tag" && the_name_before == "way"){
                if (attrs[0].DValue == "name") {
                    streets[way_id] = attrs[1].DValue;
                } else if (attrs[0].DValue == "oneway" && attrs[1].DValue == "yes") {
                    oneway = true;
                } else if (attrs[0].DValue == "maxspeed") {
                    speed = atof(attrs[1].DValue.c_str());
                } else if(attrs[0].DValue == "ref"){
                    if (streets[way_id] == "") {
                        streets[way_id] = attrs[1].DValue;
                    }
                }
            } else if(name == "nd"){
                nd_ref_ids.push_back( atof(attrs[0].DValue.c_str()) );
            } else if (name == "way"){
                the_name_before = name;
                way_id = atof(attrs[0].DValue.c_str());
                
                //clear data from last way
                speed = 25;
                oneway = false;
                nd_ref_ids.clear();
            }
        }
        
        virtual void EndElement(const std::string &name){
            //   std::cout << "end: " << name << " "<< count_v << " ";
            
             if (name == "way") {
                double dist = 0;
            if (oneway) {

                for (int i = 0; i < (int)nd_ref_ids.size() - 1; i++) {
                    if (IDToIndex.find(nd_ref_ids[i]) != IDToIndex.end() && IDToIndex.find(nd_ref_ids[i + 1]) != IDToIndex.end()) {
                        dist =  HaversineDistance(vertices[IDToIndex[nd_ref_ids[i]] ].latitude, vertices[ IDToIndex[nd_ref_ids[i]] ].longtitude, vertices[ IDToIndex[nd_ref_ids[i + 1]] ].latitude, vertices[ IDToIndex[nd_ref_ids[i + 1]] ].longtitude);
                        
                        vertices[ IDToIndex[nd_ref_ids[i]]  ].adjacentlist.push_back( edge_to( IDToIndex[nd_ref_ids[i + 1]],  dist, dist / speed ) );
                         nodes_to_street_id[ std::make_tuple(nd_ref_ids[i], nd_ref_ids[i+ 1]) ] = way_id;
                    }
                }
            } else {
                for (int i = 0; i < (int)nd_ref_ids.size() - 1; i++) {
                    if (IDToIndex.find(nd_ref_ids[i]) != IDToIndex.end() && IDToIndex.find(nd_ref_ids[i + 1]) != IDToIndex.end()) {
                        dist =  HaversineDistance(vertices[IDToIndex[nd_ref_ids[i]] ].latitude, vertices[ IDToIndex[nd_ref_ids[i]] ].longtitude, vertices[ IDToIndex[nd_ref_ids[i + 1]] ].latitude, vertices[ IDToIndex[nd_ref_ids[i + 1]] ].longtitude);
                        
                        vertices[IDToIndex[nd_ref_ids[i]]  ].adjacentlist.push_back( edge_to( IDToIndex[ nd_ref_ids[i + 1] ], dist, dist / speed ) );
                        vertices[ IDToIndex[nd_ref_ids[i + 1]]  ].adjacentlist.push_back( edge_to( IDToIndex[ nd_ref_ids[i] ], dist, dist / speed ) );
                        nodes_to_street_id[ std::make_tuple(nd_ref_ids[i], nd_ref_ids[i+ 1]) ] = way_id;
                        nodes_to_street_id[ std::make_tuple(nd_ref_ids[i + 1], nd_ref_ids[i]) ] = way_id;
                    }
                }
            }
        }
          /*
            if (name == "way") {
                double dist = 0;
                for (int i = 0; i < (int)nd_ref_ids.size(); i++) {
                    if (IDToIndex.find(nd_ref_ids[i]) == IDToIndex.end()){
                        nd_ref_ids.erase(nd_ref_ids.cbegin() + i);
                    }
                }
                if (oneway) {
                    for (int i = 0; i < (int)nd_ref_ids.size() - 1; i++) {
                        dist =  HaversineDistance(vertices[IDToIndex[nd_ref_ids[i]] ].latitude, vertices[ IDToIndex[nd_ref_ids[i]] ].longtitude, vertices[ IDToIndex[nd_ref_ids[i + 1]] ].latitude, vertices[ IDToIndex[nd_ref_ids[i + 1]] ].longtitude);
                        vertices[ IDToIndex[nd_ref_ids[i]]  ].adjacentlist.push_back( edge_to( IDToIndex[nd_ref_ids[i + 1]],  dist, dist / speed ) );
                        nodes_to_street_id[ std::make_tuple(nd_ref_ids[i], nd_ref_ids[i+ 1]) ] = way_id;
                    }
                } else {
                    for (int i = 0; i < (int)nd_ref_ids.size() - 1; i++) {
                        dist =  HaversineDistance(vertices[IDToIndex[nd_ref_ids[i]] ].latitude, vertices[ IDToIndex[nd_ref_ids[i]] ].longtitude, vertices[ IDToIndex[nd_ref_ids[i + 1]] ].latitude, vertices[ IDToIndex[nd_ref_ids[i + 1]] ].longtitude);
                        
                        vertices[IDToIndex[nd_ref_ids[i]]  ].adjacentlist.push_back( edge_to( IDToIndex[ nd_ref_ids[i + 1] ], dist, dist / speed ) );
                        vertices[ IDToIndex[nd_ref_ids[i + 1]]  ].adjacentlist.push_back( edge_to( IDToIndex[ nd_ref_ids[i] ],  dist, dist / speed ) );
                        nodes_to_street_id[ std::make_tuple(nd_ref_ids[i], nd_ref_ids[i+ 1]) ] = way_id;
                        nodes_to_street_id[ std::make_tuple(nd_ref_ids[i + 1], nd_ref_ids[i]) ] = way_id;
                    }
                }
                
            }
           */
            
        }
        
        void radix_sort(){
            int max = 0;
            for ( int i =0 ; i< (int)original.size();i++){
                if (max < std::get<1>(original[i])){
                    max = std::get<1>(original[i]);}
            }
            std::vector< std::vector<std::tuple<int, int>> > radix_arr(10);
            for (auto i : original) {
                radix_arr.at(std::get<1>(i)%10).push_back(i);
            }
            
            int exp = 100;
            for (; max/exp > 0; exp *= 10) {
                auto copy(radix_arr);
                radix_arr.clear(); radix_arr.resize(10);
                for (auto vec : copy) {
                    for (auto val: vec) {
                        radix_arr.at( (std::get<1>(val) % exp) / (exp/10) ).push_back(val);
                    }
                }
            }
            exp /= 10;
            auto copy(radix_arr);
            radix_arr.clear(); radix_arr.resize(10);
            for (auto vec : copy) {
                for (auto val: vec) {
                    radix_arr.at( std::get<1>(val)/exp ).push_back(val);
                }
            }
            original.clear();
            for (auto vec : radix_arr) {
                for (auto val: vec) {
                    sorted.push_back(std::make_tuple(std::get<0>(val), (double)std::get<1>(val) / 10000000 + 38.5) );
                }
            }
        } // to sort the "sorted" vector of tuples
        
        // for loading the map
        int speed;
        TNodeID way_id;
        bool oneway;
        std::string the_name_before;
        std::vector<TNodeID> nd_ref_ids;
        
        class edge_to{
        public:
            int index;
            double distance;
            double time;
            edge_to(int next_id, double distance, double time) : index(next_id), distance(distance), time(time) { }
            edge_to(){ index = -1; distance = 0; time = 0;}
            double get_weight(bool dist){
                if (dist) {  return distance;  }
                return time;
            }
        };
        
        class Vertex{
        public:
            double latitude;
            double longtitude;
            TNodeID id;
            std::vector<edge_to> adjacentlist;
            Vertex(double latitude, double longtitude, TNodeID id) : latitude(latitude), longtitude(longtitude), id(id) {
            }
            Vertex(){ latitude = 0; longtitude = 0; id = -1;}
        };
        
        //  all the data stored
        int count_v = 0;
        std::vector<Vertex> vertices; // the core list
        std::vector< std::tuple<int, int> > original;  // a partial copy of vertices, index and latitude * 10000000
        std::unordered_map<TNodeID, std::string> streets; // convert a way_id into a name
        std::unordered_map<std::tuple<TNodeID, TNodeID >, TNodeID> nodes_to_street_id;
        std::unordered_map<TNodeID, int> IDToIndex;  // convert a node id to index of vertices, used in loading way
        
        std::vector< std::tuple<int, double> > sorted; //already sorted, can do binary search for input nodes
        
        class CompareDist{
        public:
            bool operator()(std::tuple<int,double> n1, std::tuple<int,double> n2) {
                return std::get<1>(n1)> std::get<1>(n2);
            }
        };
            //retrieved from https://stackoverflow.com/questions/12685787/pair-inside-priority-queue
        std::priority_queue< std::tuple<int, double>, std::vector<std::tuple<int, double>>, CompareDist> prQ;
        
        double Dijkstra_SSSP(const TNodeID& src, const TNodeID& dest, std::vector< TNodeID > &path, bool Tim_Dis){
            int S = IDToIndex[src], D = IDToIndex[dest];
            std::vector<bool> visited(count_v, false);
            std::vector<int> prev_i(count_v, -1);
            std::vector<double> weight(count_v, std::numeric_limits<double>::max());
            prev_i[S] = S;
            weight[S] = 0;
            visited[S] = true;
            prQ.push(std::make_tuple(S, 0));
            int curr = -1;
            double altweight = 0;
            //initialize data;
            while (!prQ.empty()) {
                curr = std::get<0>(prQ.top());
                prQ.pop();
                if (curr == D &&  prev_i[curr] != -1) { break;  }
                for (auto other : vertices[curr].adjacentlist) {
                    double Wei = other.get_weight(Tim_Dis);
                    altweight = weight[curr] + Wei;
                    if (altweight < weight[other.index]){
                        weight[other.index] = altweight;
                        prev_i[other.index] = curr;
                        if (!visited[other.index]) {
                            visited[other.index] = true;
                            prQ.push(std::make_tuple(other.index, altweight));
                        }
                    }
                }
            }
            double finalweight = weight[D];
            if (prev_i[D] == -1) {
                return 0;
            }
            //load in the shortest path
            
            std::vector< TNodeID > temp;
            while (curr != S) {
                temp.push_back(curr);
                curr = prev_i[curr];
            }
            
            temp.push_back(S);
            for (auto i = (int)temp.size() - 1; i >= 0; i--) {
                path.push_back(vertices[temp[i]].id );
            }  //get path down
            return finalweight;
        }
        
    };
    
    CMapRouter::CMapRouter() : DData(std::make_unique< CImplementation>()){}
    CMapRouter::~CMapRouter(){}
    
    // Modified from https://rosettacode.org/wiki/Haversine_formula#C.2B.2B
    double CMapRouter::HaversineDistance(double lat1, double lon1, double lat2, double lon2){
        double LatRad1 = DEGREES_TO_RADIANS(lat1);
        double LatRad2 = DEGREES_TO_RADIANS(lat2);
        double LonRad1 = DEGREES_TO_RADIANS(lon1);
        double LonRad2 = DEGREES_TO_RADIANS(lon2);
        double DeltaLat = LatRad2 - LatRad1;
        double DeltaLon = LonRad2 - LonRad1;
        double DeltaLatSin = sin(DeltaLat/2);
        double DeltaLonSin = sin(DeltaLon/2);
        double Computation = asin(sqrt(DeltaLatSin * DeltaLatSin + cos(LatRad1) * cos(LatRad2) * DeltaLonSin * DeltaLonSin));
        const double EarthRadiusMiles = 3959.88;
        
        return 2 * EarthRadiusMiles * Computation;
    }
    
    void CMapRouter::GetMapExtents(double &minlat, double &minlon, double &maxlat, double &maxlon) const{
        minlat = 1000; minlon = 1000; maxlat = -1000; maxlon = -1000;
        for (const auto& i : DData->vertices) {
            if (maxlat < i.latitude) { maxlat = i.latitude;}
            else if (minlat > i.latitude){ minlat = i.latitude;  }
            
            if (maxlon < i.longtitude) { maxlon = i.longtitude; }
            else if (minlon > i.longtitude){   minlon = i.longtitude;}
        }
        DData->radix_sort();
    }
    
    void CMapRouter::LoadMap(std::istream &is){
        DData->CXMLParser::Parse(is);
        //   std:: cout << DData->IDToIndex[2637746337] << " " << DData->IDToIndex[62208369] << "   a a a a a                 ";
        /*for (auto k : DData->vertices) {
         std:: cout << "<node id=" << std::setprecision(10) <<k.id << " lat=" << k.latitude << " lon=" << k.longtitude << "/> " << std::endl;
         }*/
    }
    
    CMapRouter::TNodeID CMapRouter::FindClosestNode(double lat, double lon){
        /* CMapRouter::TNodeID min_index = 0;
         //    double latitude =  DData->vertices[ std::get<0>(DData->sorted[low]) ].latitude;
         double min_dist = 100000;
         double cur_dist = 0;;
         for (auto k : DData->vertices) {
         cur_dist = HaversineDistance(lat, lon, k.latitude, k.longtitude);
         
         if (cur_dist < min_dist) {
         min_dist = cur_dist;
         min_index = k.id;
         }
         }
         // std::cout << min_index;
         return min_index;*/
        int low = 0, high = (int)DData->sorted.size() - 1, mid = 0;
        while (low <= high) {
            mid = low + (high - low) / 2;
            if ( lat < std::get<1>(DData->sorted.at(mid)) ) {
                high = mid - 1;
            } else if( lat > std::get<1>(DData->sorted.at(mid)) ){
                low = mid + 1;
            } else {
                break;
            }
        }
        
        int min_index = 0;
        double latitude =  DData->vertices[ std::get<0>(DData->sorted[low]) ].latitude;
        double min_dist = HaversineDistance(lat, lon, latitude, DData->vertices[std::get<0>(DData->sorted.at(low))].longtitude);
        double cur_dist;
        
        double maxdeltalatitude = DData->MaxDeltaLatitude(min_dist);
        
        for (int i = low; i < (int)DData->sorted.size() - 1; i++) {
            if (std::get<1>(DData->sorted[i]) - latitude > maxdeltalatitude) {
                break;
            }
            cur_dist = HaversineDistance(lat, lon, std::get<1>(DData->sorted[i]), DData->vertices[ std::get<0>(DData->sorted[i]) ].longtitude);
            if (min_dist > cur_dist) {
                min_dist = cur_dist;
                min_index = std::get<0>(DData->sorted[i]);
            }
        }
        
        for (int i = low; i >= 0; i--) {
            if (latitude - std::get<1>(DData->sorted[i]) > maxdeltalatitude) {
                break;
            }
            cur_dist = HaversineDistance(lat, lon, std::get<1>(DData->sorted[i]), DData->vertices[ std::get<0>(DData->sorted[i]) ].longtitude);
            if (min_dist > cur_dist) {
                min_dist = cur_dist;
                min_index = std::get<0>(DData->sorted[i]);
            }
        }
        return DData->vertices[min_index].id;
    }
    
    double CMapRouter::FindShortestPath(TNodeID src, TNodeID dest, std::vector< TNodeID > &path){
        return DData->Dijkstra_SSSP(src, dest, path, true);  //need change
    }
    
    double CMapRouter::FindFastestPath(TNodeID src, TNodeID dest, std::vector< TNodeID > &path){
        return DData->Dijkstra_SSSP(src, dest, path, false); //dif
    }
    
    bool CMapRouter::GetPathStreetNames(const std::vector< TNodeID > &path, std::vector< std::string > &streetnames) const{
        TNodeID last_one = 1;
        std::vector<TNodeID> passed_streets;
        for (int i = 0; i < (int)path.size() - 1; i++) {
            if (DData->nodes_to_street_id[std::make_tuple(path[i], path[i + 1])] != last_one) {
                last_one = DData->nodes_to_street_id[std::make_tuple(path[i], path[i + 1])];
                passed_streets.push_back(last_one);
            }
        }
        if (passed_streets.size() == 0) {
            return false;
        }
        std::string a;
        for (auto k : passed_streets) {
            if (DData->streets[k] != "" && DData->streets[k] != a) {
                a = DData->streets[k];
                streetnames.push_back(DData->streets[k]);
            }
        }
        return true;
    }

/*
for (auto i : DData->vertices) {
    std::cout << "id: " << i.id << " la: " << std::setprecision(9) << i.latitude << " lo: " << std::setprecision(9) << i.longtitude <<  " adj: " << std::endl;
    for (auto k : i.adjacentlist) {
        std::cout << " index: "<< k.index << " way_id: " << k.way_id << " wayname: " << DData->streets[k.way_id]<< " dist: " << k.distance << " time: " << k.time << std::endl;
    }
    std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << std::endl;

for (auto i : path) {
    std::cout << i << "" << std::endl;
    }
    for (auto i  : street_id) {
        std::cout << i << " : " << streets[i] << std::endl;
    }
*/
