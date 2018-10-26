#include <gtest/gtest.h>
#include "MapRouter.h"
#include <sstream>
#include <iostream>

TEST(MapRouter, BasicTest){

    std::stringstream OSMStream("<?xml version='1.0' encoding='UTF-8'?>"
                                "<osm version=\"0.6\" generator=\"osmconvert 0.8.5\">"
                                "  <node id=\"0\" lat=\"38.5250807\" lon=\"-121.7626513\"/>"
                                "  <node id=\"1\" lat=\"38.5253341\" lon=\"-121.7577639\"/>"
                                "  <node id=\"2\" lat=\"38.5253341\" lon=\"-121.7603067\"/>"
                                "  <node id=\"3\" lat=\"38.5253341\" lon=\"-121.7606424\"/>"
                               /* "<node id=\"62224345\" lat=\"38.5250807\" lon=\"-121.7632631\"/> "
                                "<node id=\"62224347\" lat=\"38.5253341\" lon=\"-121.7626513\"/>"
                                "<node id=\"62224348\" lat=\"38.525476\" lon=\"-121.762391\"/>"
                                "<node id=\"62224351\" lat=\"38.5253341\" lon=\"-121.7618937\"/>"
                                "<node id=\"62224354\" lat=\"38.5260106\" lon=\"-121.7616517\"/>"
                                "<node id=\"62224356\" lat=\"38.5267845\" lon=\"-121.7606424\"/>"
                                "<node id=\"62224358\" lat=\"38.5270246\" lon=\"-121.7603067\"/>"
                                "<node id=\"62224360\" lat=\"38.5275044\" lon=\"-121.7595218\"/>"
                                "<node id=\"62224362\" lat=\"38.5276462\" lon=\"-121.7592564\"/>"
                                "<node id=\"62224365\" lat=\"38.5253341\" lon=\"-121.7588672\"/>"
                                "<node id=\"62224369\" lat=\"38.5284058\" lon=\"-121.7577639\"/>"
                                "<node id=\"62224641\" lat=\"38.5222679\" lon=\"-121.7673869\">"*/
                                "  <way id=\"0\">"
                                "    <nd ref=\"0\"/>"
                                "    <nd ref=\"1\"/>"
                                "    <nd ref=\"3\"/>"
                                "    <nd ref=\"2\"/>"
                                "    <tag k=\"highway\" v=\"residential\"/>"
                                "    <tag k=\"name\" v=\"Main Street\"/>"
                                "  </way>"
                                "  <way id=\"1\">"
                                "    <nd ref=\"0\"/>"
                                "    <nd ref=\"2\"/>"
                                "    <tag k=\"oneway\" v=\"yes\"/>"
                                "    <tag k=\"maxspeed\" v=\"15 mph\"/>"
                                "    <tag k=\"highway\" v=\"residential\"/>"
                                "    <tag k=\"name\" v=\"Shortcut Way\"/>"
                                "  </way>"    
                                "</osm>");
    CMapRouter Router;
    std::vector< CMapRouter::TNodeID > Path;
    std::vector< std::string > Streets;
    /*
    Router.LoadMap(OSMStream);
    EXPECT_EQ(Router.FindClosestNode(38.5253341,-121.7597639), 10);
//    EXPECT_EQ(2,5);
    EXPECT_EQ(Router.FindClosestNode(38.5222679,-121.7673869), 10);
//    EXPECT_EQ(2,5);
    EXPECT_EQ(Router.FindShortestPath(0,2,Path), 54.461481057457596933);
    EXPECT_EQ(Path.size(), 2);
    if(Path.size() == 2){
        EXPECT_EQ(Path[0],0);
        EXPECT_EQ(Path[1],2);
    }
    Path.clear();
    EXPECT_EQ(Router.FindFastestPath(0,2,Path), 3.6307654038305066102);
    EXPECT_EQ(Path.size(), 2);
    if(Path.size() == 2){
        EXPECT_EQ(Path[0],0);
        EXPECT_EQ(Path[1],2);
    }
    Path.clear();
    EXPECT_TRUE(Router.GetPathStreetNames(Path,Streets));
    EXPECT_EQ(Streets.size(), 1);
    if(Path.size() == 1){
        EXPECT_EQ(Streets[0], "Shortcut Way");
    }
//    EXPECT_EQ(2,3);
    Path.clear();
    EXPECT_EQ(Router.FindShortestPath(2,0,Path), 191.93646327377945227);
    EXPECT_EQ(Path.size(), 4);
    if(Path.size() == 4){
        EXPECT_EQ(Path[0],2);
        EXPECT_EQ(Path[1],3);
        EXPECT_EQ(Path[2],1);
        EXPECT_EQ(Path[3],0);
    }
    Streets.clear();
    EXPECT_TRUE(Router.GetPathStreetNames(Path,Streets));
    EXPECT_EQ(Streets.size(), 1);
    if(Path.size() == 1){
        EXPECT_EQ(Streets[0], "Main Street");
    }*/
//    double minlat, minlon, maxlat, maxlon;
//    Router.GetMapExtents(minlat, minlon, maxlat, maxlon);
//    EXPECT_EQ(1,2);
}
