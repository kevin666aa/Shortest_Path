#ifndef XMLPARSER_H
#define XMLPARSER_H

#include <string>
#include <vector>
#include <istream>
#include <expat.h>

class CXMLParser{
    public:
        using TAttribute = struct   {
                                        std::string DKey;
                                        std::string DValue;
                                    };
    private:
        XML_Parser DExpatParser;
        static void XMLCALL StaticStartElement(void *data, const XML_Char *name, const XML_Char **attrs);
        static void XMLCALL StaticEndElement(void *data, const XML_Char *name);
        
    public:
        CXMLParser();
        virtual ~CXMLParser();
        
        void Parse(std::istream &is);
        virtual void StartElement(const std::string &name, const std::vector< TAttribute > &attrs) = 0;
        virtual void EndElement(const std::string &name) = 0;

};

#endif
