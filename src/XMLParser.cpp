#include "XMLParser.h"
#include <iostream>

CXMLParser::CXMLParser(){
    DExpatParser = XML_ParserCreate(NULL);
    XML_SetUserData(DExpatParser, this);
    XML_SetElementHandler(DExpatParser, StaticStartElement, StaticEndElement);
}

CXMLParser::~CXMLParser(){
    XML_ParserFree(DExpatParser);
}

void XMLCALL CXMLParser::StaticStartElement(void *data, const XML_Char *name, const XML_Char **attrs){
    CXMLParser *Parser = static_cast< CXMLParser *>(data);
    std::vector< TAttribute > Attributes;
    
    for(int Index = 0; attrs[Index]; Index += 2){
        Attributes.push_back({attrs[Index], attrs[Index + 1]});    
    }
    Parser->StartElement(name, Attributes);
}

void XMLCALL CXMLParser::StaticEndElement(void *data, const XML_Char *name){
    CXMLParser *Parser = static_cast< CXMLParser *>(data);
    Parser->EndElement(name);
}

void CXMLParser::Parse(std::istream &is){
    const int BufferSize = 512;
    char Buffer[BufferSize];
    
    while(!is.eof()){
        is.read(Buffer,BufferSize);
        if(XML_STATUS_ERROR == XML_Parse(DExpatParser, Buffer, (int)is.gcount(), is.eof())){
            std::cerr<<XML_ErrorString(XML_GetErrorCode(DExpatParser))<<" @"<<XML_GetCurrentLineNumber(DExpatParser)<<std::endl;     
            return;
        }
    }
}

