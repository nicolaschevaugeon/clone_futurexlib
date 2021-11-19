/* 
    This file is a part of eXlibris C++ Library
    under the GNU Lesser General Public License.
    See the NOTICE.md & LICENSE.md files for terms 
    and conditions.
*/

#include "xParseData.h"
#include <sstream>

using std::cout;
using std::endl;

void xParseData::registerInt(std::string name, int defaultval){
    intdata.insert(std::make_pair(name, defaultval));
    isdefined.insert(std::make_pair(name,0));
}

void xParseData::registerDouble(std::string name, double defaultval){
    doubledata.insert(std::make_pair(name, defaultval));
    isdefined.insert(std::make_pair(name,0));
}

void xParseData::registerString(std::string name, std::string defaultval){
    stringdata.insert(std::make_pair(name,defaultval));
    isdefined.insert(std::make_pair(name,0));
}

void  xParseData::registerVector(std::string name, xtensor::xVector<> defaultval){
    xvectordata.insert(std::make_pair(name, defaultval));
    isdefined.insert(std::make_pair(name,0));
}

void  xParseData::registerListString(std::string name, std::list<std::string > defaultval){
    liststringdata.insert(std::make_pair(name, defaultval));
    isdefined.insert(std::make_pair(name,0));
}

void xParseData::registerListVector(std::string name, std::list<xtensor::xVector<> > defaultval){
    listvectordata.insert(std::make_pair(name, defaultval));
    isdefined.insert(std::make_pair(name,0));
}

void  xParseData::registerTensor2(std::string name, xtensor::xTensor2<> defaultval){
    xtensor2data.insert(std::make_pair(name, defaultval));
    isdefined.insert(std::make_pair(name,0));
}

void  xParseData::registerMapStringDouble(std::string name, std::map<std::string, double> defaultval){//GREG
    mapstringdoubledata.insert(std::make_pair(name, defaultval));
    isdefined.insert(std::make_pair(name,0));
}

double xParseData::getDouble(std::string name) const{
    std::map<std::string,double >::const_iterator it = doubledata.find(name);
    if (it !=doubledata.end()) return (*it).second  ;
    else{
        cout << "the double of name " << name << " is not registered in xParseData" << endl;
        throw;
    }
}

void xParseData::setDouble(std::string name, double d){
    std::map<std::string,double >::iterator it = doubledata.find(name);
    if (it !=doubledata.end()) (*it).second=d;
    else{
        cout << "the double of name " << name << " is not registered in xParseData" << endl;
        throw;
    }
}

int xParseData::getInt(std::string name) const{
    std::map<std::string, int >::const_iterator it = intdata.find(name);
    if (it != intdata.end()) return (*it).second  ;
    else{
        std::cout << "the integer of name " << name << " is not registered in xParseData" << std::endl;
        throw;
    }
}

void xParseData::setInt(std::string name, int i) {
    std::map<std::string, int >::iterator it = intdata.find(name);
    if (it != intdata.end()) (*it).second = i;
    else{
        std::cout << "the integer of name " << name << " is not registered in xParseData" << std::endl;
        throw;
    }
}

std::list <std::string  > xParseData::getListString(std::string name) const{
    std::map<std::string, std::list<std::string >  >::const_iterator it = liststringdata.find(name);
    if (it != liststringdata.end()) return (*it).second  ;
    else{
        std::cout << "the list string of name " << name << " is not registered in xParseData" << std::endl;
        throw;
    }
}

std::list <xtensor::xVector<> > xParseData::getListVector(std::string name) const{
    std::map<std::string, std::list<xtensor::xVector<> > >::const_iterator it = listvectordata.find(name);
    if (it != listvectordata.end()) return (*it).second  ;
    else{
        std::cout << "the list vector of name " << name << " is not registered in xParseData" << std::endl;
        throw;
    }
}

std::string  xParseData::getString(std::string name) const{
    std::map<std::string, std::string   >::const_iterator it = stringdata.find(name);
    if (it != stringdata.end()) return (*it).second  ;
    else{
        std::cout << "the string of name " << name << " is not registered in xParseData" << std::endl;
        throw;
    }
}
void xParseData::setString(std::string name,const std::string &s ) 
{
    std::map<std::string,  std::string >::iterator it = stringdata.find(name);
    if (it != stringdata.end()) (*it).second = s;
    else{
        std::cout << "the string of name " << name << " is not registered in xParseData" << std::endl;
        throw;
    }
}

xtensor::xVector<> xParseData::getVector(std::string name) const{
    std::map<std::string, xtensor::xVector<> >::const_iterator it = xvectordata.find(name);
    if (it != xvectordata.end()) return (*it).second  ;
    else{
        std::cout << "the xtensor::xVector of name " << name << " is not registered in xParseData" << std::endl;
        throw;
    }
}

xtensor::xTensor2<> xParseData::getTensor2(std::string name) const{
    std::map<std::string, xtensor::xTensor2<> >::const_iterator it = xtensor2data.find(name);
    if (it != xtensor2data.end()) return (*it).second  ;
    else{
        std::cout << "the xtensor::xTensor2 of name " << name << " is not registered in xParseData" << std::endl;
        throw;
    }
}

std::map<std::string, double> xParseData::getMapStringDouble(std::string name) const{
    std::map<std::string, std::map<std::string, double>  >::const_iterator it = mapstringdoubledata.find(name);
    if (it != mapstringdoubledata.end()) return (*it).second  ;
    else{
        std::cout << "the map<string, double> of name " << name << " is not registered in xParseData" << std::endl;
        throw;
    }
}

xParseData::datatype xParseData::getType(std::string name) const{
    if (intdata.find(name) !=  intdata.end())  return INT;
    if (doubledata.find(name) !=doubledata.end())  return DOUBLE;
    if (stringdata.find(name) !=stringdata.end())  return STRING;
    if (liststringdata.find(name) !=liststringdata.end())  return LISTSTRING;
    if (xvectordata.find(name)!=xvectordata.end()) return XVECTOR;
    if (listvectordata.find(name)!=listvectordata.end()) return LISTXVECTOR;
    if (xtensor2data.find(name)!=xtensor2data.end()) return XTENSOR2;
    if (mapstringdoubledata.find(name)!=mapstringdoubledata.end()) return MAPSTRINGDOUBLE;

    else return UNDEFINED;
}

void nextvalidword(std::istream &is){
    int go=1;
    char bufc[256];
    while (go){
        switch (is.peek()){
        case '/':
            is.ignore();
            if (is.peek()=='/') is.getline(bufc, 256);
            break;
        case '#':
            is.getline(bufc, 256);
            break;
        case ' ':
	case '}':
        case '\n':
        case '\t':
            is.ignore();
            break;
        default :{
            go=0;
            break;
        }
        }
    }
}

int xParseData::parse(std::istream &is){
    std::string buf, name;

    nextvalidword(is);
    while ((is.peek()=='{') ) is.ignore();
    is>> name;
    while (name!="}"&& !is.fail()){
//         std::cout <<  "before sw " << name<< std::endl;
        switch (getType(name)){
        case INT: {
            // std::cout <<  "case int  " << std::endl;
            while ((is.peek()=='=') ||(is.peek()==' ' )) is.ignore();
            int a;
            is >> a;
            if (!is) {
                std::cout << "expecting a int for " <<  name <<  std::endl;
                throw;
            }
            intdata[name] = a;
            isdefined[name] = isdefined[name]+1 ;
            break;
        }
        case DOUBLE: {
            // std::cout <<  "case double  " << std::endl;
            while ((is.peek()=='=') ||(is.peek()==' ' )) is.ignore();
            double a;
            is >> a;
            if (!is) {
                std::cout << "expecting a double for " <<  name <<  std::endl;
                throw;
            }
            doubledata[name] = a;
            isdefined[name] = isdefined[name]+1;
            break;
        }
        case STRING:{
            //std::cout <<  "case string  " << std::endl;
            std::string a;
            while ((is.peek()=='=') ||(is.peek()==' ' )) is.ignore();
            is >> a;
            if (!is) {
                std::cout << "expecting a string for " << name <<" "  << a <<  std::endl;
                throw;
            }
            stringdata[name] = a;
            isdefined[name] = isdefined[name]+1 ;
            break;
        }
        case LISTSTRING:{
            //std::cout <<  "case list string  " << std::endl;
            std::string a;
            std::list<std::string > liststring;
            while ((is.peek()=='=') ||(is.peek()==' ' )||(is.peek()==',' )||(is.peek()==';' ) ) is.ignore();
            if (is.peek() == '{' ) is.ignore();
            else std::cout << "stringlist must start with { and end with }" << std::endl;
            while(1){
                while ((is.peek()=='=') ||(is.peek()==' ' )||(is.peek()==',' )||(is.peek()==';' ) ) is.ignore();
                if (is.peek()=='}') break;
                is >> a;
                if (!is) {
                    std::cout << "liststring must end with } " << name <<  std::endl;
                    throw;
                }
                liststring.push_back(a);
            }
            liststringdata[name] = liststring;
            isdefined[name] = isdefined[name]+1 ;
            break;
        }
        case XVECTOR:{
            // std::cout <<  "case xvector  " << std::endl;
            std::string buf;
            while ((is.peek()=='=') ||(is.peek()==' ' )) is.ignore();
            xtensor::xVector<> a;
            is >> a;
            if (!is) {
                std::cout << "expecting a vector for " << name << std::endl;
                throw;
            }
            xvectordata[name] = a;
            isdefined[name] = isdefined[name]+1 ;
            is.ignore();
            break;
        }
        case LISTXVECTOR:{
	    //std::cout << "case listvector" << std::endl;
	    xtensor::xVector<> vec;
	    std::list<xtensor::xVector<> > listvector;
	    while ((is.peek()=='=') || (is.peek()==' ') || (is.peek()==',') || (is.peek()==';') ) is.ignore();
	    if (is.peek()=='{') is.ignore();
	    else std::cout << "vectorlist must start with { and end with }" << std::endl;
	    while(1){
		while ((is.peek()=='=') || (is.peek()==' ') || (is.peek()==',') || (is.peek()==';') ) is.ignore();
		if (is.peek()=='}') break;
		is >> vec;
		if (!is) {
		    std::cout << "listvector must end with } " << name << std::endl;
		    throw;
		}
		listvector.push_back(vec);
	    }
	    listvectordata[name] = listvector;
	    isdefined[name] = isdefined[name]+1;
	    is.ignore();
	    break;
        }
        case XTENSOR2:{
            // std::cout <<  "case xtensor2  " << std::endl;
            std::string buf;
            while ((is.peek()=='=') ||(is.peek()==' ' )) is.ignore();
            xtensor::xTensor2<> a;
            read_tensor2(is, a);
            if (!is) {
                std::cout << "expecting a tensor2 for " << name << std::endl;
                throw;
            }
            xtensor2data[name] = a;
            isdefined[name] = isdefined[name]+1 ;
            is.ignore();
            break;
        }

        case MAPSTRINGDOUBLE:{
            std::cout <<  "case map<string, double>  " << std::endl;
            std::string k;
            double v;
            std::map<std::string, double> mapstringdouble;

            while ((is.peek()=='=') ||(is.peek()==' ' )||(is.peek()==',' )||(is.peek()==';' ) ) is.ignore();

            if (is.peek() == '{' ) is.ignore();
            else std::cout << "map<string, double> must start with { and end with }" << std::endl;

            while(1){

                while ((is.peek()=='=') ||(is.peek()==' ' )||(is.peek()==',' )||(is.peek()==';' ) ) is.ignore();
                if (is.peek()=='}') break;
                //Key:
                is >> k;
                while ((is.peek()==' ' )||(is.peek()==':')) is.ignore();

                //Value:
                is >> v;

                //Store in the map
                mapstringdouble[k] = v;


                cout<<"map:"<<k<<" -> "<<v<<endl;


                if (!is) {
                    std::cout << "map<string, double> must end with } " << name <<  std::endl;
                    throw;
                }
            }
            mapstringdoubledata[name] = mapstringdouble;
            isdefined[name] = isdefined[name]+1 ;
            break;
        }


        case UNDEFINED:{
            std::cout<<"_" <<name << "_ was not registered"<< std::endl;
            throw buf;
        }
        }

        nextvalidword(is);

        //if (is.fail()) throw;
        is >> name;
        //std::cout<< "name " << name << std::endl;

    }

    std::map<std::string, int>::iterator it = isdefined.begin();
    while (it!=isdefined.end()){
        int n = (*it).second;
        if (n==0) std::cout << "Warning : parameter " << (*it).first  << " was not set" << std::endl;
        if (n>1) std::cout  <<  "Warning : parameter " << (*it).first << " was set more than once" << std::endl;
        ++it;
    }


    return 1;
}
