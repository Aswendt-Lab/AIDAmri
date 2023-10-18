""" Dictionary to XML - Library to convert a python dictionary to XML output
    Copyleft (C) 2007 Pianfetti Maurizio <boymix81@gmail.com>
    Package site : http://boymix81.altervista.org/files/dict2xml.tar.gz

    Revision 1.0  2007/12/15 11:57:20  Maurizio
    - First stable version

"""

__author__ = "Pianfetti Maurizio <boymix81@gmail.com>"
__contributors__ = []
__date__    = "$Date: 2007/12/15 11:57:20  $"
__credits__ = """..."""
__version__ = "$Revision: 1.0.0 $"

class Dict2XML:
    #XML output
    xml = ""

    #Tab level
    level = 0

    def __init__(self):
        self.xml = ""
        self.level = 0
    #end def

    def __del__(self):
        pass
    #end def

    def setXml(self,Xml):
        self.xml = Xml
    #end if

    def setLevel(self,Level):
        self.level = Level
    #end if

    def dict2xml(self,map): # reserved assignment
        if (str(type(map)) == "<class 'object_dict.object_dict'>" or str(type(map)) == "<type 'dict'>"):
            for key, value in map.items():
                if (str(type(value)) == "<class 'object_dict.object_dict'>" or str(type(value)) == "<type 'dict'>"):
                    if(len(value) > 0):
                        self.xml += "\t"*self.level
                        self.xml += "<%s>\n" % (key)
                        self.level += 1
                        self.dict2xml(value)
                        self.level -= 1
                        self.xml += "\t"*self.level
                        self.xml += "</%s>\n" % (key)
                    else:
                        self.xml += "\t"*(self.level)
                        self.xml += "<%s></%s>\n" % (key,key)
                    #end if
                else:
                    self.xml += "\t"*(self.level)
                    self.xml += "<%s>%s</%s>\n" % (key,value, key)
                #end if
            else:
                self.xml += "\t"*self.level
                self.xml += "<%s>%s</%s>\n" % (key,value, key)
        #end if
        return self.xml
    #end def

#end class

def createXML(dict,xml): # reserved assignment
    xmlout = Dict2XML()
    xmlout.setXml(xml)
    return xmlout.dict2xml(dict)
#end def

dict2Xml = createXML

if __name__ == "__main__":

    #Define the dict
    d={}
    d['root'] = {}
    d['root']['v1'] = "";
    d['root']['v2'] = "hi";
    d['root']['v3'] = {};
    d['root']['v3']['v31']="hi";

    #xml='<?xml version="1.0"?>\n'
    xml = ""
    print(dict2Xml(d,xml))

#end if
