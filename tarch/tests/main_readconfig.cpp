// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico

#include "tarch/configuration/ParseConfiguration.h"

int main() {
  std::string filename = "config.xml";
  tinyxml2::XMLDocument file;
  file.LoadFile(filename.c_str());
  tinyxml2::XMLElement *node = file.FirstChildElement("mytag");
  double mydouble;
  int myint;
  tarch::la::Vector<3, int> myIntVec;
  tarch::la::Vector<2, double> myDblVec;
  bool mybool = false;

  tarch::configuration::ParseConfiguration::readDoubleMandatory(mydouble, node,
                                                                "mydouble");
  std::cout << "mydouble=" << mydouble << std::endl;
  tarch::configuration::ParseConfiguration::readIntMandatory(myint, node,
                                                             "myint");
  std::cout << "myint=" << myint << std::endl;
  tarch::configuration::ParseConfiguration::readVector<2, double>(
      myDblVec, node, "myvecdbl");
  std::cout << "myvecdbl=" << myDblVec << std::endl;
  tarch::configuration::ParseConfiguration::readVector<3, int>(myIntVec, node,
                                                               "myvecint");
  std::cout << "myvecint=" << myIntVec << std::endl;
  tarch::configuration::ParseConfiguration::readBoolMandatory(mybool, node,
                                                              "mybool");
  std::cout << "mybool=" << mybool << std::endl;
}
