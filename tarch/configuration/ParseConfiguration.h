// Copyright (C) 2015 Technische Universitaet Muenchen
// This file is part of the Mamico project. For conditions of distribution
// and use, please see the copyright notice in Mamico's main folder, or at
// www5.in.tum.de/mamico
#ifndef _TARCH_CONFIGURATION_PARSECONFIGURATION_H_
#define _TARCH_CONFIGURATION_PARSECONFIGURATION_H_
#include "tarch/la/Vector.h"
#include "tarch/tinyxml2/tinyxml2.h"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>

namespace tarch {
namespace configuration {

/**	interface for configuration using tinyxml2
 *  @author Philipp Neumann
 */
class ParseConfiguration {
public:
  /**
   * Container for the XML configuration for MaMiCo
   * (Manages the memory of the associated XML document and provides a pointer to the configuration root)
   */
  class XMLConfiguration {
  private:
    XMLConfiguration(tinyxml2::XMLNode* _root, tinyxml2::XMLError _error) : root(_root), error(_error) {}

  public:
    /**
     * Root of the loaded MaMiCo configuration
     */
    tinyxml2::XMLNode* const root;

    /**
     * Status of loading the XML document
     */
    tinyxml2::XMLError const error;

    /**
     * Loads the MaMiCo configuration from a file (with or without XML root node)
     * @param filename
     */
    static XMLConfiguration load(const std::string filename) {
      tinyxml2::XMLDocument* doc = new tinyxml2::XMLDocument();
      tinyxml2::XMLError error = doc->LoadFile(filename.c_str());
      tinyxml2::XMLNode* root = doc->FirstChildElement("scenario-configuration");
      if (root == NULL) {
        root = doc;
        std::cout << "No root node <scenario-configuration> found in configuration file. (Using legacy format without XML root node.)" << std::endl;
      }
      return XMLConfiguration(root, error);
    }

    ~XMLConfiguration() { delete root->GetDocument(); }
  };

  /** parses an xml configuration file, and stores the information in the object
   *config
   *	@tparam Configuration
   */
  template <class Configuration> static void parseConfiguration(const std::string filename, const std::string topleveltag, Configuration& config) {
    XMLConfiguration xmlConfig = XMLConfiguration::load(filename);
    tinyxml2::XMLElement* node = xmlConfig.root->FirstChildElement(topleveltag.c_str());
    if (node == NULL) {
      std::cout << "Could not read input file " << filename << " (missing or invalid)" << std::endl;
      exit(EXIT_FAILURE);
    }
    config.parseSubtag(node);
  }

  /** reads a double value at node "node" with tag "tag" and stores the result
   *in "storage". the program is exited if the value cannot be read.
   *	@param storage
   *	@param node
   *	@param tag
   */
  static void readDoubleMandatory(double& storage, tinyxml2::XMLElement* node, std::string tag) {
    double value;
    if (node->QueryDoubleAttribute(tag.c_str(), &value) != tinyxml2::XML_SUCCESS) {
      std::cout << "Error while reading mandatory argument " << tag << " of XML element " << node->Name() << std::endl;
      exit(EXIT_FAILURE);
    } else {
      storage = value;
    }
  }

  /** reads a double value at node "node" with tag "tag" and stores the result
   *in "storage". the program continues and only stores the result in "storage"
   *if "tag" is existent
   *	@param storage
   *	@param node
   *	@param tag
   */
  static void readDoubleOptional(double& storage, tinyxml2::XMLElement* node, std::string tag) {
    double value;
    int result = node->QueryDoubleAttribute(tag.c_str(), &value);
    if (result == tinyxml2::XML_NO_ATTRIBUTE) {
      // nop
    } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
      std::cout << "Error while reading optional argument " << tag << " of XML element " << node->Name() << std::endl;
      exit(EXIT_FAILURE);
    } else {
      storage = value;
    }
  }

  /** reads a int value at node "node" with tag "tag" and stores the result in
   *"storage". the program is exited if the value cannot be read.
   *	@param storage
   *	@param node
   *	@param tag
   */
  static void readIntMandatory(int& storage, tinyxml2::XMLElement* node, std::string tag) {
    int value;
    if (node->QueryIntAttribute(tag.c_str(), &value) != tinyxml2::XML_SUCCESS) {
      std::cout << "Error while reading mandatory argument " << tag << " of XML element " << node->Name() << std::endl;
      exit(EXIT_FAILURE);
    } else {
      storage = value;
    }
  }

  /** reads a int value at node "node" with tag "tag" and stores the result in
   *"storage". the program continues and only stores the result in "storage" if
   *"tag" is existent
   *	@param storage
   *	@param node
   *	@param tag
   */
  static void readIntOptional(int& storage, tinyxml2::XMLElement* node, std::string tag) {
    int value;
    int result = node->QueryIntAttribute(tag.c_str(), &value);
    if (result == tinyxml2::XML_NO_ATTRIBUTE) {
      // nop
    } else if (result == tinyxml2::XML_WRONG_ATTRIBUTE_TYPE) {
      std::cout << "Error while reading optional argument " << tag << " of XML element " << node->Name() << std::endl;
      exit(EXIT_FAILURE);
    } else {
      storage = value;
    }
  }

  /** reads a bool value at node "node" with tag "tag" and stores the result in
   *"storage". the program is exited if the value cannot be read.
   *	@param storage
   *	@param node
   *	@param tag
   */
  static void readBoolMandatory(bool& storage, tinyxml2::XMLElement* node, std::string tag) {
    const char* myTextChar = node->Attribute(tag.c_str());
    if (myTextChar == NULL) {
      std::cout << "Error: mandatory bool " << tag << " could not be found!" << std::endl;
      exit(EXIT_FAILURE);
    }
    std::string myText(myTextChar);
    if (myText == "yes") {
      storage = true;
    } else if (myText == "no") {
      storage = false;
    } else {
      std::cout << "Error while reading bool optional argument: Argument can "
                   "only be yes or no!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  /** reads a bool value at node "node" with tag "tag" and stores the result in
   *"storage". the program continues and only stores the result in "storage" if
   *"tag" is existent
   *	@param storage
   *	@param node
   *	@param tag
   */
  static void readBoolOptional(bool& storage, tinyxml2::XMLElement* node, std::string tag) {
    const char* myTextChar = node->Attribute(tag.c_str());
    if (myTextChar == NULL) {
      return;
    }
    std::string myText(myTextChar);
    if (myText == "yes") {
      storage = true;
    } else if (myText == "no") {
      storage = false;
    } else {
      std::cout << "Error while reading bool optional argument: Argument can "
                   "only be yes or no!"
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  /** reads a string value at node "node" with tag "tag" and stores the result
   *in "storage". the program is exited if the value cannot be read.
   *	@param storage
   *	@param node
   *	@param tag
   */
  static void readStringMandatory(std::string& storage, tinyxml2::XMLElement* node, std::string tag) {
    const char* myText = node->Attribute(tag.c_str());
    if (myText == NULL) {
      std::cout << "Error while reading mandatory argument " << tag << " of XML element " << node->Name() << std::endl;
      exit(EXIT_FAILURE);
    } else {
      storage = std::string(myText);
    }
  }

  /** reads a string value at node "node" with tag "tag" and stores the result
   *in "storage". the program continues and only stores the result in "storage"
   *if "tag" is existent
   *	@param storage
   *	@param node
   *	@param tag
   */
  static void readStringOptional(std::string& storage, tinyxml2::XMLElement* node, std::string tag) {
    const char* myText = node->Attribute(tag.c_str());
    if (myText != NULL) {
      storage = std::string(myText);
    }
  }

  /** reads with ";"-separated entry "tag" with "size" entries and returns the
   *corresponding vector on success
   *	@tparam T Data type
   *	@tparam Size size of the entry
   *	@param result
   *	@param myText
   */
  template <unsigned int size, class T> static void readVector(tarch::la::Vector<size, T>& result, const char* myText) {
    std::string input(myText);
    for (unsigned int i = 0; i < size; i++) {
      // search for first non-whitespace entry in this vector entry
      std::size_t first = input.find_first_not_of(" ");
      // search for end of this vector component, typically denoted by ";"
      // -> if this is the last entry, then npos is accepted as well
      std::size_t last = input.find_first_of(";");
      // for debugging
      // std::cout << first << ", " << last << std::endl;
      if ((i == size - 1) && (last == std::string::npos)) {
        last = input.size();
      }

      std::stringstream ss(input.substr(first, last - first));
      // for debugging
      // std::cout << ss.str() << std::endl;
      ss >> result[i];
      if (i < size - 1) {
        input = input.substr(last + 1, input.size() - last - 1);
      }
    }
  }

  /** reads with ";"-separated entry "tag" with "size" entries and returns the
   *corresponding vector on success
   *	@tparam T Data type
   *	@tparam Size size of the entry
   *	@param result
   *	@param node
   *	@param tag
   */
  template <unsigned int size, class T> static void readVectorMandatory(tarch::la::Vector<size, T>& result, tinyxml2::XMLElement* node, std::string tag) {
    const char* myText = node->Attribute(tag.c_str());
    if (myText == NULL) {
      std::cout << "Error while reading mandatory argument " << tag << " of XML element " << node->Name() << std::endl;
      exit(EXIT_FAILURE);
    }
    readVector<size, T>(result, myText);
  }

  /** reads with ";"-separated entry "tag" with "size" entries and returns the
   *corresponding vector on success
   *if "tag" is existent
   *	@tparam T Data type
   *	@tparam Size size of the entry
   *	@param result
   *	@param node
   *	@param tag
   */
  template <unsigned int size, class T> static void readVectorOptional(tarch::la::Vector<size, T>& result, tinyxml2::XMLElement* node, std::string tag) {
    const char* myText = node->Attribute(tag.c_str());
    if (myText != NULL) {
      readVector<size, T>(result, myText);
    }
  }
};
} // namespace configuration
} // namespace tarch
#endif
