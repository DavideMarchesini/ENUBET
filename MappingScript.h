
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TString.h"

std::map<std::vector<int>, std::vector<int>> MappingAnodeCord(TString fileName) {

  std::ifstream infile(fileName);

  std::vector<int> board_anode;
  std::vector<int> coordinates;

  std::map<std::vector<int>, std::vector<int>> myMap;

  std::string board;
  std::string anode;
  std::string r;
  std::string phi;
  std::string z;

  while (infile >> board >> anode >> r >> phi >> z) {

    int i;

    if (r == "/" || phi == "/" || z == "/")
      continue;
    else if (r == "0d")
      i = -1;
    else if (r == "0u")
      i = 0;
    else
      std::istringstream(r) >> i;
	
    coordinates.push_back(i + 1);

    std::istringstream(board) >> i;
    board_anode.push_back(i);

    std::istringstream(anode) >> i;
    board_anode.push_back(i);

    std::istringstream(phi) >> i;
    coordinates.push_back(i);

    std::istringstream(z) >> i;
    coordinates.push_back(i);

    myMap[board_anode] = coordinates;

    board_anode.clear();
    coordinates.clear();
  }

  return myMap;
}
std::map<std::vector<int>, std::vector<int>> MappingCordAnode(TString fileName) {

  std::ifstream infile(fileName);

  std::vector<int> board_anode;
  std::vector<int> coordinates;

  std::map<std::vector<int>, std::vector<int>> myMap;

  std::string board;
  std::string anode;
  std::string r;
  std::string phi;
  std::string z;

  while (infile >> board >> anode >> r >> phi >> z) {

    int i;
    if (r == "/" || phi == "/" || z == "/")
      continue;
    else if (r == "0d")
      i = -1;
    else if (r == "0u")
      i = 0;
    else
      std::istringstream(r) >> i;
    coordinates.push_back(i + 1);

    std::istringstream(board) >> i;
    board_anode.push_back(i);

    std::istringstream(anode) >> i;
    board_anode.push_back(i);

    std::istringstream(phi) >> i;
    coordinates.push_back(i);

    std::istringstream(z) >> i;
    coordinates.push_back(i);

    myMap[coordinates] = board_anode;

    board_anode.clear();
    coordinates.clear();
  }

  return myMap;
}
std::map<int, std::vector<int>> MappingDigiCord(TString fileName) {

  std::ifstream infile(fileName);

  std::vector<int> coordinates;

  std::map<int, std::vector<int>> myMap;

  std::string digiboard;
  std::string ch;
  std::string r;
  std::string phi;
  std::string z;

  while (infile >> digiboard >> ch >> r >> phi >> z) {

    int i;
    int board_ch;

    if (r == "/" || phi == "/" || z == "/")
      continue;
    else if (r == "0d")
      i = -1;
    else if (r == "0u")
      i = 0;
    else
      std::istringstream(r) >> i;
    coordinates.push_back(i + 1);

    std::istringstream(ch) >> board_ch;

    std::istringstream(phi) >> i;
    coordinates.push_back(i);

    std::istringstream(z) >> i;
    coordinates.push_back(i);

    myMap[board_ch] = coordinates;

    coordinates.clear();
  }

  return myMap;
}

void printMap(std::map<std::vector<int>, std::vector<int>> map) {
  for (auto elem = map.begin(); elem != map.end(); elem++) {
    for (int i = 0; i < elem->first.size(); i++)
      std::cout << elem->first.at(i) << ", ";
    std::cout << " :::: ";
    for (int i = 0; i < elem->second.size(); i++)
      std::cout << elem->second.at(i) << ", ";
    std::cout << std::endl;
  }
  return;
}
