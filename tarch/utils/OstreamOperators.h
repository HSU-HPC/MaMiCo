#pragma once

#include <map>
#include <set>
#include <vector>

template <class T> std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
  os << "[";
  for (auto& item : v) {
    os << " " << item;
  }
  os << "]";
  return os;
}

template <class T> std::ostream& operator<<(std::ostream& os, const std::set<T>& s) {
  os << "{";
  for (auto& item : s) {
    os << item << ", ";
  }
  os << "}";
  return os;
}

template <class K, class V> std::ostream& operator<<(std::ostream& os, const std::map<K, V>& m) {
  os << "{";
  for (auto& item : m) {
    os << item.first << ": " << item.second << ",";
  }
  os << "}";
  return os;
}