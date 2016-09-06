#include "ClassLabel.hpp"

std::string getClassColorString(const ClassLabel c) {
  switch (c) {
    case 0: return "black";
    case 1: return "green";
    case 2: return "red";
    case 3: return "yellow";
    case 4: return "blue";
    default: return "white";
  }
}

void getClassColorRGB(const ClassLabel c, unsigned int &R, unsigned int &G, unsigned int &B) {
  switch (c) {
    case 0: R=0; G=0; B=0; break;
    case 1: R=0; G=255; B=0; break;
    case 2: R=255; G=0; B=0; break;
    case 3: R=255; G=255; B=0; break;
    case 4: R=0; G=0; B=255; break;
    default: R=255; G=255; B=255; break;
  }
}
