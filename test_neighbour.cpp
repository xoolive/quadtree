#include <cstdio>
#include <cstdlib>

#include "logger.h"
#include "neighbour.h"


int main() {

  Logger log(__FILE__);

  log.testint(__LINE__, Neighbour::samelevel(1, NORTH, 1), 3,
              "samelevel(1, NORTH, 1)");
  log.testint(__LINE__, Neighbour::samelevel(1, NORTHWEST, 1), 2,
              "samelevel(1, NORTHWEST, 1)");
  log.testint(__LINE__, Neighbour::samelevel(1, WEST, 1), 0,
              "samelevel(1, WEST, 1)");
  log.testint(__LINE__, Neighbour::samelevel(1, EAST, 1), 4,
              "samelevel(1, EAST, 1)");

  log.testhex(__LINE__, Neighbour::samelevel(0x3a, EAST, 3)      , 0x3b,
              "samelevel(0x3a, EAST, 3)     ");
  log.testhex(__LINE__, Neighbour::samelevel(0x3a, WEST, 3)      , 0x2f,
              "samelevel(0x3a, WEST, 3)     ");
  log.testhex(__LINE__, Neighbour::samelevel(0x3a, SOUTHWEST, 3) , 0x2d,
              "samelevel(0x3a, SOUTHWEST, 3)");
  log.testhex(__LINE__, Neighbour::samelevel(0x3a, SOUTH, 3)     , 0x38,
              "samelevel(0x3a, SOUTH, 3)    ");
  log.testhex(__LINE__, Neighbour::samelevel(0x3a, SOUTHEAST, 3) , 0x39,
              "samelevel(0x3a, SOUTHEAST, 3)");

  log.testhex(__LINE__, Neighbour::samelevel(0x66, WEST, 4),      0x63,
              "samelevel(0x66, WEST, 4),    ");
  log.testhex(__LINE__, Neighbour::samelevel(0x66, SOUTHWEST, 4), 0x61,
              "samelevel(0x66, SOUTHWEST, 4)");
  log.testhex(__LINE__, Neighbour::samelevel(0x66, SOUTH, 4),     0x64,
              "samelevel(0x66, SOUTH, 4),   ");
  log.testhex(__LINE__, Neighbour::samelevel(0x66, SOUTHEAST, 4), 0x65,
              "samelevel(0x66, SOUTHEAST, 4)");
  log.testhex(__LINE__, Neighbour::samelevel(0x66, EAST, 4),      0x67,
              "samelevel(0x66, EAST, 4),    ");
  log.testhex(__LINE__, Neighbour::samelevel(0x66, NORTHEAST, 4), 0x6d,
              "samelevel(0x66, NORTHEAST, 4)");
  log.testhex(__LINE__, Neighbour::samelevel(0x66, NORTH, 4),     0x6c,
              "samelevel(0x66, NORTH, 4),   ");
  log.testhex(__LINE__, Neighbour::samelevel(0x66, NORTHWEST, 4), 0x69,
              "samelevel(0x66, NORTHWEST, 4)");

  log.testhex(__LINE__, Neighbour::samelevel(0x3a, EAST, 3)      , 0x3b,
              "samelevel(0x3a, EAST, 3)      ");
  log.testhex(__LINE__, Neighbour::samelevel(0x3a, WEST, 3)      , 0x2f,
              "samelevel(0x3a, WEST, 3)      ");
  log.testhex(__LINE__, Neighbour::samelevel(0x3a, SOUTHWEST, 3) , 0x2d,
              "samelevel(0x3a, SOUTHWEST, 3) ");
  log.testhex(__LINE__, Neighbour::samelevel(0x3a, SOUTH, 3)     , 0x38,
              "samelevel(0x3a, SOUTH, 3)     ");
  log.testhex(__LINE__, Neighbour::samelevel(0x3a, SOUTHEAST, 3) , 0x39,
              "samelevel(0x3a, SOUTHEAST, 3) ");

  log.testint(__LINE__, Neighbour::samelevel(1, NORTH, 1), 3,
              "samelevel(1, NORTH, 1)");
  log.testint(__LINE__, Neighbour::samelevel(1, NORTHWEST, 1), 2,
              "samelevel(1, NORTHWEST, 1)");
  log.testint(__LINE__, Neighbour::samelevel(1, WEST, 1), 0,
              "samelevel(1, WEST, 1)");
  log.testint(__LINE__, Neighbour::samelevel(1, EAST, 1), 4,
              "samelevel(1, EAST, 1)");

  return log.reportexit();

}
