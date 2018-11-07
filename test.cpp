/**********************************************************************
  Test script for XtalComp

  Copyright (C) 2011 by David C. Lonie

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
***********************************************************************/

#include "xtalcomp.h"

#include <stdio.h>
#include <stdlib.h>

void runTest(bool (*testFunc)(), const char * testName, int &successes,
             int &failures)
{
  if ((*testFunc)()) {
    fprintf(stdout, "Test '%s' passes.\n", testName);
    ++successes;
  }
  else {
    fprintf(stdout, "Test '%s' fails!\n", testName);
    ++failures;
  }
}

bool simpleCase()
{
  XcMatrix cell1 ( 3.0, 0.0, 0.0, 2.0, 4.0, 0.0, 2.0, 5.0, 3.0 );
  XcMatrix cell2 (cell1);

  std::vector<XcVector> pos1;
  pos1.reserve(4);

  pos1.push_back(XcVector(0.0, 0.25, 0.25));
  pos1.push_back(XcVector(0.25, 0.25, 0.25));
  pos1.push_back(XcVector(0.0, 0.5, 0.25));
  pos1.push_back(XcVector(0.0, 0.25, 0.75));
  std::vector<XcVector> pos2 (pos1);

  std::vector<unsigned int> types1;
  types1.reserve(4);

  types1.push_back(1);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(1);
  std::vector<unsigned int> types2 (types1);

  bool match = XtalComp::compare(cell1, types1, pos1,
                                 cell2, types2, pos2,
                                 NULL, 0.05, 0.25);

  if (!match)
    return false;

  // Displace an atom, ensure that comparison fails.
  pos2[0] += XcVector(0.5,0,0);
  match = XtalComp::compare(cell1, types1, pos1,
                            cell2, types2, pos2,
                            NULL, 0.05, 0.25);
  if (match)
    return false;

  return true;
}

bool simpleNiggli()
{
  XcMatrix cell1 ( 3.0, 0.0, 0.0, 2.0, 4.0, 0.0, 2.0, 5.0, 3.0 );
  XcMatrix cell2 (cell1);

  std::vector<XcVector> pos1;
  pos1.reserve(4);

  pos1.push_back(XcVector(0.0, 0.25, 0.25));
  pos1.push_back(XcVector(0.25, 0.25, 0.25));
  pos1.push_back(XcVector(0.0, 0.5, 0.25));
  pos1.push_back(XcVector(0.0, 0.25, 0.75));
  std::vector<XcVector> pos2 (pos1);

  std::vector<unsigned int> types1;
  types1.reserve(4);

  types1.push_back(1);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(1);
  std::vector<unsigned int> types2 (types1);

  // Modify the structure 2
  const XcMatrix linComb (1.0, 1.0, 0.0,
                          1.0, 0.0, 1.0,
                          0.0, 0.0, 1.0);
  const XcMatrix xform (0.0, -1.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 0.0, -1.0);
  XcMatrix fcoordUpdate (xform * cell2.transpose());
  cell2 = linComb * cell2 * xform.transpose();
  fcoordUpdate = cell2.transpose().inverse() * fcoordUpdate;
  for (std::vector<XcVector>::iterator it = pos2.begin(), it_end = pos2.end();
       it != it_end; ++it) {
    *it = fcoordUpdate * (*it);
  }

  bool match = XtalComp::compare(cell1, types1, pos1,
                                 cell2, types2, pos2,
                                 NULL, 0.05, 0.25);
  if (!match)
    return false;

  // Displace an atom, ensure that comparison fails.
  pos2[0] += XcVector(0.5,0,0);
  match = XtalComp::compare(cell1, types1, pos1,
                            cell2, types2, pos2,
                            NULL, 0.05, 0.25);
  if (match)
    return false;

  return true;
}

bool simpleUniformTranslation()
{
  XcMatrix cell1 ( 3.0, 0.0, 0.0, 2.0, 4.0, 0.0, 2.0, 5.0, 3.0 );
  XcMatrix cell2 (cell1);

  std::vector<XcVector> pos1;
  pos1.reserve(4);

  pos1.push_back(XcVector(0.0, 0.25, 0.25));
  pos1.push_back(XcVector(0.25, 0.25, 0.25));
  pos1.push_back(XcVector(0.0, 0.5, 0.25));
  pos1.push_back(XcVector(0.0, 0.25, 0.75));
  std::vector<XcVector> pos2 (pos1);

  std::vector<unsigned int> types1;
  types1.reserve(4);

  types1.push_back(1);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(1);
  std::vector<unsigned int> types2 (types1);

  // Displace all coordinates in types2 by a random vector
  srand(0); // Intentionally seeding with zero
  XcVector disp (rand() / static_cast<double>(RAND_MAX),
                 rand() / static_cast<double>(RAND_MAX),
                 rand() / static_cast<double>(RAND_MAX));
  for (std::vector<XcVector>::iterator it = pos2.begin(), it_end = pos2.end();
       it != it_end; ++it) {
    (*it) += disp;
  }

  bool match = XtalComp::compare(cell1, types1, pos1,
                                 cell2, types2, pos2,
                                 NULL, 0.05, 0.25);
  if (!match)
    return false;

  // Displace an atom, ensure that comparison fails.
  pos2[0] += XcVector(0.5,0,0);
  match = XtalComp::compare(cell1, types1, pos1,
                            cell2, types2, pos2,
                            NULL, 0.05, 0.25);
  if (match)
    return false;

  return true;
}

bool simpleRandomNoise()
{
  XcMatrix cell1 ( 3.0, 0.0, 0.0, 2.0, 4.0, 0.0, 2.0, 5.0, 3.0 );
  XcMatrix cell2 (cell1);

  std::vector<XcVector> pos1;
  pos1.reserve(4);

  pos1.push_back(XcVector(0.0, 0.25, 0.25));
  pos1.push_back(XcVector(0.25, 0.25, 0.25));
  pos1.push_back(XcVector(0.0, 0.5, 0.25));
  pos1.push_back(XcVector(0.0, 0.25, 0.75));
  std::vector<XcVector> pos2 (pos1);

  std::vector<unsigned int> types1;
  types1.reserve(4);

  types1.push_back(1);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(1);
  std::vector<unsigned int> types2 (types1);

  // Displace all coordinates in types2 by a random vector
  const double cartesianNoiseMax = 0.005;
  srand(0); // Intentionally seeding with zero
  for (std::vector<XcVector>::iterator it = pos2.begin(), it_end = pos2.end();
       it != it_end; ++it) {
    XcVector disp (rand() / static_cast<double>(RAND_MAX),
                   rand() / static_cast<double>(RAND_MAX),
                   rand() / static_cast<double>(RAND_MAX));
    // Convert displacement to cartesian units
    disp = cell2.transpose() * disp;
    // Normalize and set length
    disp *= (rand() / static_cast<double>(RAND_MAX))
        * cartesianNoiseMax / disp.norm();
    // Convert back to fractional units
    disp = cell2.transpose().inverse() * disp;
    // Displace vector
    (*it) += disp;
  }

  bool match = XtalComp::compare(cell1, types1, pos1,
                                 cell2, types2, pos2,
                                 NULL, 0.05, 0.25);
  if (!match)
    return false;

  // Displace an atom, ensure that comparison fails.
  pos2[0] += XcVector(0.5,0,0);
  match = XtalComp::compare(cell1, types1, pos1,
                            cell2, types2, pos2,
                            NULL, 0.05, 0.25);
  if (match)
    return false;

  return true;
}

bool allOfTheAbove()
{
  XcMatrix cell1 ( 3.0, 0.0, 0.0, 2.0, 4.0, 0.0, 2.0, 5.0, 3.0 );
  XcMatrix cell2 (cell1);

  std::vector<XcVector> pos1;
  pos1.reserve(4);

  pos1.push_back(XcVector(0.0, 0.25, 0.25));
  pos1.push_back(XcVector(0.25, 0.25, 0.25));
  pos1.push_back(XcVector(0.0, 0.5, 0.25));
  pos1.push_back(XcVector(0.0, 0.25, 0.75));
  std::vector<XcVector> pos2 (pos1);

  std::vector<unsigned int> types1;
  types1.reserve(4);

  types1.push_back(1);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(1);
  std::vector<unsigned int> types2 (types1);

  //******************************************************************
  // From niggli test:
  // Modify the structure 2
  const XcMatrix linComb (1.0, 1.0, 0.0,
                          1.0, 0.0, 1.0,
                          0.0, 0.0, 1.0);
  const XcMatrix xform (0.0, -1.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 0.0, -1.0);
  XcMatrix fcoordUpdate (xform * cell2.transpose());
  cell2 = linComb * cell2 * xform.transpose();
  fcoordUpdate = cell2.transpose().inverse() * fcoordUpdate;
  for (std::vector<XcVector>::iterator it = pos2.begin(), it_end = pos2.end();
       it != it_end; ++it) {
    *it = fcoordUpdate * (*it);
  }
  //******************************************************************

  //******************************************************************
  // From uniform displacement test:
  // Displace all coordinates in types2 by a random vector
  srand(0); // Intentionally seeding with zero
  XcVector disp (rand() / static_cast<double>(RAND_MAX),
                 rand() / static_cast<double>(RAND_MAX),
                 rand() / static_cast<double>(RAND_MAX));
  for (std::vector<XcVector>::iterator it = pos2.begin(), it_end = pos2.end();
       it != it_end; ++it) {
    (*it) += disp;
  }
  //******************************************************************

  //******************************************************************
  // From noise test:
  // Displace all coordinates in types2 by a random vector
  const double cartesianNoiseMax = 0.005;
  srand(0); // Intentionally seeding with zero
  for (std::vector<XcVector>::iterator it = pos2.begin(), it_end = pos2.end();
       it != it_end; ++it) {
    XcVector disp (rand() / static_cast<double>(RAND_MAX),
                   rand() / static_cast<double>(RAND_MAX),
                   rand() / static_cast<double>(RAND_MAX));
    // Convert displacement to cartesian units
    disp = cell2.transpose() * disp;
    // Normalize and set length
    disp *= (rand() / static_cast<double>(RAND_MAX))
        * cartesianNoiseMax / disp.norm();
    // Convert back to fractional units
    disp = cell2.transpose().inverse() * disp;
    // Displace vector
    (*it) += disp;
  }
  //******************************************************************

  bool match = XtalComp::compare(cell1, types1, pos1,
                                  cell2, types2, pos2,
                                  NULL, 0.05, 0.25);
  if (!match)
    return false;

  // Displace an atom, ensure that comparison fails.
  pos2[0] += XcVector(0.5,0,0);
  match = XtalComp::compare(cell1, types1, pos1,
                            cell2, types2, pos2,
                            NULL, 0.05, 0.25);
  if (match)
    return false;

  return true;
}

bool hexagonalCellTest()
{
  XcMatrix cell1 (3.8398, 0.0, 0.0, -1.9199, 3.32536, 0.0, 0.0, 0.0, 5.93459);

  std::vector<XcVector> pos1;
  pos1.reserve(12);
  pos1.push_back(XcVector( 0.33333, 0.66667, 0.56072));
  pos1.push_back(XcVector( 0.66667, 0.33333, 0.43928));
  pos1.push_back(XcVector( 0.66667, 0.33333, 0.06072));
  pos1.push_back(XcVector( 0.33333, 0.66667, 0.93928));
  pos1.push_back(XcVector( 0.16448, 0.83552, 0.25000));
  pos1.push_back(XcVector( 0.83552, 0.16448, 0.75000));
  pos1.push_back(XcVector( 0.00000, 0.00000, 0.00000));
  pos1.push_back(XcVector( 0.00000, 0.00000, 0.50000));
  pos1.push_back(XcVector( 0.16448, 0.32896, 0.25000));
  pos1.push_back(XcVector( 0.83552, 0.67104, 0.75000));
  pos1.push_back(XcVector( 0.67104, 0.83552, 0.25000));
  pos1.push_back(XcVector( 0.32896, 0.16448, 0.75000));

  std::vector<unsigned int> types1;
  types1.reserve(12);
  types1.push_back(1);
  types1.push_back(1);
  types1.push_back(1);
  types1.push_back(1);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(3);
  types1.push_back(3);
  types1.push_back(3);
  types1.push_back(3);

  XcMatrix cell2 (cell1);
  std::vector<XcVector> pos2 (pos1);
  std::swap(pos2[4], pos2[8]);
  std::swap(pos2[5], pos2[9]);
  std::vector<unsigned int> types2 (types1);

  bool match = XtalComp::compare(cell1, types1, pos1,
                                 cell2, types2, pos2,
                                 NULL, 0.05, 0.25);

  if (!match)
    return false;

  // Displace an atom, ensure that comparison fails.
  pos2[0] += XcVector(0.5,0,0);
  match = XtalComp::compare(cell1, types1, pos1,
                            cell2, types2, pos2,
                            NULL, 0.05, 0.25);
  if (match)
    return false;

  return true;
}

// Test for a bug in which XtalComp rotated atoms so that, after wrapping,
// some of them overlapped. Then, both atoms were matched to the same atom
// in the reference xtal which resulted in a false positive.
bool atomsOverlapping()
{
  XcMatrix cell1 ( 5.79828, 0.0, 0.0, 0.0, 5.79828, 0.0, 0.0, 0.0, 8.2 );
  XcMatrix cell2 (cell1);

  std::vector<XcVector> pos1;
  pos1.reserve(20);

  pos1.push_back(XcVector(0.0, 0.0, 0.0));
  pos1.push_back(XcVector(0.0, 0.0, 0.5));
  pos1.push_back(XcVector(0.5, 0.5, 0.0));
  pos1.push_back(XcVector(0.5, 0.5, 0.5));
  pos1.push_back(XcVector(0.5, 0.0, 0.25));
  pos1.push_back(XcVector(0.5, 0.0, 0.75));
  pos1.push_back(XcVector(0.0, 0.5, 0.25));
  pos1.push_back(XcVector(0.0, 0.5, 0.75));
  pos1.push_back(XcVector(0.25, 0.25, 0.25));
  pos1.push_back(XcVector(0.25, 0.25, 0.75));
  pos1.push_back(XcVector(0.25, 0.75, 0.25));
  pos1.push_back(XcVector(0.25, 0.75, 0.75));
  pos1.push_back(XcVector(0.75, 0.25, 0.25));
  pos1.push_back(XcVector(0.75, 0.25, 0.75));
  pos1.push_back(XcVector(0.75, 0.75, 0.25));
  pos1.push_back(XcVector(0.75, 0.75, 0.75));
  pos1.push_back(XcVector(0.5, 0.0, 0.0));
  pos1.push_back(XcVector(0.5, 0.0, 0.5));
  pos1.push_back(XcVector(0.0, 0.5, 0.0));
  pos1.push_back(XcVector(0.0, 0.5, 0.5));

  std::vector<XcVector> pos2;
  pos2.reserve(20);

  pos2.push_back(XcVector(0.0, 0.0, 0.0));
  pos2.push_back(XcVector(0.0, 0.0, 0.5));
  pos2.push_back(XcVector(0.5, 0.5, 0.0));
  pos2.push_back(XcVector(0.5, 0.5, 0.5));
  pos2.push_back(XcVector(0.5, 0.0, 0.25));
  pos2.push_back(XcVector(0.5, 0.0, 0.75));
  pos2.push_back(XcVector(0.0, 0.5, 0.25));
  pos2.push_back(XcVector(0.0, 0.5, 0.75));
  pos2.push_back(XcVector(0.25, 0.25, 0.25));
  pos2.push_back(XcVector(0.25, 0.25, 0.75));
  pos2.push_back(XcVector(0.25, 0.75, 0.25));
  pos2.push_back(XcVector(0.75, 0.25, 0.25));
  pos2.push_back(XcVector(0.75, 0.75, 0.25));
  pos2.push_back(XcVector(0.75, 0.75, 0.75));
  pos2.push_back(XcVector(0.5, 0.0, 0.0));
  pos2.push_back(XcVector(0.0, 0.5, 0.0));
  pos2.push_back(XcVector(0.25, 0.75, 0.75));
  pos2.push_back(XcVector(0.75, 0.25, 0.75));
  pos2.push_back(XcVector(0.5, 0.0, 0.5));
  pos2.push_back(XcVector(0.0, 0.5, 0.5));

  std::vector<unsigned int> types1;
  types1.reserve(20);

  types1.push_back(1);
  types1.push_back(1);
  types1.push_back(1);
  types1.push_back(1);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(2);
  types1.push_back(3);
  types1.push_back(3);
  types1.push_back(3);
  types1.push_back(3);
  types1.push_back(3);
  types1.push_back(3);
  types1.push_back(3);
  types1.push_back(3);
  types1.push_back(4);
  types1.push_back(4);
  types1.push_back(4);
  types1.push_back(4);

  std::vector<unsigned int> types2 (types1);

  bool match = XtalComp::compare(cell2, types2, pos2,
                                 cell1, types1, pos1,
                                 NULL, 0.05, 0.25);

  if (match)
    return false;

  return true;
}

// This tests some cases where Craig Fisher reported that reversing the order
// of the POSCARs into XtalComp gives different results.
bool craigFisherTest()
{
  // The formula for all 4 of these is Al1Si35
  std::vector<unsigned int> types(36, 14);
  types[0] = 13;

  std::vector<XcVector> pos1, pos2, pos3, pos4;
  pos1.reserve(36);
  pos2.reserve(36);
  pos3.reserve(36);
  pos4.reserve(36);

  // structure05.vasp
  XcMatrix cell1(18.1260000000000012,  0.0000000000000000, 0.0000000000000000,
                 -9.0629999999999971, 15.6975764689967363, 0.0000000000000000,
                  0.0000000000000000,  0.0000000000000000, 7.5670000000000002);

  pos1.push_back(XcVector(0.2635000000000000, 0.9041000000000000, 0.5000000));
  pos1.push_back(XcVector(0.1648000000000000, 0.4982000000000000, 0.2030000));
  pos1.push_back(XcVector(0.0959000000000000, 0.3594000000000000, 0.5000000));
  pos1.push_back(XcVector(0.5018000000000000, 0.6666666666666666, 0.2030000));
  pos1.push_back(XcVector(0.6405999999999999, 0.7365000000000000, 0.5000000));
  pos1.push_back(XcVector(0.3333333333333333, 0.8352000000000001, 0.2030000));
  pos1.push_back(XcVector(0.8352000000000001, 0.5018000000000000, 0.2030000));
  pos1.push_back(XcVector(0.9041000000000000, 0.6405999999999999, 0.5000000));
  pos1.push_back(XcVector(0.4982000000000000, 0.3333333333333333, 0.2030000));
  pos1.push_back(XcVector(0.3594000000000000, 0.2635000000000000, 0.5000000));
  pos1.push_back(XcVector(0.6666666666666666, 0.1648000000000000, 0.2030000));
  pos1.push_back(XcVector(0.7365000000000000, 0.0959000000000000, 0.5000000));
  pos1.push_back(XcVector(0.5018000000000000, 0.8352000000000001, 0.2030000));
  pos1.push_back(XcVector(0.6405999999999999, 0.9041000000000000, 0.5000000));
  pos1.push_back(XcVector(0.3333333333333333, 0.4982000000000000, 0.2030000));
  pos1.push_back(XcVector(0.2635000000000000, 0.3594000000000000, 0.5000000));
  pos1.push_back(XcVector(0.1648000000000000, 0.6666666666666666, 0.2030000));
  pos1.push_back(XcVector(0.0959000000000000, 0.7365000000000000, 0.5000000));
  pos1.push_back(XcVector(0.4982000000000000, 0.1648000000000000, 0.2030000));
  pos1.push_back(XcVector(0.3594000000000000, 0.0959000000000000, 0.5000000));
  pos1.push_back(XcVector(0.6666666666666666, 0.5018000000000000, 0.2030000));
  pos1.push_back(XcVector(0.7365000000000000, 0.6405999999999999, 0.5000000));
  pos1.push_back(XcVector(0.8352000000000001, 0.3333333333333333, 0.2030000));
  pos1.push_back(XcVector(0.9041000000000000, 0.2635000000000000, 0.5000000));
  pos1.push_back(XcVector(0.8352000000000001, 0.5018000000000000, 0.7970000));
  pos1.push_back(XcVector(0.4982000000000000, 0.3333333333333333, 0.7970000));
  pos1.push_back(XcVector(0.6666666666666666, 0.1648000000000000, 0.7970000));
  pos1.push_back(XcVector(0.1648000000000000, 0.4982000000000000, 0.7970000));
  pos1.push_back(XcVector(0.5018000000000000, 0.6666666666666666, 0.7970000));
  pos1.push_back(XcVector(0.3333333333333333, 0.8352000000000001, 0.7970000));
  pos1.push_back(XcVector(0.4982000000000000, 0.1648000000000000, 0.7970000));
  pos1.push_back(XcVector(0.6666666666666666, 0.5018000000000000, 0.7970000));
  pos1.push_back(XcVector(0.8352000000000001, 0.3333333333333333, 0.7970000));
  pos1.push_back(XcVector(0.5018000000000000, 0.8352000000000001, 0.7970000));
  pos1.push_back(XcVector(0.3333333333333333, 0.4982000000000000, 0.7970000));
  pos1.push_back(XcVector(0.1648000000000000, 0.6666666666666666, 0.7970000));

  // structure05-niggli.vasp
  XcMatrix cell2(0.00000000000000,   0.00000000000000, 7.56700000000000,
                 9.06300000000000, -15.69757600000000, 0.00000000000000,
                 9.06300000000000,  15.69757600000000, 0.00000000000000);

  pos2.push_back(XcVector(0.50000, -0.64060, 0.26350));
  pos2.push_back(XcVector(0.20300, -0.33340, 0.16480));
  pos2.push_back(XcVector(0.50000, -0.26350, 0.09590));
  pos2.push_back(XcVector(0.20300, -0.16487, 0.50180));
  pos2.push_back(XcVector(0.50000, -0.09590, 0.64060));
  pos2.push_back(XcVector(0.20300, -0.50187, 0.33333));
  pos2.push_back(XcVector(0.20300,  0.33340, 0.83520));
  pos2.push_back(XcVector(0.50000,  0.26350, 0.90410));
  pos2.push_back(XcVector(0.20300,  0.16487, 0.49820));
  pos2.push_back(XcVector(0.50000,  0.09590, 0.35940));
  pos2.push_back(XcVector(0.20300,  0.50187, 0.66667));
  pos2.push_back(XcVector(0.50000,  0.64060, 0.73650));
  pos2.push_back(XcVector(0.20300, -0.33340, 0.50180));
  pos2.push_back(XcVector(0.50000, -0.26350, 0.64060));
  pos2.push_back(XcVector(0.20300, -0.16487, 0.33333));
  pos2.push_back(XcVector(0.50000, -0.09590, 0.26350));
  pos2.push_back(XcVector(0.20300, -0.50187, 0.16480));
  pos2.push_back(XcVector(0.50000, -0.64060, 0.09590));
  pos2.push_back(XcVector(0.20300,  0.33340, 0.49820));
  pos2.push_back(XcVector(0.50000,  0.26350, 0.35940));
  pos2.push_back(XcVector(0.20300,  0.16487, 0.66667));
  pos2.push_back(XcVector(0.50000,  0.09590, 0.73650));
  pos2.push_back(XcVector(0.20300,  0.50187, 0.83520));
  pos2.push_back(XcVector(0.50000,  0.64060, 0.90410));
  pos2.push_back(XcVector(0.79700,  0.33340, 0.83520));
  pos2.push_back(XcVector(0.79700,  0.16487, 0.49820));
  pos2.push_back(XcVector(0.79700,  0.50187, 0.66667));
  pos2.push_back(XcVector(0.79700, -0.33340, 0.16480));
  pos2.push_back(XcVector(0.79700, -0.16487, 0.50180));
  pos2.push_back(XcVector(0.79700, -0.50187, 0.33333));
  pos2.push_back(XcVector(0.79700,  0.33340, 0.49820));
  pos2.push_back(XcVector(0.79700,  0.16487, 0.66667));
  pos2.push_back(XcVector(0.79700,  0.50187, 0.83520));
  pos2.push_back(XcVector(0.79700, -0.33340, 0.50180));
  pos2.push_back(XcVector(0.79700, -0.16487, 0.33333));
  pos2.push_back(XcVector(0.79700, -0.50187, 0.16480));

  // structure07.vasp
  XcMatrix cell3(18.1260000000000012,  0.0000000000000000, 0.0000000000000000,
                 -9.0629999999999971, 15.6975764689967363, 0.0000000000000000,
                  0.0000000000000000,  0.0000000000000000, 7.5670000000000002);

  pos3.push_back(XcVector(0.9041000000000000, 0.6405999999999999, 0.500000000));
  pos3.push_back(XcVector(0.2635000000000000, 0.9041000000000000, 0.500000000));
  pos3.push_back(XcVector(0.1648000000000000, 0.4982000000000000, 0.203000000));
  pos3.push_back(XcVector(0.0959000000000000, 0.3594000000000000, 0.500000000));
  pos3.push_back(XcVector(0.5018000000000000, 0.6666666666666666, 0.203000000));
  pos3.push_back(XcVector(0.6405999999999999, 0.7365000000000000, 0.500000000));
  pos3.push_back(XcVector(0.3333333333333333, 0.8352000000000001, 0.203000000));
  pos3.push_back(XcVector(0.8352000000000001, 0.5018000000000000, 0.203000000));
  pos3.push_back(XcVector(0.4982000000000000, 0.3333333333333333, 0.203000000));
  pos3.push_back(XcVector(0.3594000000000000, 0.2635000000000000, 0.500000000));
  pos3.push_back(XcVector(0.6666666666666666, 0.1648000000000000, 0.203000000));
  pos3.push_back(XcVector(0.7365000000000000, 0.0959000000000000, 0.500000000));
  pos3.push_back(XcVector(0.5018000000000000, 0.8352000000000001, 0.203000000));
  pos3.push_back(XcVector(0.6405999999999999, 0.9041000000000000, 0.500000000));
  pos3.push_back(XcVector(0.3333333333333333, 0.4982000000000000, 0.203000000));
  pos3.push_back(XcVector(0.2635000000000000, 0.3594000000000000, 0.500000000));
  pos3.push_back(XcVector(0.1648000000000000, 0.6666666666666666, 0.203000000));
  pos3.push_back(XcVector(0.0959000000000000, 0.7365000000000000, 0.500000000));
  pos3.push_back(XcVector(0.4982000000000000, 0.1648000000000000, 0.203000000));
  pos3.push_back(XcVector(0.3594000000000000, 0.0959000000000000, 0.500000000));
  pos3.push_back(XcVector(0.6666666666666666, 0.5018000000000000, 0.203000000));
  pos3.push_back(XcVector(0.7365000000000000, 0.6405999999999999, 0.500000000));
  pos3.push_back(XcVector(0.8352000000000001, 0.3333333333333333, 0.203000000));
  pos3.push_back(XcVector(0.9041000000000000, 0.2635000000000000, 0.500000000));
  pos3.push_back(XcVector(0.8352000000000001, 0.5018000000000000, 0.797000000));
  pos3.push_back(XcVector(0.4982000000000000, 0.3333333333333333, 0.797000000));
  pos3.push_back(XcVector(0.6666666666666666, 0.1648000000000000, 0.797000000));
  pos3.push_back(XcVector(0.1648000000000000, 0.4982000000000000, 0.797000000));
  pos3.push_back(XcVector(0.5018000000000000, 0.6666666666666666, 0.797000000));
  pos3.push_back(XcVector(0.3333333333333333, 0.8352000000000001, 0.797000000));
  pos3.push_back(XcVector(0.4982000000000000, 0.1648000000000000, 0.797000000));
  pos3.push_back(XcVector(0.6666666666666666, 0.5018000000000000, 0.797000000));
  pos3.push_back(XcVector(0.8352000000000001, 0.3333333333333333, 0.797000000));
  pos3.push_back(XcVector(0.5018000000000000, 0.8352000000000001, 0.797000000));
  pos3.push_back(XcVector(0.3333333333333333, 0.4982000000000000, 0.797000000));
  pos3.push_back(XcVector(0.1648000000000000, 0.6666666666666666, 0.797000000));

  // structure07-niggli.vasp
  XcMatrix cell4(0.00000000000000,   0.00000000000000, 7.56700000000000,
                 9.06300000000000, -15.69757600000000, 0.00000000000000,
                 9.06300000000000,  15.69757600000000, 0.00000000000000);

  pos4.push_back(XcVector(0.50000,  0.26350, 0.90410));
  pos4.push_back(XcVector(0.20300, -0.33340, 0.16480));
  pos4.push_back(XcVector(0.50000, -0.26350, 0.09590));
  pos4.push_back(XcVector(0.20300, -0.16487, 0.50180));
  pos4.push_back(XcVector(0.50000, -0.09590, 0.64060));
  pos4.push_back(XcVector(0.20300, -0.50187, 0.33333));
  pos4.push_back(XcVector(0.50000, -0.64060, 0.26350));
  pos4.push_back(XcVector(0.20300,  0.33340, 0.83520));
  pos4.push_back(XcVector(0.20300,  0.16487, 0.49820));
  pos4.push_back(XcVector(0.50000,  0.09590, 0.35940));
  pos4.push_back(XcVector(0.20300,  0.50187, 0.66667));
  pos4.push_back(XcVector(0.50000,  0.64060, 0.73650));
  pos4.push_back(XcVector(0.20300, -0.33340, 0.50180));
  pos4.push_back(XcVector(0.50000, -0.26350, 0.64060));
  pos4.push_back(XcVector(0.20300, -0.16487, 0.33333));
  pos4.push_back(XcVector(0.50000, -0.09590, 0.26350));
  pos4.push_back(XcVector(0.20300, -0.50187, 0.16480));
  pos4.push_back(XcVector(0.50000, -0.64060, 0.09590));
  pos4.push_back(XcVector(0.20300,  0.33340, 0.49820));
  pos4.push_back(XcVector(0.50000,  0.26350, 0.35940));
  pos4.push_back(XcVector(0.20300,  0.16487, 0.66667));
  pos4.push_back(XcVector(0.50000,  0.09590, 0.73650));
  pos4.push_back(XcVector(0.20300,  0.50187, 0.83520));
  pos4.push_back(XcVector(0.50000,  0.64060, 0.90410));
  pos4.push_back(XcVector(0.79700,  0.33340, 0.83520));
  pos4.push_back(XcVector(0.79700,  0.16487, 0.49820));
  pos4.push_back(XcVector(0.79700,  0.50187, 0.66667));
  pos4.push_back(XcVector(0.79700, -0.33340, 0.16480));
  pos4.push_back(XcVector(0.79700, -0.16487, 0.50180));
  pos4.push_back(XcVector(0.79700, -0.50187, 0.33333));
  pos4.push_back(XcVector(0.79700,  0.33340, 0.49820));
  pos4.push_back(XcVector(0.79700,  0.16487, 0.66667));
  pos4.push_back(XcVector(0.79700,  0.50187, 0.83520));
  pos4.push_back(XcVector(0.79700, -0.33340, 0.50180));
  pos4.push_back(XcVector(0.79700, -0.16487, 0.33333));
  pos4.push_back(XcVector(0.79700, -0.50187, 0.16480));

  // Store pointers to these to simplify test code, since we need to
  // test every combination.
  std::vector<XcMatrix*> cells;
  cells.push_back(&cell1);
  cells.push_back(&cell2);
  cells.push_back(&cell3);
  cells.push_back(&cell4);

  std::vector<std::vector<XcVector>*> positions;
  positions.push_back(&pos1);
  positions.push_back(&pos2);
  positions.push_back(&pos3);
  positions.push_back(&pos4);

  // We need to test every combination, including swapping the order
  bool anyFailures = false;
  for (size_t i = 0; i < 4; ++i) {
    for (size_t j = 0; j < 4; ++j) {
      if (i == j)
        continue;

      bool match = XtalComp::compare(*cells[i], types, *positions[i],
                                     *cells[j], types, *positions[j],
                                     NULL, 0.05, 0.25);

      if (!match) {
        fprintf(stdout, "'%s' failed for comparison: '%lu' with '%lu'\n",
                __FUNCTION__, i + 1, j + 1);
        anyFailures = true;
      }
    }
  }

  if (anyFailures)
    return false;

  return true;
}

int main()
{
  int failures = 0;
  int successes = 0;

  runTest(&simpleCase, "Simple Case", successes, failures);
  runTest(&simpleNiggli, "Simple Niggli", successes, failures);
  runTest(&simpleUniformTranslation, "Simple Uniform Translation",
          successes, failures);
  runTest(&simpleUniformTranslation, "Simple Random Noise (max = 0.005 A)",
          successes, failures);
  runTest(&allOfTheAbove, "All of the above test", successes, failures);
  runTest(&hexagonalCellTest, "Hexagonal cell test", successes, failures);
  runTest(&atomsOverlapping, "Atoms Overlapping Test", successes, failures);
  runTest(&craigFisherTest, "Craig Fisher Test", successes, failures);

  return failures;
}
