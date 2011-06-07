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

void runTest(bool (*testFunc)(), const char * testName, int &successes, int &failures)
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

  bool match = XtalComp::XtalComp::compare(cell1, types1, pos1, cell2, types2, pos2, 0.05, 0.25);
  if (!match)
    return false;

  // Displace an atom, ensure that comparison fails.
  pos2[0] += XcVector(0.5,0,0);
  match = XtalComp::XtalComp::compare(cell1, types1, pos1, cell2, types2, pos2, 0.05, 0.25);
  if (match)
    return false;

  return true;
}

int main()
{
  int failures = 0;
  int successes = 0;

  runTest(&simpleCase, "Simple Case", successes, failures);

  return failures;
}
