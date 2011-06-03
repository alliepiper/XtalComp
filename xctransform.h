/**********************************************************************
  XcTransform - Transformation class for transforming XcVectors

  WARNING: This is not your typical transform class -- it has been
  specialized for XtalComp. It stores the rotation and translation
  separately, and applies the translation followed by the rotation
  when multiplied by a vector

  Copyright (C) 2011 by David C. Lonie

  This source code is released under the New BSD License, (the "License").

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
 ***********************************************************************/

#ifndef XCTRANSFORM_H
#define XCTRANSFORM_H

class XcTransform
{
 public:
  XcTransform() {};
  XcTransform(const XcTransform &other) : rot(other.rot), trans(other.trans) {}

  XcTransform & setIdentity()
  {
    this->rot.fillFromScalar(1.0);
    this->trans.fill(0.0);
    return *this;
  }

  XcTransform &    rotate(const XcMatrix &mat) {this->rot *= mat;}
  XcTransform & prerotate(const XcMatrix &mat) {this->rot = mat * this->rot;}

  XcTransform & translate(const XcVector &vec) {this->trans += vec;}

  XcVector & translation() {return this->trans;}
  const XcVector & translation() const {return this->trans;}

  XcMatrix & rotation() {return this->rot;}
  const XcMatrix & rotation() const {return this->rot;}

  XcVector operator*(const XcVector & vec) const
  {
    return this->rot * (vec + this->trans);
  }

 private:
  XcMatrix rot;
  XcVector trans;
};

#endif
