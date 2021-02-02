// MIT License
//
// Copyright (c) 2021 Oleksandr Tkachenko
// Cryptography and Privacy Engineering Group (ENCRYPTO)
// TU Darmstadt, Germany
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <span>
#include <vector>

#include "protocols/gate.h"

namespace encrypto::motion {

class Share;
using SharePointer = std::shared_ptr<Share>;

class ShareWrapper;

/// \brief yields a set of shares that correspond to single "SIMD layers" of the parent, e.g., if
/// parent contains 10 SIMD values, the output will contain 10 shares with 1 SIMD value and the
/// same number of wires as parent.
class UnsimdifyGate : public OneGate {
 public:
  UnsimdifyGate(const SharePointer& parent);

  ~UnsimdifyGate() final = default;

  void EvaluateSetup() final override;

  void EvaluateOnline() final override;

  std::vector<SharePointer> GetOutputAsVectorOfShares();

  UnsimdifyGate() = delete;

  UnsimdifyGate(const Gate&) = delete;
};

}  // namespace encrypto::motion