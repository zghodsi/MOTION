// MIT License
//
// Copyright (c) 2019 Oleksandr Tkachenko
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

#include "matmult.h"

#include "algorithm/algorithm_description.h"
#include "protocols/share_wrapper.h"
#include "secure_type/secure_unsigned_integer.h"
#include "statistics/analysis.h"
#include "statistics/run_time_statistics.h"
#include "utility/config.h"

encrypto::motion::RunTimeStatistics EvaluateProtocol(
        encrypto::motion::PartyPointer& party, std::size_t number_of_simd, std::size_t bit_size,
        std::size_t dim, encrypto::motion::MpcProtocol protocol) {
    const std::vector<encrypto::motion::BitVector<>> temporary_bool(
            bit_size, encrypto::motion::BitVector<>(number_of_simd));

    encrypto::motion::SecureUnsignedInteger a, b;
    std::vector<encrypto::motion::SecureUnsignedInteger> v1;
    std::vector<encrypto::motion::SecureUnsignedInteger> v2;
    std::vector<encrypto::motion::SecureUnsignedInteger> tmp;
    std::vector<encrypto::motion::SecureUnsignedInteger> v;

    switch (protocol) {
        case encrypto::motion::MpcProtocol::kBooleanGmw: {
             for (int i=0; i<dim*dim; i++) {
                 v1.push_back(party->In<encrypto::motion::MpcProtocol::kBooleanGmw>(temporary_bool, 0));
                 v2.push_back(party->In<encrypto::motion::MpcProtocol::kBooleanGmw>(temporary_bool, 0));
             }
             break;
        }
        case encrypto::motion::MpcProtocol::kBmr: {
            for (int i=0; i<dim*dim; i++) {
                v1.push_back(party->In<encrypto::motion::MpcProtocol::kBmr>(temporary_bool, 0));
                v2.push_back(party->In<encrypto::motion::MpcProtocol::kBmr>(temporary_bool, 0));
            }
            break;
        }
        case encrypto::motion::MpcProtocol::kArithmeticGmw: {
            switch (bit_size) {
                case 8: {
                    std::vector<std::uint8_t> temporary_arithmetic(number_of_simd);
                    for (int i=0; i<dim*dim; i++) {
                        v1.push_back(party->In<encrypto::motion::MpcProtocol::kArithmeticGmw>(temporary_arithmetic, 0));
                        v2.push_back(party->In<encrypto::motion::MpcProtocol::kArithmeticGmw>(temporary_arithmetic, 0));
                    }
                    break;
                }
                case 16: {
                    std::vector<std::uint16_t> temporary_arithmetic(number_of_simd);
                    for (int i=0; i<dim*dim; i++) {
                        v1.push_back(party->In<encrypto::motion::MpcProtocol::kArithmeticGmw>(temporary_arithmetic, 0));
                        v2.push_back(party->In<encrypto::motion::MpcProtocol::kArithmeticGmw>(temporary_arithmetic, 0));
                    }
                    break;
                }
                case 32: {
                    std::vector<std::uint32_t> temporary_arithmetic(number_of_simd);
                    for (int i=0; i<dim*dim; i++) {
                        v1.push_back(party->In<encrypto::motion::MpcProtocol::kArithmeticGmw>(temporary_arithmetic, 0));
                        v2.push_back(party->In<encrypto::motion::MpcProtocol::kArithmeticGmw>(temporary_arithmetic, 0));
                    }
                    break;
                }
                case 64: {
                    std::vector<std::uint64_t> temporary_arithmetic(number_of_simd);
                    for (int i=0; i<dim*dim; i++) {
                        v1.push_back(party->In<encrypto::motion::MpcProtocol::kArithmeticGmw>(temporary_arithmetic, 0));
                        v2.push_back(party->In<encrypto::motion::MpcProtocol::kArithmeticGmw>(temporary_arithmetic, 0));
                    }
                    break;
                }   
                default:
                    throw std::invalid_argument("Unknown bit size");
            }
            break;
        }
        default: {
            throw std::invalid_argument("Invalid MPC protocol");
        }
    }

    for (int i=0; i<dim; i++) {
        for (int j=0; j<dim; j++) {
            for (int k=0; k<dim; k++) {
                // tmp <- m1[i][k]*m2[k][j] = v1[i*dim+k]*v2[k*dim+j]
                tmp.push_back(v1[i*dim+k]*v2[k*dim+j]);
            }
            // add all tmp values
            v.push_back(tmp.back());
            tmp.pop_back();
            for (int k=0; k<dim-1; k++) {
                v.back() = v.back()+tmp.back();
                tmp.pop_back();
            }
            if (tmp.size()>0) 
                throw std::invalid_argument("Vector tmp not zero");
        }
    }
    if (v.size()!=dim*dim) 
        throw std::invalid_argument("Vector v size not correct");

    party->Run();
    party->Finish();
    const auto& statistics = party->GetBackend()->GetRunTimeStatistics();
    return statistics.front();
}
