// Copyright 2015 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// TODO(schwehr): Add logging back in or remove it.

#include "decode_body.h"

#include <memory>

#include "ais.h"

using libais::AisMsg;
using std::unique_ptr;

namespace libais {

// TODO(schwehr): Switch to make_unique when C++14 is available on Travis-CI.
template <typename T, typename... Args>
std::unique_ptr<T> MakeUnique(Args &&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

unique_ptr<AisMsg> CreateAisMsg6(const string &body, const int fill_bits) {
  libais::Ais6 msg(body.c_str(), fill_bits);
  switch (msg.dac) {
    // International Maritime Organization (IMO).
    case libais::AIS_DAC_1_INTERNATIONAL:
      switch (msg.fi) {
        case 0:
          return MakeUnique<libais::Ais6_1_0>(body.c_str(), fill_bits);
        case 1:
          return MakeUnique<libais::Ais6_1_1>(body.c_str(), fill_bits);
        case 2:
          return MakeUnique<libais::Ais6_1_2>(body.c_str(), fill_bits);
        case 3:
          return MakeUnique<libais::Ais6_1_3>(body.c_str(), fill_bits);
        case 4:
          return MakeUnique<libais::Ais6_1_4>(body.c_str(), fill_bits);
        case 12:
          return MakeUnique<libais::Ais6_1_12>(body.c_str(), fill_bits);
        case 14:
          return MakeUnique<libais::Ais6_1_14>(body.c_str(), fill_bits);
        case 18:
          return MakeUnique<libais::Ais6_1_18>(body.c_str(), fill_bits);
        case 20:
          return MakeUnique<libais::Ais6_1_20>(body.c_str(), fill_bits);
        case 25:
          return MakeUnique<libais::Ais6_1_25>(body.c_str(), fill_bits);
        // TODO(schwehr): 28.
        // TODO(schwehr): 30.
        case 32:
          return MakeUnique<libais::Ais6_1_32>(body.c_str(), fill_bits);
        case 40:
          return MakeUnique<libais::Ais6_1_40>(body.c_str(), fill_bits);
      }
  }
  return nullptr;
}

unique_ptr<AisMsg> CreateAisMsg8(const string &body, const int fill_bits) {
  libais::Ais8 msg(body.c_str(), fill_bits);
  switch (msg.dac) {
    // International Maritime Organization (IMO).
    case libais::AIS_DAC_1_INTERNATIONAL:
      switch (msg.fi) {
        case 0:
          return MakeUnique<libais::Ais8_1_0>(body.c_str(), fill_bits);
        case 11:
          return MakeUnique<libais::Ais8_1_11>(body.c_str(), fill_bits);
        case 13:
          return MakeUnique<libais::Ais8_1_13>(body.c_str(), fill_bits);
        case 15:
          return MakeUnique<libais::Ais8_1_15>(body.c_str(), fill_bits);
        case 16:
          return MakeUnique<libais::Ais8_1_16>(body.c_str(), fill_bits);
        case 17:
          return MakeUnique<libais::Ais8_1_17>(body.c_str(), fill_bits);
        case 19:
          return MakeUnique<libais::Ais8_1_19>(body.c_str(), fill_bits);
        case 21:
          return MakeUnique<libais::Ais8_1_21>(body.c_str(), fill_bits);
        case 22:
          return MakeUnique<libais::Ais8_1_22>(body.c_str(), fill_bits);
        case 24:
          return MakeUnique<libais::Ais8_1_24>(body.c_str(), fill_bits);
        case 26:
          return MakeUnique<libais::Ais8_1_26>(body.c_str(), fill_bits);
        case 27:
          return MakeUnique<libais::Ais8_1_27>(body.c_str(), fill_bits);
        case 29:
          return MakeUnique<libais::Ais8_1_29>(body.c_str(), fill_bits);
        case 31:
          return MakeUnique<libais::Ais8_1_31>(body.c_str(), fill_bits);
      }
    // European River Information System (RIS).
    case libais::AIS_DAC_200_RIS:
      switch (msg.fi) {
        case 10:
          return MakeUnique<libais::Ais8_200_10>(body.c_str(), fill_bits);
        case 23:
          return MakeUnique<libais::Ais8_200_23>(body.c_str(), fill_bits);
        case 24:
          return MakeUnique<libais::Ais8_200_24>(body.c_str(), fill_bits);
        case 40:
          return MakeUnique<libais::Ais8_200_40>(body.c_str(), fill_bits);
        case 55:
          return MakeUnique<libais::Ais8_200_55>(body.c_str(), fill_bits);
      }
    // TODO(schwehr): 366 US Coast Guard.
    case 367:  // US Coast Guard.
      switch (msg.fi) {
        case 22:
          return MakeUnique<libais::Ais8_367_22>(body.c_str(), fill_bits);
      }
  }
  return nullptr;
}

unique_ptr<AisMsg> CreateAisMsg(const string &body, const int fill_bits) {
  if (body.empty()) {
    return nullptr;
  }
  if (fill_bits < 0 || fill_bits > 5) {
    return nullptr;
  }

  switch (body[0]) {
    case '1':  // FALLTHROUGH
    case '2':  // FALLTHROUGH
    case '3':  // 1-3: Class A position report.
      return MakeUnique<libais::Ais1_2_3>(body.c_str(), fill_bits);

    case '4':  // FALLTHROUGH - 4 - Basestation report
    case ';':  // 11 - UTC date response
      return MakeUnique<libais::Ais4_11>(body.c_str(), fill_bits);

    case '5':  // 5 - Ship and Cargo
      return MakeUnique<libais::Ais5>(body.c_str(), fill_bits);

    case '6':  // 6 - Addressed binary message
      return CreateAisMsg6(body, fill_bits);

    case '7':  // FALLTHROUGH - 7 - ACK for addressed binary message
    case '=':  // 13 - ASRM Ack  (safety message)
      return MakeUnique<libais::Ais7_13>(body.c_str(), fill_bits);

    case '8':  // 8 - Binary broadcast message (BBM)
      return CreateAisMsg8(body, fill_bits);

    case '9':  // 9 - SAR Position
      return MakeUnique<libais::Ais9>(body.c_str(), fill_bits);

    case ':':  //  10 - UTC Query
      return MakeUnique<libais::Ais10>(body.c_str(), fill_bits);

    // ';' 11 - See 4

    case '<':  // 12 - Addressed Safety Related Messages (ASRM)
      return MakeUnique<libais::Ais12>(body.c_str(), fill_bits);

    // '=' 13 - See 7

    case '>':  // 14 - Safety Related Broadcast Message (SRBM)
      return MakeUnique<libais::Ais14>(body.c_str(), fill_bits);

    case '?':  // 15 - Interrogation
      return MakeUnique<libais::Ais15>(body.c_str(), fill_bits);

    case '@':  // 16 - Assigned mode command
      return MakeUnique<libais::Ais16>(body.c_str(), fill_bits);

    case 'A':  // 17 - GNSS broadcast
      return MakeUnique<libais::Ais17>(body.c_str(), fill_bits);

    case 'B':  // 18 - Position, Class B
      return MakeUnique<libais::Ais18>(body.c_str(), fill_bits);

    case 'C':  // 19 - Position and ship, Class B
      return MakeUnique<libais::Ais19>(body.c_str(), fill_bits);

    case 'D':  // 20 - Data link management
      return MakeUnique<libais::Ais20>(body.c_str(), fill_bits);

    case 'E':  // 21 - Aids to navigation report
      return MakeUnique<libais::Ais21>(body.c_str(), fill_bits);

    case 'F':  // 22 - Channel Management
      return MakeUnique<libais::Ais22>(body.c_str(), fill_bits);

    case 'G':  // 23 - Group Assignment Command
      return MakeUnique<libais::Ais23>(body.c_str(), fill_bits);

    case 'H':  // 24 - Static data report
      return MakeUnique<libais::Ais24>(body.c_str(), fill_bits);

    case 'I':  // 25 - Single slot binary message
      return MakeUnique<libais::Ais25>(body.c_str(), fill_bits);

    case 'J':  // 26 - Multi slot binary message with comm state
      return MakeUnique<libais::Ais26>(body.c_str(), fill_bits);

    case 'K':  // 27 - Long-range AIS broadcast message
      return MakeUnique<libais::Ais27>(body.c_str(), fill_bits);

    default:
      return nullptr;
  }

  return nullptr;
}

}  // namespace libais
