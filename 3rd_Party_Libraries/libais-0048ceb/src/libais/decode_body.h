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

// Decode Automatic Identification System (AIS) messages.

#ifndef LIBAIS_DECODE_BODY_H_
#define LIBAIS_DECODE_BODY_H_

#include <memory>
#include <string>

#include "ais.h"

namespace libais {

// Decodes the payload of an AIS message and returns an AisMsg instance.
// Returns a nullptr on failure.
// The body is the armored text from 1 or more sentences that compose
// the encoded bits for an AIS message.
// The fill_bits are the number of pad bits in the last character of the
// body.  AIS messages are 8-bit aligned and the characters in the armored
// body are 6-bit aligned.
std::unique_ptr<libais::AisMsg> CreateAisMsg(const string &body,
                                             const int fill_bits);

}  // namespace libais

#endif  // LIBAIS_DECODE_BODY_H_
