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

// Parse NMEA AIS VDM strings without extra metadata.
//
// Handles a sequence of AIS VDM NMEA messages.  Messages are composed of
// 1 to 10 lines.  These lines are connected by a sequence number if there
// are more than 1 line required to make a message.  NMEA AIS VDM lines are of
// the following form as a pseudo regular expression, where each line has 7
// comma separated fields:
//
//   ![A-Z]{5},[1-9],[1-9],[0-9]?,[AB],armored_data,[0-5]\*[0-9A-F][0-9A-F]
//
// The fields are:
//
//   0. Start bang, 2 character talker and 3 character sentence.
//   1. Total number of sentences.
//   2. Current sentence number.
//   3. Sequence number.
//   4. VHF radio channel.
//   5. 6-bit per character encoded binary data (not Base64 encoding).  The
//      first character maps to the AIS message type.
//   6. The number of fill bits included in the last character, a star, and
//      a 2 character hex xor checksum.
//
// An example single line message:
//
//   !AIVDM,1,1,,A,14VIk0002sMM04vE>V9jGimn08RP,0*0D
//
// A multi-line message:
//
// clang-format off
//   !SAVDM,2,1,4,B,55NGlfP00001L@GO??0lU=>0hUaaV2222222220O1p>454wV07PhE82Dk0CQ,0*64
//   // NOLINT
//   !SAVDM,2,2,4,B,84i@H0SmFH0,2*54
// clang-format on
//
// The library is structured around NmeaSentence, which manages each line of
// text.  It validates the checksum and separates the field components.
//
// VdmStream takes a series of sentences, converts them to NmeaScentences, and
// then assembles the parts into complete messages.  When it has all of the
// constituent sentences, it uses libais to decode the armored 6-bit payload
// into AisMsg instances.
//
// The VdmStream is not thread safe.
//
// See Also:
//   http://catb.org/gpsd/AIVDM.html
//   http://www.itu.int/rec/R-REC-M.1371/en

#ifndef LIBAIS_VDM_H_
#define LIBAIS_VDM_H_

#include <deque>
#include <memory>
#include <string>
#include <vector>

// #include "base/logging.h"
#include "ais.h"

namespace libais {

// No more than 10 sentences per message.  The total is [1..9].
// The most observed in real data has been 4.
static const int kMaxSentences = 10;

// AIS receivers use the sequence to group multi-line messages.
static const int kNumSequenceChannels = 10;

// Computes the xor checksum of a string.
uint8_t Checksum(const string &content);
// Convert a number to a two character hex string.
string ToHex2(int32_t val);
// Returns the 2 upper case character xor checksum of a string.
string ChecksumHexString(const string &base);
bool ValidateChecksum(const string &line);

// Manages single lines of NMEA AIS VDM text.
class NmeaSentence {
 public:
  NmeaSentence(const string &talker, const string &sentence_type,
               int sentence_total, int sentence_number, int sequence_number,
               char channel, const string &body, int fill_bits,
               int64_t line_number)
      : talker_(talker),
        sentence_type_(sentence_type),
        sentence_total_(sentence_total),
        sentence_number_(sentence_number),
        sequence_number_(sequence_number),
        channel_(channel),
        body_(body),
        fill_bits_(fill_bits),
        line_number_(line_number) {}

  // Returns a smart pointer for an instance of a NmeaSentence created from a
  // line of text.  Returns nullptr on failure.  The line_number argument
  // tracks the source location in a file or the count of lines pushed through
  // a channel.
  static std::unique_ptr<NmeaSentence> Create(const string &line,
                                              int64_t line_number);

  // Returns a composite NmeaSentence from multiple prior sentences that make up
  // a multi-line message.  The resulting sentence is such that the new payload
  // can be decoded as one 6-bit encoded unit.  Sets the sentence total and
  // sentence number both to 1 so that it looks like a single line sentence.
  // Returns a nullptr if the merge failed.
  // Merging sentences from different messages or out of order messages will
  // fail.
  std::unique_ptr<NmeaSentence> Merge(
      const vector<std::unique_ptr<NmeaSentence>> &prior_sentences) const;

  // Reconstructs the NMEA AIS VDM sentence representation of this instance.
  string ToString() const;
  // Returns the hex digest of the 6-bit encoded body.  Used for detection of
  // multiple instances of a line in an input stream.
  // TODO(schwehr): Consider mixing in the VHF channel into the text.
  string ToMd5Digest() const;

  // Returns true if another sentence is a part of the same multiline message.
  // The given sentence must have the same channel, sequence number and total
  // number of sentences as the instance sentence to be a part of the same AIS
  // NMEA multi-line message.
  //
  // In the same message:
  //   !SAVDM,2,1,1,A,54a...
  //   !SAVDM,2,2,1,A,888...
  //
  // Different channels:
  //   !SAVDM,2,1,1,B,54a...
  //   !SAVDM,2,2,1,A,888...
  //
  // Different sequence:
  //   !SAVDM,2,1,8,A,54a...
  //   !SAVDM,2,2,2,A,888...
  //
  // Different totals:
  //   !SAVDM,2,1,1,A,54a...
  //   !SAVDM,3,2,1,A,888...
  //
  // TODO(schwehr): What should this return when called with the same sentence?
  //   !SAVDM,3,2,1,A,54a...
  //   !SAVDM,3,3,1,A,54a...
  bool VerifyInSameMessage(const NmeaSentence &sentence) const;

  // Accessors to get each of the fields or sub-fields of the sentence.
  string talker() const { return talker_; }
  string sentence_type() const { return sentence_type_; }
  int sentence_total() const { return sentence_total_; }
  int sentence_number() const { return sentence_number_; }
  int sequence_number() const { return sequence_number_; }
  char channel() const { return channel_; }
  string body() const { return body_; }
  int fill_bits() const { return fill_bits_; }

  int64_t line_number() const { return line_number_; }

 private:
  // Which device type sent the message.  AI is an AIS transceiver.
  const string talker_;
  // This will be VDM for an AIS message received and VDO for "own ship."
  const string sentence_type_;
  // Total number of lines that compose the message this is a part of.
  const int sentence_total_;
  // Which sentence is this?  Ranges from 1 to sentence_total_.
  const int sentence_number_;
  // Logical channel for the AIS receiver.  A receiver will put a multi-line
  // message on a single channel and then move to the next channel.  Ranges from
  // 0 to 9, but most receivers will not use 0.
  const int sequence_number_;  // -1 if not valid.
  // VHF channel.  Will be either A or B.  Some receivers do not record this.
  const char channel_;
  // The VDM 6-bit encoded/armored payload of the message.
  const string body_;
  // The number of padding.  Should be zero for all messages except the last.
  // Must be 0, 2, or 4.
  const int fill_bits_;

  const int64_t line_number_;
};

// This class processes a sequence of lines to find the AIS messages across
// lines.  AIS messages come in groups of 1 or more lines.  Its job is
// to return decoded AIS messages as libais::AisMsg instances as they are found
// in the stream.
// Sentences are assembled into messages and output in FIFO order.
class VdmStream {
 public:
  VdmStream() : line_number_(0), incoming_sentences_(kNumSequenceChannels) {}

  // Returns true if the sentence was used or false if the line was ignored.
  // A line will be ignored if it is not a valid VDM line or if it is a later
  // part of a multi-line message but missing one or more initial lines.
  bool AddLine(const string &line);  // Was push
  // Returns nullptr if there are not decoded messages currently available.
  std::unique_ptr<libais::AisMsg> PopOldestMessage();

  int size() const { return messages_.size(); }
  bool empty() const { return messages_.empty(); }

 private:
  // Line number starts at 0 and is incremented to 1 with the first line.
  int64_t line_number_;

  // Decoded messages ready for pickup.
  std::deque<std::unique_ptr<libais::AisMsg>> messages_;
  // Sentences for each sequence number that have yet to get all the required
  // parts to be complete.
  vector<vector<std::unique_ptr<NmeaSentence>>> incoming_sentences_;
};

}  // namespace libais

#endif  // LIBAIS_VDM_H_
