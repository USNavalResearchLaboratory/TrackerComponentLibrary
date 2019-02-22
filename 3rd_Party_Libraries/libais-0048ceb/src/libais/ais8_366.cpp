// USCG binary application messages (BBM).

#include "ais.h"

namespace libais {

// US Coast Guard (USCG) blue force encrypted position report.
//
// Also known as: Blue Force Tracking (BFT) or encrypted AIS (EAIS).
// Messages use AES 256 or Blowfish, but the decryption is not done here.
//
// "NAIS Performance Specification" says:
//
// Technical Characteristics for USCG Encrypted Automatic
// Identification System (EAIS) Very High Frequency Data Link (VDL)
// Standard (v4.0)
//
// (FIPS) 140-2 and 197 for data communications encryption
Ais8_366_56::Ais8_366_56(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad) {
  assert(dac == 366);
  assert(fi == 56);

  if (num_bits < 56 || num_bits > 1192) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(56);
  int num_full_bytes = bs.GetRemaining() / 8;

  for (int i = 0; i < num_full_bytes; i++) {
    encrypted.push_back(bs.ToUnsignedInt(56 + i * 8, 8));
  }

  if (bs.GetRemaining() > 0) {
    assert(bs.GetRemaining() < 8);
    encrypted.push_back(
        bs.ToUnsignedInt(bs.GetPosition(), bs.GetRemaining()));
  }

  status = AIS_OK;
}

}  // namespace libais
