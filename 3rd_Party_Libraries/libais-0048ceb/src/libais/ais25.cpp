// Msg 25 - I - Single slot binary message

// See also: http://www.e-navigation.nl/asm

#include "ais.h"

namespace libais {

Ais25::Ais25(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), use_app_id(false),  dest_mmsi_valid(false),
      dest_mmsi(false), dac(0), fi(0) {
  assert(message_id == 25);

  if (num_bits < 40 || num_bits > 168) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(38);
  const bool addressed = bs[38];
  use_app_id = bs[39];
  if (addressed) {
    dest_mmsi_valid = true;
    dest_mmsi = bs.ToUnsignedInt(40, 30);
    if (use_app_id) {
      dac = bs.ToUnsignedInt(70, 10);
      fi = bs.ToUnsignedInt(80, 6);
    }
    // TODO(schwehr): Handle the payloads.
  } else {
    // broadcast
    if (use_app_id) {
      dac = bs.ToUnsignedInt(40, 10);
      fi = bs.ToUnsignedInt(50, 6);
    }
    // TODO(schwehr): Handle the payloads.
  }

  // TODO(schwehr): Add assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

}  // namespace libais
