// 'K' - 27 - Long-range AIS broadcast message

#include "ais.h"

namespace libais {

Ais27::Ais27(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), position_accuracy(0), raim(false),
      nav_status(0), sog(0), cog(0), gnss(false), spare(0) {
  assert(message_id == 27);

  if (pad != 0 || num_bits != 96) {
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
  position_accuracy = bs[38];
  raim = bs[39];
  nav_status = bs.ToUnsignedInt(40, 4);
  position = bs.ToAisPoint(44, 35);
  sog = bs.ToUnsignedInt(79, 6);  // Knots.
  cog = bs.ToUnsignedInt(85, 9);  // Degrees.
  // 0 is a current GNSS position.  1 is NOT the current GNSS position
  gnss = !bs[94];
  spare = bs[95];

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

}  // namespace libais
