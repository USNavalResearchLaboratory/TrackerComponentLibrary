// Msg 23 - G - Channel Management

#include "ais.h"

namespace libais {

Ais23::Ais23(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), spare(0), station_type(0), type_and_cargo(0),
      spare2(3), txrx_mode(0), interval_raw(0), quiet(0), spare3(0) {

  assert(message_id == 23);

  if (pad != 2 || num_chars != 27) {
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
  spare = bs.ToUnsignedInt(38, 2);

  position1 = bs.ToAisPoint(40, 35);
  position2 = bs.ToAisPoint(75, 35);

  station_type = bs.ToUnsignedInt(110, 4);
  type_and_cargo = bs.ToUnsignedInt(114, 8);
  spare2 = bs.ToUnsignedInt(122, 22);

  txrx_mode = bs.ToUnsignedInt(144, 2);
  interval_raw = bs.ToUnsignedInt(146, 4);
  quiet = bs.ToUnsignedInt(150, 4);
  spare3 = bs.ToUnsignedInt(154, 6);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

}  // namespace libais
