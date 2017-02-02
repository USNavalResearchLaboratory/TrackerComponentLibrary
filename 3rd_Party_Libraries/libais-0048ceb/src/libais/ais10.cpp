// : UTC and date query

#include "ais.h"

namespace libais {

Ais10::Ais10(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), spare(0), dest_mmsi(0), spare2(0) {

  assert(message_id == 10);

  if (pad != 0 || num_chars != 12) {
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
  dest_mmsi = bs.ToUnsignedInt(40, 30);
  spare2 = bs.ToUnsignedInt(70, 2);

  assert(bs.GetRemaining() == 0);

  status = AIS_OK;
}

ostream& operator<< (ostream &o, const Ais10 &msg) {
  return o << msg.message_id << ": " << msg.mmsi
           << " dest=" << msg.dest_mmsi
           << " " << msg.spare << " " << msg.spare2;
}

}  // namespace libais
