// @ - Assigned mode command
// TODO(schwehr): Use valid flags rather than negative numbers.

#include "ais.h"

namespace libais {

Ais16::Ais16(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), spare(0), dest_mmsi_a(0), offset_a(0), inc_a(0),
      dest_mmsi_b(0), offset_b(0), inc_b(0), spare2(0) {
  assert(message_id == 16);

  // 96 or 144 bits
  // 168 bits violates the spec but is common
  // TODO(schwehr): check the pad
  if (num_chars != 16 && num_chars != 24 && num_chars != 28) {
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

  dest_mmsi_a = bs.ToUnsignedInt(40, 30);
  offset_a = bs.ToUnsignedInt(70, 12);
  inc_a = bs.ToUnsignedInt(82, 10);
  if (num_chars == 16) {
    dest_mmsi_b = -1;
    offset_b = -1;
    inc_b = -1;
    spare2 = bs.ToUnsignedInt(92, 4);

    assert(bs.GetRemaining() == 0);
    status = AIS_OK;
    return;
  }

  dest_mmsi_b = bs.ToUnsignedInt(92, 30);
  offset_b = bs.ToUnsignedInt(122, 12);
  inc_b = bs.ToUnsignedInt(134, 10);
  spare2 = -1;

  // Currently crashes with the check.
  // TODO(schwehr): Add assert(bs.GetRemaining() == 0);

  status = AIS_OK;
}

}  // namespace libais
