// Class B static data report.  Can be one of 4 different parts.  Only
// A and B defined for ITU 1371-3

#include "ais.h"

namespace libais {

Ais24::Ais24(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), part_num(0), type_and_cargo(0),
      dim_a(0), dim_b(0), dim_c(0), dim_d(0), spare(0) {

  assert(message_id == 24);

  if (num_bits != 160 && num_bits != 168) {
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
  part_num = bs.ToUnsignedInt(38, 2);

  switch (part_num) {
  case 0:  // Part A
    name = bs.ToString(40, 120);
    if (num_bits == 168) {
      // Accept the invalid size.
      spare = bs.ToUnsignedInt(160, 8);
    }
    break;
  case 1:  // Part B
    if (num_bits == 160) {
      // Some devices incorrectly use part 1 as 0.
      name = bs.ToString(40, 120);
      part_num = 0;
      break;
    }
    type_and_cargo = bs.ToUnsignedInt(40, 8);
    vendor_id = bs.ToString(48, 42);
    callsign = bs.ToString(90, 42);
    dim_a = bs.ToUnsignedInt(132, 9);
    dim_b = bs.ToUnsignedInt(141, 9);
    dim_c = bs.ToUnsignedInt(150, 6);
    dim_d = bs.ToUnsignedInt(156, 6);
    spare = bs.ToUnsignedInt(162, 6);
    break;
  case 2:  // FALLTHROUGH - Not defined by ITU 1371-5
  case 3:  // FALLTHROUGH - Not defined by ITU 1371-5
  default:
    status = AIS_ERR_BAD_MSG_CONTENT;
    return;
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

}  // namespace libais
