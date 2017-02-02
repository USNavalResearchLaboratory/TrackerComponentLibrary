// Msg 21 - ATON status - 'E'

#include "ais.h"

namespace libais {

Ais21::Ais21(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), aton_type(0), position_accuracy(0), dim_a(0),
      dim_b(0), dim_c(0), dim_d(0), fix_type(0), timestamp(0), off_pos(false),
      aton_status(0), raim(false), virtual_aton(false), assigned_mode(false),
      spare(0), spare2(0) {
  assert(message_id == 21);

  // TODO(schwehr): make this more careful than 272-360
  if (num_bits != 268 && (num_bits < 272 || num_bits > 360)) {
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
  aton_type = bs.ToUnsignedInt(38, 5);
  name = bs.ToString(43, 120);
  position_accuracy = bs[163];
  position = bs.ToAisPoint(164, 55);
  dim_a = bs.ToUnsignedInt(219, 9);
  dim_b = bs.ToUnsignedInt(228, 9);
  dim_c = bs.ToUnsignedInt(237, 6);
  dim_d = bs.ToUnsignedInt(243, 6);
  fix_type = bs.ToUnsignedInt(249, 4);
  timestamp = bs.ToUnsignedInt(253, 6);
  off_pos = bs[259];
  aton_status = bs.ToUnsignedInt(260, 8);
  if (num_bits == 268) {
    // Non-standard small message.
    assert(bs.GetRemaining() == 0);
    status = AIS_OK;
    return;
  }
  raim = bs[268];
  virtual_aton = bs[269];
  assigned_mode = bs[270];
  spare = bs[271];

  const size_t extra_chars = bs.GetRemaining() / 6;
  const size_t extra_bits = bs.GetRemaining() % 6;

  if (extra_chars > 0) {
    name += bs.ToString(272, extra_chars * 6);
  }

  if (extra_bits > 0) {
    spare2 = bs.ToUnsignedInt(272 + extra_chars * 6, extra_bits);
  } else {
    spare2 = 0;
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

}  // namespace libais
