// Msg 19 C - Extended Class B equipment position report

#include "ais.h"

namespace libais {

Ais19::Ais19(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), spare(0), sog(0.0), position_accuracy(0),
      cog(0.0), true_heading(0), timestamp(0), spare2(0), type_and_cargo(0),
      dim_a(0), dim_b(0), dim_c(0), dim_d(0), fix_type(0), raim(false), dte(0),
      assigned_mode(0), spare3(0) {
  assert(message_id == 19);

  if (pad != 0 || num_chars != 52) {
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
  spare = bs.ToUnsignedInt(38, 8);
  sog = bs.ToUnsignedInt(46, 10) / 10.;

  position_accuracy = bs[56];
  position = bs.ToAisPoint(57, 55);

  cog = bs.ToUnsignedInt(112, 12) / 10.;
  true_heading = bs.ToUnsignedInt(124, 9);
  timestamp = bs.ToUnsignedInt(133, 6);
  spare2 = bs.ToUnsignedInt(139, 4);

  name = bs.ToString(143, 120);

  type_and_cargo = bs.ToUnsignedInt(263, 8);
  dim_a = bs.ToUnsignedInt(271, 9);
  dim_b = bs.ToUnsignedInt(280, 9);
  dim_c = bs.ToUnsignedInt(289, 6);
  dim_d = bs.ToUnsignedInt(295, 6);

  fix_type = bs.ToUnsignedInt(301, 4);
  raim = bs[305];
  dte = bs[306];
  assigned_mode = bs[307];
  spare3 = bs.ToUnsignedInt(308, 4);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

}  // namespace libais
