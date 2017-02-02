// Class A shipdata

#include "ais.h"

namespace libais {

Ais5::Ais5(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), ais_version(0), imo_num(0),
      type_and_cargo(0), dim_a(0), dim_b(0), dim_c(0), dim_d(0),
      fix_type(0), eta_month(0), eta_day(0), eta_hour(0), eta_minute(0),
      draught(0.0), dte(0), spare(0) {

  assert(message_id == 5);

  if (pad != 2 || num_chars != 71) {
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
  ais_version = bs.ToUnsignedInt(38, 2);
  imo_num = bs.ToUnsignedInt(40, 30);
  callsign = bs.ToString(70, 42);

  name = bs.ToString(112, 120);

  type_and_cargo = bs.ToUnsignedInt(232, 8);
  dim_a = bs.ToUnsignedInt(240, 9);
  dim_b = bs.ToUnsignedInt(249, 9);
  dim_c = bs.ToUnsignedInt(258, 6);
  dim_d = bs.ToUnsignedInt(264, 6);
  fix_type = bs.ToUnsignedInt(270, 4);
  eta_month = bs.ToUnsignedInt(274, 4);
  eta_day = bs.ToUnsignedInt(278, 5);
  eta_hour = bs.ToUnsignedInt(283, 5);
  eta_minute = bs.ToUnsignedInt(288, 6);
  draught = bs.ToUnsignedInt(294, 8) / 10.;
  destination = bs.ToString(302, 120);
  dte = bs[422];
  spare = bs[423];

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

ostream& operator<< (ostream& o, const Ais5 &msg) {
  return o << 5 << ": " << msg.mmsi << " \"" << msg.name << "\" "
           << msg.type_and_cargo << " " << msg.dim_a + msg.dim_b
           << "x" << msg.dim_c + msg.dim_d << "x" << msg.draught << "m";
}

}  // namespace libais
