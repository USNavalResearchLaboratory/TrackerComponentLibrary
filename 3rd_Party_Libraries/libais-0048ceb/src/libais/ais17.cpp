// A - GNSS broacast -
// TODO(schwehr): only partially coded - need to finish
// http://www.itu.int/rec/R-REC-M.823/en
// http://www.iala-aism.org/iala/publications/documentspdf/doc_348_eng.pdf

// In 823, 30 bit words = 24 bits data followed by 6 parity bits.
// Parity bits are left out of the AIS payload?

#include "ais.h"

namespace libais {

Ais17::Ais17(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), spare(0), spare2(0), gnss_type(0), z_cnt(0),
      station(0), seq(0), health(0) {
  assert(message_id == 17);

  if (num_bits != 80 && (num_bits < 120 || num_bits > 816)) {
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

  position = bs.ToAisPoint(40, 35);
  spare2 = bs.ToUnsignedInt(75, 5);

  // Spec states that there might be no data.
  if (num_bits == 80) {
    gnss_type = -1;
    station = -1;
    z_cnt = -1;
    seq = -1;
    return;
  }

  gnss_type = bs.ToUnsignedInt(80, 6);
  station = bs.ToUnsignedInt(86, 10);
  z_cnt = bs.ToUnsignedInt(96, 13);
  seq = bs.ToUnsignedInt(109, 3);
  bs.SeekRelative(5);
  health = bs.ToUnsignedInt(117, 3);

  // TODO(schwehr): Implement parsing the payload.

  // TODO(schwehr): Add assert(bs.GetRemaining() == 0);

  status = AIS_OK;  // TODO(schwehr): Not really okay yet.
}


ostream& operator<< (ostream &o, const Ais17 &m) {
    return o << "[" << m.message_id << "]: " << m.mmsi
             << " " << m.position << " t:"
             << m.gnss_type << ", z:" << m.z_cnt
             << ", d s:" << m.station << ", seq:"
             << m.seq << ", h:" << m.health;
}

}  // namespace libais
