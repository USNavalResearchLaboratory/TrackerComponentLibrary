// < - ASRM

#include "ais.h"

namespace libais {

Ais12::Ais12(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), seq_num(0), dest_mmsi(0), retransmitted(false),
      spare(0), spare2(0) {

  assert(message_id == 12);

  // WARNING: The ITU 1371 specifications says the maximum number of bits is
  // 1008, but it appears that the maximum should be 1192.
  if (num_bits < 72 || num_bits > 1192)  {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;  // Spec says 1008
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(38);
  seq_num = bs.ToUnsignedInt(38, 2);
  dest_mmsi = bs.ToUnsignedInt(40, 30);
  retransmitted = bs[70];
  spare = bs[71];
  const int num_txt = (num_bits - 72) / 6;
  const int num_txt_bits = num_txt * 6;
  text = bs.ToString(72, num_txt_bits);
  if (bs.GetRemaining() > 0) {
    spare2 = bs.ToUnsignedInt(72 + num_txt_bits, bs.GetRemaining());
  }

  status = AIS_OK;
}

}  // namespace libais
