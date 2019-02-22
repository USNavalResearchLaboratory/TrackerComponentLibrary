// Msg 22 - F - Channel Management

#include "ais.h"

namespace libais {

Ais22::Ais22(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), spare(0), chan_a(0), chan_b(0), txrx_mode(0),
      power_low(false), pos_valid(false), dest_valid(false), dest_mmsi_1(0),
      dest_mmsi_2(0), chan_a_bandwidth(0), chan_b_bandwidth(0), zone_size(0),
      spare2(0) {

  assert(message_id == 22);

  if (pad != 0 || num_chars != 28) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  spare = bs.SeekTo(38).ToUnsignedInt(38, 2);

  chan_a = bs.ToUnsignedInt(40, 12);
  chan_b = bs.ToUnsignedInt(52, 12);
  txrx_mode = bs.ToUnsignedInt(64, 4);
  power_low = bs[68];

  // WARNING: OUT OF ORDER DECODE
  bs.SeekTo(139);
  bool addressed = bs[139];

  bs.SeekTo(69);
  if (!addressed) {
    // geographic position
    // TODO(schwehr): For all implementations in libais, set the valid flag
    //   after setting all of the data members.
    pos_valid = true;
    // TODO(schwehr): Confirm this is correct!
    position1 = bs.ToAisPoint(69, 35);
    position2 = bs.ToAisPoint(104, 35);
  } else {
    dest_valid = true;
    dest_mmsi_1 = bs.ToUnsignedInt(69, 30);
    // 5 spare bits
    bs.SeekRelative(5);
    dest_mmsi_2 = bs.ToUnsignedInt(104, 30);
    // 5 spare bits
    bs.SeekRelative(5);
  }

  // OUT OF ORDER: addressed is earlier before
  bs.SeekRelative(1);
  chan_a_bandwidth = bs[140];
  chan_b_bandwidth = bs[141];
  zone_size = bs.ToUnsignedInt(142, 3);

  spare2 = bs.ToUnsignedInt(145, 23);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

}  // namespace libais
