// AIS message 4 or 11

#include "ais.h"

namespace libais {

Ais4_11::Ais4_11(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), year(0), month(0), day(0), hour(0), minute(0),
      second(0), position_accuracy(0), fix_type(0),
      transmission_ctl(0), spare(0), raim(false), sync_state(0),
      slot_timeout(0), received_stations_valid(false), received_stations(0),
      slot_number_valid(false), slot_number(0), utc_valid(false), utc_hour(0),
      utc_min(0), utc_spare(0), slot_offset_valid(false), slot_offset(0) {
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

  assert(message_id == 4 || message_id == 11);

  bs.SeekTo(38);
  year = bs.ToUnsignedInt(38, 14);
  month = bs.ToUnsignedInt(52, 4);
  day = bs.ToUnsignedInt(56, 5);
  hour = bs.ToUnsignedInt(61, 5);
  minute = bs.ToUnsignedInt(66, 6);
  second = bs.ToUnsignedInt(72, 6);

  position_accuracy = bs[78];
  position = bs.ToAisPoint(79, 55);

  fix_type = bs.ToUnsignedInt(134, 4);
  transmission_ctl = bs[138];
  spare = bs.ToUnsignedInt(139, 9);
  raim = bs[148];

  // SOTDMA commstate
  sync_state = bs.ToUnsignedInt(149, 2);
  slot_timeout = bs.ToUnsignedInt(151, 3);

  switch (slot_timeout) {
  case 0:
    slot_offset = bs.ToUnsignedInt(154, 14);
    slot_offset_valid = true;
    break;
  case 1:
    utc_hour = bs.ToUnsignedInt(154, 5);
    utc_min = bs.ToUnsignedInt(159, 7);
    utc_spare = bs.ToUnsignedInt(166, 2);
    utc_valid = true;
    break;
  case 2:  // FALLTHROUGH
  case 4:  // FALLTHROUGH
  case 6:
    slot_number = bs.ToUnsignedInt(154, 14);
    slot_number_valid = true;
    break;
  case 3:  // FALLTHROUGH
  case 5:  // FALLTHROUGH
  case 7:
    received_stations = bs.ToUnsignedInt(154, 14);
    received_stations_valid = true;
    break;
  default:
    assert(false);
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

ostream& operator<< (ostream &o, const Ais4_11 &msg) {
  return o << msg.message_id << ": " << msg.mmsi;
}

}  // namespace libais
