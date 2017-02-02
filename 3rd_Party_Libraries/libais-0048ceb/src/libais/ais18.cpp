// Class B position report - 18 "B"

#include "ais.h"

namespace libais {

Ais18::Ais18(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad),
      spare(0),
      sog(0.0),
      position_accuracy(0),
      cog(0.0),
      true_heading(0),
      timestamp(0),
      spare2(0),
      unit_flag(0),
      display_flag(0),
      dsc_flag(0),
      band_flag(0),
      m22_flag(0),
      mode_flag(0),
      raim(false),
      commstate_flag(0),
      sync_state(0),
      slot_timeout_valid(false),
      slot_timeout(0),
      received_stations_valid(false),
      received_stations(0),
      slot_number_valid(false),
      slot_number(0),
      utc_valid(false),
      utc_hour(0),
      utc_min(0),
      utc_spare(0),
      slot_offset_valid(false),
      slot_offset(0),
      slot_increment_valid(false),
      slot_increment(0),
      slots_to_allocate_valid(false),
      slots_to_allocate(0),
      keep_flag_valid(false),
      keep_flag(0),
      commstate_cs_fill_valid(false),
      commstate_cs_fill(0) {
  assert(message_id == 18);

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

  bs.SeekTo(38);
  spare = bs.ToUnsignedInt(38, 8);
  sog = bs.ToUnsignedInt(46, 10) / 10.;

  position_accuracy = bs[56];
  position = bs.ToAisPoint(57, 55);

  cog = bs.ToUnsignedInt(112, 12) / 10.;
  true_heading = bs.ToUnsignedInt(124, 9);
  timestamp = bs.ToUnsignedInt(133, 6);
  spare2 = bs.ToUnsignedInt(139, 2);
  unit_flag = bs[141];
  display_flag = bs[142];
  dsc_flag = bs[143];
  band_flag = bs[144];
  m22_flag = bs[145];
  mode_flag = bs[146];
  raim = bs[147];
  commstate_flag = bs[148];  // 0 SOTDMA, 1 ITDMA

  if (unit_flag == 0) {
    sync_state = bs.ToUnsignedInt(149, 2);
    if (commstate_flag == 0) {
      // SOTDMA
      slot_timeout = bs.ToUnsignedInt(151, 3);
      slot_timeout_valid = true;

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

    } else {
      // ITDMA
      slot_increment = bs.ToUnsignedInt(151, 13);
      slot_increment_valid = true;

      slots_to_allocate = bs.ToUnsignedInt(164, 3);
      slots_to_allocate_valid = true;

      keep_flag = bs[167];
      keep_flag_valid = true;
    }
  } else {
    // Carrier Sense (CS) with unit_flag of 1.
    commstate_cs_fill = bs.ToUnsignedInt(149, 19);
    commstate_cs_fill_valid = true;
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

}  // namespace libais
