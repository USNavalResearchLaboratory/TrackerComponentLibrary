// Address Binary Message (ABM) 6

#include <cmath>
#include <iomanip>

#include "ais.h"

namespace libais {

Ais6::Ais6(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), seq(0), mmsi_dest(0), retransmit(false),
      spare(0), dac(0), fi(0) {
  assert(message_id == 6);

  // TODO(olafsinram): 46 or rather 56??
  const int payload_len = num_bits - 46;  // in bits w/o DAC/FI
  if (num_bits < 88 || payload_len < 0 || payload_len > 952) {
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
  seq = bs.ToUnsignedInt(38, 2);
  mmsi_dest = bs.ToUnsignedInt(40, 30);
  retransmit = !bs[70];
  spare = bs[71];
  dac = bs.ToUnsignedInt(72, 10);
  fi = bs.ToUnsignedInt(82, 6);
}

// http://www.e-navigation.nl/content/monitoring-aids-navigation
Ais6_0_0::Ais6_0_0(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad),
      sub_id(1),
      voltage(0.0),
      current(0.0),
      dc_power_supply(true),
      light_on(true),
      battery_low(false),
      off_position(false),
      spare2(0) {
  assert(dac == 0);
  assert(fi == 0);

  if (num_bits != 136) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  sub_id = bs.ToUnsignedInt(88, 16);
  voltage = bs.ToUnsignedInt(104, 12) / 10.0;
  current = bs.ToUnsignedInt(116, 10) / 10.0;
  dc_power_supply = bs[126];
  light_on = bs[127];
  battery_low = bs[128];
  off_position = bs[129];

  spare2 = bs.ToUnsignedInt(130, 6);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

Ais6_1_0::Ais6_1_0(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), ack_required(false), msg_seq(0),
      spare2(0) {

  assert(dac == 1);
  assert(fi == 0);

  if (num_bits < 88 || num_bits > 936) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;  // TODO(schwehr): what is the real max size?
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  ack_required = bs[88];
  msg_seq = bs.ToUnsignedInt(89, 11);

  const size_t text_size = 6 * ((num_bits - 100) / 6);
  const size_t spare2_size = num_bits - 100 - text_size;
  text =  bs.ToString(100, text_size);

  if (!spare2_size)
    spare2 = 0;
  else
    spare2 = bs.ToUnsignedInt(100 + text_size, spare2_size);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

Ais6_1_1::Ais6_1_1(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), ack_dac(0), msg_seq(0), spare2(0) {
  assert(dac == 1);
  assert(fi == 1);

  if (num_bits != 112) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  ack_dac = bs.ToUnsignedInt(88, 10);
  msg_seq = bs.ToUnsignedInt(98, 11);
  spare2 = bs.ToUnsignedInt(109, 3);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

Ais6_1_2::Ais6_1_2(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), req_dac(0), req_fi(0) {
  assert(dac == 1);
  assert(fi == 2);

  if (num_bits != 104) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  req_dac = bs.ToUnsignedInt(88, 10);
  req_fi = bs.ToUnsignedInt(98, 6);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// IFM 3: Capability interrogation - OLD ITU 1371-1
Ais6_1_3::Ais6_1_3(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), req_dac(0), spare2(0) {
  assert(dac == 1);
  assert(fi == 3);

  if (num_bits != 104) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  req_dac = bs.ToUnsignedInt(88, 10);
  spare2 = bs.ToUnsignedInt(94, 6);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// IFM 4: Capability reply - OLD ITU 1371-4
// TODO(schwehr): WTF?  10 + 128 + 6 == 80  Is this 168 or 232 bits?
Ais6_1_4::Ais6_1_4(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), ack_dac(0), capabilities(),
      cap_reserved(), spare2(0) {
  assert(dac == 1);
  assert(fi == 4);

  // TODO(schwehr): num_bits for 6_1_4.  226 bits?
  if (num_bits != 232) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  ack_dac = bs.ToUnsignedInt(88, 10);
  for (size_t cap_num = 0; cap_num < 128/2; cap_num++) {
    size_t start = 98 + cap_num * 2;
    capabilities[cap_num] = bs[start];
    cap_reserved[cap_num] = bs[start + 1];
  }
  // spare2 = bs.ToUnsignedInt(226, 6);  // OR NOT
  // TODO(schwehr): add in the offset of the dest mmsi

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// IMO 1371-5 Ack
Ais6_1_5::Ais6_1_5(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), ack_dac(0), ack_fi(0), seq_num(0),
      ai_available(false), ai_response(0), spare(0) {
  assert(dac == 1);
  assert(fi == 5);

  if (num_bits != 168) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  ack_dac = bs.ToUnsignedInt(88, 10);
  ack_fi = bs.ToUnsignedInt(98, 6);
  seq_num = bs.ToUnsignedInt(104, 11);
  ai_available = static_cast<bool>(bs[115]);
  ai_response = bs.ToUnsignedInt(116, 3);
  spare = bs.ToUnsignedInt(119, 49);

  assert(bs.GetRemaining() == 0);

  status = AIS_OK;
}

// IMO Circ 289 - Dangerous cargo
// See also Circ 236
Ais6_1_12::Ais6_1_12(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), utc_month_dep(0), utc_day_dep(0),
      utc_hour_dep(0), utc_min_dep(0), utc_month_next(0),
      utc_day_next(0), utc_hour_next(0), utc_min_next(0),
      un(0), value(0), value_unit(0), spare2(0) {
  assert(dac == 1);
  assert(fi == 12);

  if (num_bits != 360) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  // TODO(schwehr): add in the offset of the dest mmsi

#if 0
  bs.SeekTo(56);
  last_port = bs.ToString(56, 30);
  utc_month_dep = bs.ToUnsignedInt(86, 4);
  utc_day_dep = bs.ToUnsignedInt(90, 5);
  utc_hour_dep = bs.ToUnsignedInt(95, 5);
  utc_min_dep = bs.ToUnsignedInt(100, 6);
  next_port = bs.ToString(106, 30);
  utc_month_next = bs.ToUnsignedInt(136, 4);  // estimated arrival
  utc_day_next = bs.ToUnsignedInt(140, 5);
  utc_hour_next = bs.ToUnsignedInt(145, 5);
  utc_min_next = bs.ToUnsignedInt(150, 6);
  main_danger = bs.ToString(156, 120);
  imo_cat = bs.ToString(276, 24);
  un = bs.ToUnsignedInt(300, 13);
  value = bs.ToUnsignedInt(313, 10);  // TODO(schwehr): units
  value_unit = bs.ToUnsignedInt(323, 2);
  spare = bs.ToUnsignedInt(325, 3);
  // 360
#endif

  // TODO(schwehr): Add assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// 6_1_13 Does not exist

// IMO Circ 289 - Tidal Window
// See also Circ 236
Ais6_1_14::Ais6_1_14(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), utc_month(0), utc_day(0) {
  // TODO(schwehr): untested - no sample of the correct length yet
  assert(dac == 1);
  assert(fi == 14);

  if (num_bits != 376) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  utc_month = bs.ToUnsignedInt(88, 4);
  utc_day = bs.ToUnsignedInt(92, 5);

  for (size_t window_num = 0; window_num < 3; window_num++) {
    Ais6_1_14_Window w;
    const size_t start = 97 + window_num * 93;
    // Reversed order for lng/lat.
    float y = bs.ToInt(start, 27) / 600000.;
    float x = bs.ToInt(start + 27, 28) / 600000.;
    w.position = AisPoint(x, y);
    w.utc_hour_from = bs.ToUnsignedInt(start + 55, 5);
    w.utc_min_from = bs.ToUnsignedInt(start + 60, 6);
    w.utc_hour_to = bs.ToUnsignedInt(start + 66, 5);
    w.utc_min_to = bs.ToUnsignedInt(start + 71, 6);
    w.cur_dir = bs.ToUnsignedInt(start + 77, 9);
    w.cur_speed  = bs.ToUnsignedInt(start + 86, 7) / 10.;

    windows.push_back(w);
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// IMO Circ 289 - Clearance time to enter port
Ais6_1_18::Ais6_1_18(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), link_id(0), utc_month(0), utc_day(0),
      utc_hour(0), utc_min(0), spare2() {
  assert(dac == 1);
  assert(fi == 18);

  if (num_bits != 360) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  link_id = bs.ToUnsignedInt(88, 10);
  utc_month = bs.ToUnsignedInt(98, 4);
  utc_day = bs.ToUnsignedInt(102, 5);
  utc_hour = bs.ToUnsignedInt(107, 5);
  utc_min = bs.ToUnsignedInt(112, 6);
  port_berth = bs.ToString(118, 120);
  dest = bs.ToString(238, 30);
  position = bs.ToAisPoint(268, 49);
  spare2[0] = bs.ToUnsignedInt(317, 32);
  spare2[1] = bs.ToUnsignedInt(349, 11);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// IMO Circ 289 - Berthing data
Ais6_1_20::Ais6_1_20(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), link_id(0), length(0), depth(0.0),
      mooring_position(0), utc_month(0), utc_day(0), utc_hour(0), utc_min(0),
      services_known(false), services() {
  assert(dac == 1);
  assert(fi == 20);

  if (num_bits != 360) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  link_id = bs.ToUnsignedInt(88, 10);
  length = bs.ToUnsignedInt(98, 9);
  depth = bs.ToUnsignedInt(107, 8);
  mooring_position = bs.ToUnsignedInt(115, 3);
  utc_month = bs.ToUnsignedInt(118, 4);
  utc_day = bs.ToUnsignedInt(122, 5);
  utc_hour = bs.ToUnsignedInt(127, 5);
  utc_min = bs.ToUnsignedInt(132, 6);
  services_known = bs[138];
  for (size_t serv_num = 0; serv_num < 26; serv_num++) {
    // TODO(schwehr): const int val = bs.ToUnsignedInt(139 + 2*serv_num, 2);
    services[serv_num]
        = static_cast<int>(bs.ToUnsignedInt(139 + 2*serv_num, 2));
  }
  name = bs.ToString(191, 120);
  position = bs.ToAisPoint(311, 49);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

Ais6_1_25_Cargo::Ais6_1_25_Cargo()
    : code_type(0), imdg_valid(false), imdg(0), spare_valid(false), spare(0),
      un_valid(false), un(0), bc_valid(false), bc(0), marpol_oil_valid(false),
      marpol_oil(0), marpol_cat_valid(false), marpol_cat(0) {}

// IMO Circ 289 - Dangerous cargo indication 2
// See also Circ 236
Ais6_1_25::Ais6_1_25(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), amount_unit(0), amount(0) {
  assert(dac == 1);
  assert(fi == 25);

  // TODO(schwehr): verify multiple of the size of cargos + header
  //   or padded to a slot boundary
  // Allowing a message with no payloads
  // TODO(schwehr): (num_bits - 100) % 17 != 0) is okay
  if (num_bits < 100 || num_bits > 576) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }
  if ((num_bits - 100) % 17 != 0) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  amount_unit = bs.ToUnsignedInt(88, 2);
  amount = bs.ToUnsignedInt(90, 10);
  const size_t total_cargos = static_cast<int>(floor((num_bits - 100) / 17.));
  for (size_t cargo_num = 0; cargo_num < total_cargos; cargo_num++) {
    Ais6_1_25_Cargo cargo;
    const size_t start = 100 + 17*cargo_num;
    cargo.code_type = bs.ToUnsignedInt(start, 4);

    // TODO(schwehr): Is this the correct behavior?
    switch (cargo.code_type) {
      // No 0
      case 1:  // IMDG Code in packed form
        cargo.imdg = bs.ToUnsignedInt(start + 4, 7);
        cargo.imdg_valid = true;
        cargo.spare = bs.ToUnsignedInt(start + 11, 6);
        cargo.spare_valid = true;
        break;
      case 2:  // IGC Code
        cargo.un = bs.ToUnsignedInt(start + 4, 13);
        cargo.un_valid = true;
        break;
      case 3:  // BC Code
        cargo.bc = bs.ToUnsignedInt(start + 4, 3);
        cargo.bc_valid = true;
        cargo.imdg = bs.ToUnsignedInt(start + 7, 7);
        cargo.imdg_valid = true;
        cargo.spare = bs.ToUnsignedInt(start + 14, 3);
        cargo.spare_valid = true;
        break;
      case 4:  // MARPOL Annex I
        cargo.marpol_oil = bs.ToUnsignedInt(start + 4, 4);
        cargo.marpol_oil_valid = true;
        cargo.spare = bs.ToUnsignedInt(start + 8, 9);
        cargo.spare_valid = true;
        break;
      case 5:  // MARPOL Annex II IBC
        cargo.marpol_cat = bs.ToUnsignedInt(start + 4, 3);
        cargo.marpol_cat_valid = true;
        cargo.spare = bs.ToUnsignedInt(start + 7, 10);
        cargo.spare_valid = true;
        break;
      // 6: Regional use
      // 7: 7-15 reserved for future
      default:
        break;  // Just push in an all blank record?
    }
    cargos.push_back(cargo);
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// TODO(schwehr): 6_1_28 - Modify 8_1_28 once that is debugged

Ais6_1_32_Window::Ais6_1_32_Window()
    : from_utc_hour(0), from_utc_min(0), to_utc_hour(0), to_utc_min(0),
      cur_dir(0), cur_speed(0.0) {}

// IMO Circ 289 - Tidal window
// See also Circ 236
Ais6_1_32::Ais6_1_32(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), utc_month(0), utc_day(0) {
  assert(dac == 1);
  assert(fi == 32);

  // TODO(schwehr): might get messages with not all windows
  if (num_bits != 350) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  utc_month = bs.ToUnsignedInt(88, 4);
  utc_day = bs.ToUnsignedInt(92, 5);

  for (size_t window_num = 0; window_num < 3; window_num++) {
    Ais6_1_32_Window w;
    const size_t start = 97 + 88*window_num;
    w.position = bs.ToAisPoint(start, 49);
    w.from_utc_hour = bs.ToUnsignedInt(start + 49, 5);
    w.from_utc_min = bs.ToUnsignedInt(start + 54, 6);
    w.to_utc_hour = bs.ToUnsignedInt(start + 60, 5);
    w.to_utc_min = bs.ToUnsignedInt(start + 65, 6);
    w.cur_dir = bs.ToUnsignedInt(start + 71, 9);
    w.cur_speed = bs.ToUnsignedInt(start + 80, 8) / 10.;
    windows.push_back(w);
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// IFM 40: people on board - OLD ITU 1371-4
Ais6_1_40::Ais6_1_40(const char *nmea_payload, const size_t pad)
    : Ais6(nmea_payload, pad), persons(0), spare2(0) {
  assert(dac == 1);
  assert(fi == 40);

  if (num_bits != 104) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(88);
  persons = bs.ToUnsignedInt(88, 13);
  spare2 = bs.ToUnsignedInt(101, 3);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

}  // namespace libais
