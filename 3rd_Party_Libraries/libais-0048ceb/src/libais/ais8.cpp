// Binary Broadcast Message (BBM) - 8

// See also: http://www.e-navigation.nl/asm

#include <iomanip>

#include "ais.h"

namespace libais {

Ais8::Ais8(const char *nmea_payload, const size_t pad)
    : AisMsg(nmea_payload, pad), spare(0), dac(0), fi(0) {
  assert(message_id == 8);

  // in bits w/o DAC/FI
  // TODO(schwehr): Verify if this is 46 or 56 (accumulated bits below).
  const int payload_len = num_bits - 56;
  if (payload_len < 0 || payload_len > 952) {
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
  dac = bs.ToUnsignedInt(40, 10);
  fi = bs.ToUnsignedInt(50, 6);
}

Ais8_1_0::Ais8_1_0(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), ack_required(false), msg_seq(0), text(),
      spare2(0) {
  assert(dac == 1);
  assert(fi == 0);

  if (num_bits < 56 || num_bits > 1024) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(56);
  ack_required = bs[56];
  msg_seq = bs.ToUnsignedInt(57, 11);

  const size_t text_size = 6 * ((num_bits - 68)/6);
  // wrong?  needs to land on 8-bit boundary
  const size_t spare2_size = num_bits - 68 - text_size;
  text = bs.ToString(68, text_size);

  // TODO(schwehr): Is this correct?
  if (!spare2_size)
    spare2 = 0;
  else
    spare2 = bs.ToUnsignedInt(68, spare2_size);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

Ais8_1_11::Ais8_1_11(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), day(0), hour(0), minute(0), wind_ave(0),
      wind_gust(0), wind_dir(0), wind_gust_dir(0), air_temp(0.0), rel_humid(0),
      dew_point(0.0), air_pres(0.0), air_pres_trend(0), horz_vis(0.0),
      water_level(0.0), water_level_trend(0), surf_cur_speed(0.0),
      surf_cur_dir(0), cur_speed_2(0.0), cur_dir_2(0), cur_depth_2(0),
      cur_speed_3(0.0), cur_dir_3(0), cur_depth_3(0), wave_height(0.0),
      wave_period(0), wave_dir(0), swell_height(0.0), swell_period(0),
      swell_dir(0), sea_state(0), water_temp(0.0), precip_type(0),
      salinity(0.0), ice(0), spare2(0), extended_water_level(0) {
  assert(dac == 1);
  assert(fi == 11);

  if (num_chars != 59) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;  // 352 + 2 spares to be 6 bit aligned
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(56);
  // Reverse order lat/lng!
  // TODO(schwehr): Reverse order, or just reverse bit count? Compare 6.1.14!
  float y = bs.ToInt(56, 24) / 60000.;
  float x = bs.ToInt(80, 25) / 60000.;
  position = AisPoint(x, y);

  day = bs.ToUnsignedInt(105, 5);
  hour = bs.ToUnsignedInt(110, 5);
  minute = bs.ToUnsignedInt(115, 6);
  wind_ave = bs.ToUnsignedInt(121, 7);
  wind_gust = bs.ToUnsignedInt(128, 7);
  wind_dir = bs.ToUnsignedInt(135, 9);
  wind_gust_dir = bs.ToUnsignedInt(144, 9);
  air_temp = bs.ToUnsignedInt(153, 11) / 10. - 60;
  rel_humid = bs.ToUnsignedInt(164, 7);
  dew_point = bs.ToUnsignedInt(171, 10) / 10. - 20;  // TODO(schwehr): verify
  air_pres = bs.ToUnsignedInt(181, 9) + 800;
  air_pres_trend = bs.ToUnsignedInt(190, 2);
  horz_vis = bs.ToUnsignedInt(192, 8) / 10.;
  // TODO(schwehr): verify for -10.0 to 30.0
  water_level = bs.ToUnsignedInt(200, 9) / 10. - 10;
  water_level_trend = bs.ToUnsignedInt(209, 2);
  surf_cur_speed = bs.ToUnsignedInt(211, 8) / 10.;
  surf_cur_dir = bs.ToUnsignedInt(219, 9);
  cur_speed_2 = bs.ToUnsignedInt(228, 8) / 10.;
  cur_dir_2 = bs.ToUnsignedInt(236, 9);
  cur_depth_2 = bs.ToUnsignedInt(245, 5);
  cur_speed_3 = bs.ToUnsignedInt(250, 8) / 10.;
  cur_dir_3 = bs.ToUnsignedInt(258, 9);
  cur_depth_3 = bs.ToUnsignedInt(267, 5);

  wave_height = bs.ToUnsignedInt(272, 8) / 10.;
  wave_period = bs.ToUnsignedInt(280, 6);
  wave_dir = bs.ToUnsignedInt(286, 9);
  swell_height = bs.ToUnsignedInt(295, 8) / 10.;
  swell_period = bs.ToUnsignedInt(303, 6);
  swell_dir = bs.ToUnsignedInt(309, 9);

  sea_state = bs.ToUnsignedInt(318, 4);
  // TODO(schwehr): verify for -10.0 to +50.0
  water_temp = bs.ToUnsignedInt(322, 10) / 10. - 10;
  precip_type = bs.ToUnsignedInt(332, 3);
  salinity = bs.ToUnsignedInt(335, 9) / 10.0;  // Part per mil (1/1000).
  ice = bs.ToUnsignedInt(344, 2);
  // There is no way to know which meaning to attach to the following 6 bits
  // TODO(schwehr): how to treat this spare vrs water level?
  spare2 = bs.ToUnsignedInt(346, 6);
  bs.SeekRelative(-6);
  extended_water_level = bs.ToUnsignedInt(346, 6);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// No 8_1_12

// IMO Circ 289 - Fairway Closed
// See also Circ 236
Ais8_1_13::Ais8_1_13(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), reason(), location_from(), location_to(),
      radius(0), units(0), day_from(0), month_from(0), hour_from(0),
      minute_from(0), day_to(0), month_to(0), hour_to(0), minute_to(0),
      spare2(0) {
  assert(dac == 1);
  assert(fi == 13);

  if (num_bits != 472) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(56);
  reason = bs.ToString(56, 120);
  location_from = bs.ToString(176, 120);
  location_to = bs.ToString(296, 120);
  radius = bs.ToUnsignedInt(416, 10);
  units = bs.ToUnsignedInt(426, 2);
  day_from = bs.ToUnsignedInt(428, 5);
  month_from = bs.ToUnsignedInt(433, 4);
  hour_from = bs.ToUnsignedInt(437, 5);
  minute_from = bs.ToUnsignedInt(442, 6);
  day_to = bs.ToUnsignedInt(448, 5);
  month_to = bs.ToUnsignedInt(453, 4);
  hour_to = bs.ToUnsignedInt(457, 5);
  minute_to = bs.ToUnsignedInt(462, 6);
  spare2 = bs.ToUnsignedInt(468, 4);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// No 8_1_14

// IMO Circ 289 - Extended Shipdata - Air gap
// See also Circ 236
Ais8_1_15::Ais8_1_15(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), air_draught(0.0), spare2(0) {
  assert(dac == 1);
  assert(fi == 15);

  if (num_bits != 72) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(56);
  air_draught = bs.ToUnsignedInt(56, 11) / 10.;
  spare2 = bs.ToUnsignedInt(67, 5);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// IMO Circ 289 - Number of persons on board
// See also Circ 236
// TODO(schwehr): there might also be an addressed version?
Ais8_1_16::Ais8_1_16(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), persons(0), spare2(0) {
  assert(dac == 1);
  assert(fi == 16);

  if (num_bits != 72) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(56);
  persons = bs.ToUnsignedInt(56, 13);
  spare2 = bs.ToUnsignedInt(69, 3);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

Ais8_1_17_Target::Ais8_1_17_Target()
    : type(0), id(), spare(0), cog(0), timestamp(0), sog(0) {}

// IMO Circ 289 - VTS Generated/Synthetic Targets
// See also Circ 236
Ais8_1_17::Ais8_1_17(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad) {
  assert(dac == 1);
  assert(fi == 17);

  const size_t num_targets = (num_bits - 56) / 120;
  const size_t extra_bits = (num_bits - 56) % 120;

  if (extra_bits || num_targets > 4) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(56);
  for (size_t target_num = 0; target_num < num_targets; target_num++) {
    Ais8_1_17_Target target;
    const size_t start = 56 + (120 * target_num);
    target.type = bs.ToUnsignedInt(start, 2);
    target.id = bs.ToString(start + 2, 42);
    target.spare = bs.ToUnsignedInt(start + 44, 4);
    // booo - lat, lon inverse order
    float y = bs.ToInt(start + 48, 24) / 60000.;
    float x = bs.ToInt(start + 72, 25) / 60000.;
    target.position = AisPoint(x, y);

    target.cog = bs.ToUnsignedInt(start + 97, 9);
    target.timestamp = bs.ToUnsignedInt(start + 106, 6);
    target.sog = bs.ToUnsignedInt(start + 112, 8);
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// No msg 8_1_18

// IMO Circ 289 - Marine traffic signal
Ais8_1_19::Ais8_1_19(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), link_id(0), name(), status(0), signal(0),
      utc_hour_next(0), utc_min_next(0), spare2() {
  assert(dac == 1);
  assert(fi == 19);

  // Some people transmit without the idiodic spare padding
  if (num_bits != 258 && num_bits != 360) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(56);
  link_id = bs.ToUnsignedInt(56, 10);
  name = bs.ToString(66, 120);
  position = bs.ToAisPoint(186, 49);
  status = bs.ToUnsignedInt(235, 2);
  signal = bs.ToUnsignedInt(237, 5);
  utc_hour_next = bs.ToUnsignedInt(242, 5);
  utc_min_next = bs.ToUnsignedInt(247, 6);
  next_signal = bs.ToUnsignedInt(253, 5);
  if (num_bits == 360) {
    spare2[0] = bs.ToUnsignedInt(258, 32);
    spare2[1] = bs.ToUnsignedInt(290, 32);
    spare2[2] = bs.ToUnsignedInt(322, 32);
    spare2[3] = bs.ToUnsignedInt(354, 6);
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// No 8_1_20

// IMO Circ 289 - Weather observation report from ship
// See also Circ 236
Ais8_1_21::Ais8_1_21(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), type_wx_report(0), location(), utc_day(0),
      utc_hour(), utc_min(), horz_viz(0.0), humidity(0), wind_speed(0),
      wind_dir(0), pressure(0.0), pressure_tendency(0), air_temp(0.0),
      water_temp(0.0), wave_period(0), wave_height(0.0), wave_dir(0),
      swell_height(0.0), swell_dir(0), swell_period(0), spare2(0),
      utc_month(0), cog(0), sog(0.0), heading(0), rel_pressure(0.0),
      wind_speed_ms(0.0), wind_dir_rel(0), wind_speed_rel(0.0),
      wind_gust_speed(0.0), wind_gust_dir(0), air_temp_raw(0),
      water_temp_raw(0), wx(), cloud_total(0), cloud_low(0),
      cloud_low_type(0), cloud_middle_type(0), cloud_high_type(0),
      alt_lowest_cloud_base(0.0), swell_dir_2(0), swell_period_2(0),
      swell_height_2(0.0), ice_thickness(0.0), ice_accretion(0),
      ice_accretion_cause(0), sea_ice_concentration(0), amt_type_ice(0),
      ice_situation(0), ice_devel(0), bearing_ice_edge(0) {
  assert(dac == 1);
  assert(fi == 21);

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

  bs.SeekTo(56);
  type_wx_report = bs[56];
  if (!type_wx_report) {
    // WX obs from ship
    location = bs.ToString(57, 120);
    position = bs.ToAisPoint(177, 49);
    utc_day = bs.ToUnsignedInt(226, 5);
    utc_hour = bs.ToUnsignedInt(231, 5);
    utc_min = bs.ToUnsignedInt(236, 6);
    wx[0] = bs.ToUnsignedInt(242, 4);  // TODO(schwehr): set wx[1] and wx[2]?
    horz_viz = bs.ToUnsignedInt(246, 8) / 10.;  // nautical miles
    humidity = bs.ToUnsignedInt(254, 7);  // %
    wind_speed = bs.ToUnsignedInt(261, 7);  // ave knots
    wind_dir = bs.ToUnsignedInt(268, 9);
    pressure = bs.ToUnsignedInt(277, 9);  // hPa
    pressure_tendency = bs.ToUnsignedInt(286, 4);
    // TODO(schwehr): is air_temp correct?
    air_temp = bs.ToInt(290, 11) / 10.;  // C
    water_temp = bs.ToUnsignedInt(301, 10) / 10. - 10;  // C
    wave_period = bs.ToUnsignedInt(311, 6);  // s
    wave_height = bs.ToUnsignedInt(317, 8) / 10.;
    wave_dir = bs.ToUnsignedInt(325, 9);
    swell_height = bs.ToUnsignedInt(334, 8) / 10.;  // m
    swell_dir = bs.ToUnsignedInt(342, 9);
    swell_period = bs.ToUnsignedInt(351, 6);  // s
    spare2 = bs.ToUnsignedInt(357, 3);
  } else {
    // Type 1: WMO OBS from ship.
    float x = (bs.ToUnsignedInt(57, 16) / 100.) - 180;
    float y = (bs.ToUnsignedInt(73, 15) / 100.) - 180;
    position = AisPoint(x, y);

    utc_month = bs.ToUnsignedInt(88, 4);
    utc_day = bs.ToUnsignedInt(92, 6);
    utc_hour = bs.ToUnsignedInt(98, 5);
    utc_min = bs.ToUnsignedInt(102, 3) * 10;
    cog = bs.ToUnsignedInt(106, 7) * 5;
    sog = bs.ToUnsignedInt(113, 5) * 0.5;
    heading = bs.ToUnsignedInt(118, 7) *5;  // Assume this is true degrees????
    pressure = bs.ToUnsignedInt(125, 11) / 10. + 900;
    rel_pressure = bs.ToUnsignedInt(136, 10) / 10. -50;
    pressure_tendency = bs.ToUnsignedInt(146, 4);
    wind_dir = bs.ToUnsignedInt(150, 7) * 5;
    wind_speed_ms = bs.ToUnsignedInt(157, 8) * 0.5;  // m/s
    wind_dir_rel = bs.ToUnsignedInt(165, 7) * 5;
    wind_speed_rel = bs.ToUnsignedInt(172, 8) * 0.5;  // m/s
    wind_gust_speed = bs.ToUnsignedInt(180, 8) * 0.5;  // m/s
    wind_gust_dir = bs.ToUnsignedInt(188, 7) * 5;
    // 0C = 273.15 Kelvin
    // TODO(schwehr): change this to celcius
    air_temp_raw = bs.ToUnsignedInt(195, 10);
    humidity = bs.ToUnsignedInt(205, 7);
    water_temp_raw = bs.ToUnsignedInt(212, 9);  // TODO(schwehr): Change to C.

    auto pow2 = [](unsigned int val) { return val * val; };
    horz_viz = pow2(bs.ToUnsignedInt(221, 6)) * 13.073;  // m
    wx[0] = bs.ToUnsignedInt(227, 9);  // current
    wx[1] = bs.ToUnsignedInt(236, 5);  // past 1
    wx[2] = bs.ToUnsignedInt(241, 5);  // past 2
    cloud_total = bs.ToUnsignedInt(246, 4) * 10;
    cloud_low = bs.ToUnsignedInt(250, 4);
    cloud_low_type = bs.ToUnsignedInt(254, 6);
    cloud_middle_type = bs.ToUnsignedInt(260, 6);
    cloud_high_type = bs.ToUnsignedInt(266, 6);
    alt_lowest_cloud_base = pow2(bs.ToUnsignedInt(272, 7)) * 0.16;
    wave_period = bs.ToUnsignedInt(279, 5);  // s
    wave_height = bs.ToUnsignedInt(284, 6) * 0.5;  // m
    swell_dir = bs.ToUnsignedInt(290, 6) * 10;
    swell_period = bs.ToUnsignedInt(296, 5);  // s
    swell_height = bs.ToUnsignedInt(301, 6) * 0.5;  // m
    swell_dir_2 = bs.ToUnsignedInt(307, 6) * 10;
    swell_period_2 = bs.ToUnsignedInt(313, 5);  // s
    swell_height_2 = bs.ToUnsignedInt(318, 6) * 0.5;  // m
    ice_thickness = bs.ToUnsignedInt(324, 7) / 100.;  // m.  Network is cm,
    ice_accretion = bs.ToUnsignedInt(331, 3);
    ice_accretion_cause = bs.ToUnsignedInt(334, 3);
    sea_ice_concentration = bs.ToUnsignedInt(337, 5);
    amt_type_ice = bs.ToUnsignedInt(342, 4);
    ice_situation = bs.ToUnsignedInt(346, 5);
    ice_devel = bs.ToUnsignedInt(351, 5);
    bearing_ice_edge = bs.ToUnsignedInt(356, 4) * 45;
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// See ais8_1_22.cpp
// No 8_1_23

// IMO Circ 289 - Extended ship static and voyage-related
// See also Circ 236
Ais8_1_24::Ais8_1_24(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), link_id(0), air_draught(0.0), last_port(),
      next_ports(), solas_status(), ice_class(0), shaft_power(0), vhf(0),
      lloyds_ship_type(), gross_tonnage(0), laden_ballast(0), heavy_oil(0),
      light_oil(0), diesel(0), bunker_oil(0), persons(0), spare2(0) {
  assert(dac == 1);
  assert(fi == 24);

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

  bs.SeekTo(56);
  link_id = bs.ToUnsignedInt(56, 10);
  air_draught = bs.ToUnsignedInt(66, 13) / 10.;  // m
  last_port = bs.ToString(79, 30);
  next_ports[0] = bs.ToString(109, 30);
  next_ports[1] = bs.ToString(139, 30);

  // TODO(schwehr): enum list of param types
  // 0 NA, 1 operational, 2 SNAFU, 3 no data
  for (size_t equip_num = 0; equip_num < 26; equip_num++) {
    solas_status[equip_num] = bs.ToUnsignedInt(169 + 2 * equip_num, 2);
  }
  ice_class = bs.ToUnsignedInt(221, 4);
  shaft_power = bs.ToUnsignedInt(225, 18);  // horses
  vhf = bs.ToUnsignedInt(243, 12);
  lloyds_ship_type = bs.ToString(255, 42);
  gross_tonnage = bs.ToUnsignedInt(297, 18);
  laden_ballast = bs.ToUnsignedInt(315, 2);
  heavy_oil = bs.ToUnsignedInt(317, 2);
  light_oil = bs.ToUnsignedInt(319, 2);
  diesel = bs.ToUnsignedInt(321, 2);
  bunker_oil = bs.ToUnsignedInt(323, 14);  // tonnes
  persons = bs.ToUnsignedInt(337, 13);
  spare2 = bs.ToUnsignedInt(350, 10);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}


// No 8_1_25
// See ais8_1_26.cpp

// IMO Circ 289 - Route information
// See also Circ 236
Ais8_1_27::Ais8_1_27(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), link_id(0), sender_type(0), route_type(0),
      utc_month(0), utc_day(0), utc_hour(0), utc_min(0), duration(0) {
  assert(dac == 1);
  assert(fi == 27);

  const size_t num_waypoints = (num_bits - 117) / 55;
  const size_t extra_bits = (num_bits - 117) % 55;

  if (extra_bits || num_waypoints > 16) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(56);
  link_id = bs.ToUnsignedInt(56, 10);
  sender_type = bs.ToUnsignedInt(66, 3);
  route_type = bs.ToUnsignedInt(69, 5);
  utc_month = bs.ToUnsignedInt(74, 4);
  utc_day = bs.ToUnsignedInt(78, 5);
  utc_hour = bs.ToUnsignedInt(83, 5);
  utc_min = bs.ToUnsignedInt(88, 6);
  duration = bs.ToUnsignedInt(94, 18);
  // TODO(schwehr): manage the case where num_waypoints does not match
  // const size_t num_waypoints_stated = bs.ToUnsignedInt(112, 5);
  for (size_t waypoint_num = 0; waypoint_num < num_waypoints; waypoint_num++) {
    const size_t start = 117 + 55*waypoint_num;
    waypoints.push_back(bs.ToAisPoint(start, 55));
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// No 8_1_28

// IMO Circ 289 - Text description
// See also Circ 236
Ais8_1_29::Ais8_1_29(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), link_id(0), spare2(0) {
  assert(dac == 1);
  assert(fi == 29);

  if (num_bits < 72 || num_bits > 1032) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  bs.SeekTo(56);
  link_id = bs.ToUnsignedInt(56, 10);
  size_t text_bits = ((num_bits - 66) / 6) * 6;
  text = bs.ToString(66, text_bits);
  const size_t spare2_bits = num_bits - 66 - text_bits;
  if (spare2_bits) {
    const size_t start = 66 + text_bits;
    spare2 = bs.ToUnsignedInt(start, spare2_bits);
  } else {
    spare2 = 0;
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}



// IMO Circ 289 - Meteorological and Hydrographic data
// See also Circ 236
Ais8_1_31::Ais8_1_31(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), position_accuracy(0), utc_day(0), utc_hour(0),
      utc_min(0), wind_ave(0), wind_gust(0), wind_dir(0), wind_gust_dir(0),
      air_temp(0.0), rel_humid(0), dew_point(0.0), air_pres(0.0),
      air_pres_trend(0), horz_vis(0.0), water_level(0.0),
      water_level_trend(0), surf_cur_speed(0.0), surf_cur_dir(0),
      cur_speed_2(0.0), cur_dir_2(0), cur_depth_2(0), cur_speed_3(0.0),
      cur_dir_3(0), cur_depth_3(0), wave_height(0.0), wave_period(0),
      wave_dir(0), swell_height(0.0), swell_period(0), swell_dir(0),
      sea_state(0), water_temp(0.0), precip_type(0), salinity(0.0),
      ice(0), spare2(0) {
  assert(dac == 1);
  assert(fi == 31);

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

  bs.SeekTo(56);
  position = bs.ToAisPoint(56, 49);
  position_accuracy = bs[105];
  utc_day = bs.ToUnsignedInt(106, 5);
  utc_hour = bs.ToUnsignedInt(111, 5);
  utc_min = bs.ToUnsignedInt(116, 6);

  wind_ave = bs.ToUnsignedInt(122, 7);  // kts
  wind_gust = bs.ToUnsignedInt(129, 7);  // kts
  wind_dir = bs.ToUnsignedInt(136, 9);
  wind_gust_dir = bs.ToUnsignedInt(145, 9);
  air_temp = bs.ToInt(154, 11) / 10.;  // C
  rel_humid = bs.ToUnsignedInt(165, 7);
  dew_point = bs.ToInt(172, 10)/ 10.;  // TODO(schwehr): How is this mapped?
  air_pres = (bs.ToUnsignedInt(182, 9) + 800) / 100.0;  // Pa
  air_pres_trend = bs.ToUnsignedInt(191, 2);

  horz_vis = bs.ToUnsignedInt(193, 8) / 10.;  // NM
  water_level = bs.ToUnsignedInt(201, 12) / 100. - 10;  // m
  water_level_trend = bs.ToUnsignedInt(213, 2);

  surf_cur_speed = bs.ToUnsignedInt(215, 8) / 10.;
  surf_cur_dir = bs.ToUnsignedInt(223, 9);
  cur_speed_2 = bs.ToUnsignedInt(232, 8) / 10.;  // kts
  cur_dir_2 = bs.ToUnsignedInt(240, 9);
  cur_depth_2 = bs.ToUnsignedInt(249, 5);  // m
  cur_speed_3 = bs.ToUnsignedInt(254, 8) / 10.;  // kts
  cur_dir_3 = bs.ToUnsignedInt(262, 9);
  cur_depth_3 = bs.ToUnsignedInt(271, 5);  // m

  wave_height = bs.ToUnsignedInt(276, 8);  // m

  wave_period = bs.ToUnsignedInt(284, 6);
  wave_dir = bs.ToUnsignedInt(290, 9);
  swell_height = bs.ToUnsignedInt(299, 8) / 10.;
  swell_period = bs.ToUnsignedInt(307, 6);
  swell_dir = bs.ToUnsignedInt(313, 9);
  sea_state = bs.ToUnsignedInt(322, 4);  // beaufort scale - Table 1.2
  water_temp = bs.ToInt(326, 10) / 10.;
  precip_type = bs.ToUnsignedInt(336, 3);
  salinity = bs.ToUnsignedInt(339, 9) / 10.;
  ice = bs.ToUnsignedInt(348, 2);  // yes/no/undef/unknown
  spare2 = bs.ToUnsignedInt(350, 10);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}


}  // namespace libais
