// The USCG RTCM SC 121 Environmental Sensor Message tested in
// Tampa, FL from 200[78] util ?

// http://www.dtic.mil/cgi-bin/GetTRDoc?AD=ADA504755
// Phase I Summary Report on AIS Transmit Project (Environmental Message)

#include "ais.h"

namespace libais {

Ais8_1_26_Location::Ais8_1_26_Location(const AisBitset &bs,
                                       const size_t offset) {
  position = bs.ToAisPoint(offset, 55);
  z = bs.ToUnsignedInt(offset + 55, 11) / 10.;
  owner = bs.ToUnsignedInt(offset + 66, 4);
  timeout = bs.ToUnsignedInt(offset + 70, 3);
  spare = bs.ToUnsignedInt(offset + 73, 12);
}

Ais8_1_26_Station::Ais8_1_26_Station(const AisBitset &bs,
                                     const size_t offset) {
  name = bs.ToString(offset, 84);
  spare = bs.ToUnsignedInt(offset + 84, 1);
}

Ais8_1_26_Wind::Ais8_1_26_Wind(const AisBitset &bs,
                               const size_t offset) {
  wind_speed = bs.ToUnsignedInt(offset, 7);
  wind_gust  = bs.ToUnsignedInt(offset + 7, 7);  // knots
  wind_dir = bs.ToUnsignedInt(offset + 14, 9);
  wind_gust_dir = bs.ToUnsignedInt(offset + 23, 9);
  sensor_type = bs.ToUnsignedInt(offset + 32, 3);
  wind_forecast = bs.ToUnsignedInt(offset + 35, 7);
  wind_gust_forecast = bs.ToUnsignedInt(offset + 42, 7);  // knots
  wind_dir_forecast = bs.ToUnsignedInt(offset + 49, 9);
  utc_day_forecast = bs.ToUnsignedInt(offset + 58, 5);
  utc_hour_forecast = bs.ToUnsignedInt(offset + 63, 5);
  utc_min_forecast = bs.ToUnsignedInt(offset + 68, 6);
  duration = bs.ToUnsignedInt(offset + 74, 8);
  spare = bs.ToUnsignedInt(offset + 82, 3);
}

Ais8_1_26_WaterLevel::Ais8_1_26_WaterLevel(const AisBitset &bs,
                                           const size_t offset) {
  type = bs[offset];
  level = bs.ToInt(offset + 1, 16) / 100.;
  trend = bs.ToUnsignedInt(offset + 17, 2);
  vdatum = bs.ToUnsignedInt(offset + 19, 5);
  sensor_type = bs.ToUnsignedInt(offset + 24, 3);
  forecast_type = bs[offset + 27];
  level_forecast = bs.ToInt(offset + 28, 16) / 100.;
  utc_day_forecast = bs.ToUnsignedInt(offset + 44, 5);
  utc_hour_forecast = bs.ToUnsignedInt(offset + 49, 5);
  utc_min_forecast = bs.ToUnsignedInt(offset + 54, 6);
  duration = bs.ToUnsignedInt(offset + 60, 8);
  spare = bs.ToUnsignedInt(offset + 68, 17);
}

Ais8_1_26_Curr2D::Ais8_1_26_Curr2D(const AisBitset &bs,
                                   const size_t offset) {
  for (size_t idx = 0; idx < 3; idx++) {
    size_t start = offset + idx * 26;
    currents[idx].speed = bs.ToUnsignedInt(start, 8) / 10.;
    currents[idx].dir = bs.ToUnsignedInt(start + 8, 9);
    currents[idx].depth = bs.ToUnsignedInt(start + 17, 9);
  }
  type = bs.ToUnsignedInt(offset + 78, 3);
  spare = bs.ToUnsignedInt(offset + 81, 4);
}

Ais8_1_26_Curr3D::Ais8_1_26_Curr3D(const AisBitset &bs,
                                   const size_t offset) {
  for (size_t idx = 0; idx < 2; idx++) {
    size_t start = offset + idx * 33;
    currents[idx].north = bs.ToUnsignedInt(start, 8) / 10.;
    currents[idx].east = bs.ToUnsignedInt(start + 8, 8) / 10.;
    currents[idx].up = bs.ToUnsignedInt(start + 16, 8) / 10.;
    currents[idx].depth = bs.ToUnsignedInt(start + 24, 9);
  }
  type = bs.ToUnsignedInt(offset + 66, 3);
  spare = bs.ToUnsignedInt(offset + 69, 16);
}

Ais8_1_26_HorzFlow::Ais8_1_26_HorzFlow(const AisBitset &bs,
                                       const size_t offset) {
  for (size_t idx = 0; idx < 2; idx++) {
    size_t start = offset + idx * 42;
    currents[idx].bearing = bs.ToUnsignedInt(start, 9);
    currents[idx].dist = bs.ToUnsignedInt(start + 9, 7);
    currents[idx].speed = bs.ToUnsignedInt(start + 16, 8) / 10.;
    currents[idx].dir = bs.ToUnsignedInt(start + 24, 9);
    currents[idx].level = bs.ToUnsignedInt(start + 33, 9);
  }
  spare = bs[offset + 84];
}

Ais8_1_26_SeaState::Ais8_1_26_SeaState(const AisBitset &bs,
                                       const size_t offset) {
  swell_height = bs.ToUnsignedInt(offset, 8) / 10.;
  swell_period = bs.ToUnsignedInt(offset + 8, 6);
  swell_dir = bs.ToUnsignedInt(offset + 14, 9);
  sea_state = bs.ToUnsignedInt(offset + 23, 4);
  swell_sensor_type = bs.ToUnsignedInt(offset + 27, 3);
  water_temp = bs.ToInt(offset + 30, 10) / 10.;
  water_temp_depth = bs.ToUnsignedInt(offset + 40, 7) / 10.;
  water_sensor_type = bs.ToUnsignedInt(offset + 47, 3);
  wave_height = bs.ToUnsignedInt(offset + 50, 8) / 10.;
  wave_period = bs.ToUnsignedInt(offset + 58, 6);
  wave_dir = bs.ToUnsignedInt(offset + 64, 9);
  wave_sensor_type = bs.ToUnsignedInt(offset + 73, 3);
  salinity = bs.ToUnsignedInt(offset + 76, 9) / 10.;
}

Ais8_1_26_Salinity::Ais8_1_26_Salinity(const AisBitset &bs,
                                       const size_t offset) {
  water_temp = bs.ToUnsignedInt(offset, 10) / 10. - 10;
  conductivity = bs.ToUnsignedInt(offset + 10, 10) / 100.;
  pressure = bs.ToUnsignedInt(offset + 20, 16) / 10.;
  salinity = bs.ToUnsignedInt(offset + 36, 9) / 10.;
  salinity_type = bs.ToUnsignedInt(offset + 45, 2);
  sensor_type = bs.ToUnsignedInt(offset + 47, 3);
  spare[0] = bs.ToUnsignedInt(offset + 50, 32);
  spare[1] = bs.ToUnsignedInt(offset + 82, 3);
}

Ais8_1_26_Wx::Ais8_1_26_Wx(const AisBitset &bs,
                           const size_t offset) {
  air_temp = bs.ToInt(offset, 11) / 10.;
  air_temp_sensor_type = bs.ToUnsignedInt(offset + 11, 3);
  precip = bs.ToUnsignedInt(offset + 14, 2);
  horz_vis = bs.ToUnsignedInt(offset + 16, 8) / 10.;
  dew_point = bs.ToInt(offset + 24, 10) / 10.;
  dew_point_type = bs.ToUnsignedInt(offset + 34, 3);
  air_pressure = (bs.ToUnsignedInt(offset + 37, 9) + 800) / 100.0;  // Pa.
  air_pressure_trend = bs.ToUnsignedInt(offset + 46, 2);
  air_pressor_type = bs.ToUnsignedInt(offset + 48, 3);
  salinity = bs.ToUnsignedInt(offset + 51, 9) / 10.;
  spare = bs.ToUnsignedInt(offset + 60, 25);
}

Ais8_1_26_AirDraught::Ais8_1_26_AirDraught(const AisBitset &bs,
                                           const size_t offset) {
  draught = bs.ToUnsignedInt(offset, 13) / 100.;
  gap = bs.ToUnsignedInt(offset + 13, 13) / 10.;
  trend = bs.ToUnsignedInt(offset + 26, 2);
  forecast_gap = bs.ToUnsignedInt(offset + 28, 13) / 10.;
  utc_day_forecast = bs.ToUnsignedInt(offset + 41, 5);
  utc_hour_forecast = bs.ToUnsignedInt(offset + 46, 5);
  utc_min_forecast = bs.ToUnsignedInt(offset + 51, 6);
  spare = bs.ToUnsignedInt(offset + 57, 28);
}

// TODO(schwehr): Refactor to be like the 8:367:22 factory.
Ais8_1_26_SensorReport*
ais8_1_26_sensor_report_factory(const AisBitset &bs,
                                const size_t offset) {
  const Ais8_1_26_SensorEnum rpt_type =
      (Ais8_1_26_SensorEnum)bs.ToUnsignedInt(offset, 4);

  // WARNING: out of order decoding
  // Only get the report header if we can decode the type
  const size_t rpt_start = offset + 27;  // skip tp after site_id
  bs.SeekTo(rpt_start);
  Ais8_1_26_SensorReport *rpt = nullptr;
  switch (rpt_type) {
  case AIS8_1_26_SENSOR_LOCATION:
    rpt = new Ais8_1_26_Location(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_STATION:
    rpt = new Ais8_1_26_Station(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_WIND:
    rpt = new Ais8_1_26_Wind(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_WATER_LEVEL:
    rpt = new Ais8_1_26_WaterLevel(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_CURR_2D:
    rpt = new Ais8_1_26_Curr2D(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_CURR_3D:
    rpt = new Ais8_1_26_Curr3D(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_HORZ_FLOW:
    rpt = new Ais8_1_26_HorzFlow(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_SEA_STATE:
    rpt = new Ais8_1_26_SeaState(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_SALINITY:
    rpt = new Ais8_1_26_Salinity(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_WX:
    rpt = new Ais8_1_26_Wx(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_AIR_DRAUGHT:
    rpt = new Ais8_1_26_AirDraught(bs, rpt_start);
    break;
  case AIS8_1_26_SENSOR_RESERVED_11: break;  // Leave rpt == 0 to indicate error
  case AIS8_1_26_SENSOR_RESERVED_12: break;
  case AIS8_1_26_SENSOR_RESERVED_13: break;
  case AIS8_1_26_SENSOR_RESERVED_14: break;
  case AIS8_1_26_SENSOR_RESERVED_15: break;
  default:
    {}  // Leave rpt == 0 to indicate error
  }

  if (!rpt)
    return rpt;

  rpt->report_type = rpt_type;
  bs.SeekTo(offset + 4);
  rpt->utc_day = bs.ToUnsignedInt(offset + 4, 5);
  rpt->utc_hr = bs.ToUnsignedInt(offset + 9, 5);
  rpt->utc_min = bs.ToUnsignedInt(offset + 14, 6);
  rpt->site_id = bs.ToUnsignedInt(offset + 20, 7);
  return rpt;
}

Ais8_1_26::Ais8_1_26(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad) {
  assert(dac == 1);
  assert(fi == 26);

  if (168 > num_bits || num_bits > 1098) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  const size_t num_sensor_reports = (num_bits - 56) / AIS8_1_26_REPORT_SIZE;

  // TODO(schwehr): what to do about extra data in sensor report msg 8_1_26?
  // if ((num_bits - 56) % AIS8_1_26_REPORT_SIZE)

  for (size_t report_idx = 0; report_idx < num_sensor_reports; report_idx++) {
    const size_t start = 56 + report_idx * AIS8_1_26_REPORT_SIZE;
    bs.SeekTo(start);
    Ais8_1_26_SensorReport *sensor = ais8_1_26_sensor_report_factory(bs, start);
    if (sensor) {
        reports.push_back(sensor);
    } else {
      status = AIS_ERR_BAD_SUB_SUB_MSG;
      return;
    }
  }

  // TODO(schwehr): Enable this assert after fixing the message.
  // assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// TODO(schwehr): Use unique_ptr to manage memory.
Ais8_1_26::~Ais8_1_26() {
  for (size_t i = 0; i < reports.size(); i++) {
    delete reports[i];
    reports[i] = nullptr;
  }
}

}  // namespace libais
