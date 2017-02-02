// River Information Systems ECE-TRANS-SC3-2006-10r-RIS.pdf
// http://www.unece.org/fileadmin/DAM/trans/doc/finaldocs/sc3/ECE-TRANS-SC3-176e.pdf

#include "ais.h"

namespace libais {

// Inland ship static and voyage related data
Ais8_200_10::Ais8_200_10(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), length(0.0), beam(0.0), ship_type(0),
      haz_cargo(0), draught(0.0), loaded(0), speed_qual(0), course_qual(0),
      heading_qual(0), spare2(0) {
  assert(dac == 200);
  assert(fi == 10);

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

  bs.SeekTo(56);
  eu_id = bs.ToString(56, 48);
  length = bs.ToUnsignedInt(104, 13) / 10.;  // m
  beam = bs.ToUnsignedInt(117, 10) / 10.;  // m
  ship_type = bs.ToUnsignedInt(127, 14);
  haz_cargo = bs.ToUnsignedInt(141, 3);
  draught = bs.ToUnsignedInt(144, 11) / 10.;  // m
  loaded = bs.ToUnsignedInt(155, 2);
  speed_qual = bs[157];
  course_qual = bs[158];
  heading_qual = bs[159];
  spare2 = bs.ToUnsignedInt(160, 8);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// River Information Systems ECE-TRANS-SC3-2006-10r-RIS.pdf
Ais8_200_23::Ais8_200_23(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), utc_year_start(0), utc_month_start(0),
      utc_day_start(0), utc_year_end(0), utc_month_end(0), utc_day_end(0),
      utc_hour_start(0), utc_min_start(0), utc_hour_end(0), utc_min_end(0),
      type(0), min(0), max(0), classification(0), wind_dir(0), spare2(0) {
  assert(dac == 200);
  assert(fi == 23);

  if (num_bits != 256) {
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

  // TODO(schwehr): Figure out the correct bits & test against actual messages.

  // The total for the bit column in table 2.11 of ECE/TRANS/SC.3/2006/10
  // add up to 256.  However, start date and end date fields are lists with
  // totals of 17 bits each, but within the details, the spec refers to 18 bits
  // for each.  It is likely the year should be one less bit.
  utc_year_start = bs.ToUnsignedInt(56, 8);
  utc_month_start = bs.ToUnsignedInt(65, 4);
  utc_day_start = bs.ToUnsignedInt(69, 5);

  utc_year_end = bs.ToUnsignedInt(73, 8);
  utc_month_end = bs.ToUnsignedInt(82, 4);
  utc_day_end = bs.ToUnsignedInt(86, 5);

  utc_hour_start = bs.ToUnsignedInt(90, 5);
  utc_min_start = bs.ToUnsignedInt(95, 6);
  utc_hour_end = bs.ToUnsignedInt(101, 5);
  utc_min_end = bs.ToUnsignedInt(106, 6);

  position1 = bs.ToAisPoint(112, 55);
  position2 = bs.ToAisPoint(167, 55);

  type = bs.ToUnsignedInt(222, 4);
  min = bs.ToUnsignedInt(226, 9);
  max = bs.ToUnsignedInt(235, 9);
  classification = bs.ToUnsignedInt(244, 2);
  wind_dir = bs.ToUnsignedInt(246, 4);
  spare2 = bs.ToUnsignedInt(250, 6);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}


// River Information Systems ECE-TRANS-SC3-2006-10r-RIS.pdf
Ais8_200_24::Ais8_200_24(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad) {
  assert(dac == 200);
  assert(fi == 24);

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

  bs.SeekTo(56);
  bs.ToString(56, 12);
  for (size_t i = 0; i < 4; i++) {
    size_t start = 68 + 25*i;
    gauge_ids[i] = bs.ToUnsignedInt(start, 11);
    const int sign = bs[start + 11] ? 1 : -1;  // 0 negative, 1 pos
    // ERROR: the spec has a bit listing mistake
    levels[i] = sign * bs.ToUnsignedInt(start + 12, 13);
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

// River Information Systems ECE-TRANS-SC3-2006-10r-RIS.pdf
Ais8_200_40::Ais8_200_40(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), form(0), dir(0), stream_dir(0), status_raw(0),
      spare2(0) {
  assert(dac == 200);
  assert(fi == 40);

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

  bs.SeekTo(56);
  position = bs.ToAisPoint(56, 55);
  form = bs.ToUnsignedInt(111, 4);
  dir = bs.ToUnsignedInt(115, 9);  // degrees
  stream_dir = bs.ToUnsignedInt(124, 3);
  status_raw = bs.ToUnsignedInt(127, 30);
  // TODO(schwehr): Convert the status_raw to the 9 signal lights.
  // Appears to be a base 10 setup where each of 9 digits range from 0 to 7
  // for the light level.
  spare2 = bs.ToUnsignedInt(157, 11);

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}


// River Information Systems ECE-TRANS-SC3-2006-10r-RIS.pdf
//
// TODO(schwehr): Search the logs for the various possible sizes and make
//   tests based on those messages.
Ais8_200_55::Ais8_200_55(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), crew(0), passengers(0), yet_more_personnel(0) {
  assert(dac == 200);
  assert(fi == 55);

  // The specification says that there are 51 spare bits, but it is possible
  // that some transmitters may leave off the spare bits.
  if (num_bits != 88 && num_bits != 136 && num_bits != 168) {
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
  crew = bs.ToUnsignedInt(56, 8);
  passengers = bs.ToUnsignedInt(64, 13);
  yet_more_personnel = bs.ToUnsignedInt(77, 8);

  if (num_bits == 88) {
    spare2[0] = bs.ToUnsignedInt(85, 3);
  } else if (num_bits == 136) {
    spare2[0] = bs.ToUnsignedInt(85, 32);
    spare2[1] = bs.ToUnsignedInt(117, 19);
  } else {
    spare2[0] = bs.ToUnsignedInt(85, 32);
    spare2[1] = bs.ToUnsignedInt(117, 32);
    spare2[2] = bs.ToUnsignedInt(149, 19);
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

}  // namespace libais
