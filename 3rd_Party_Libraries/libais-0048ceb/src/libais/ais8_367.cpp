// 8:367:22 Defined by an email from Greg Johnson representing the
// USCG, Fall 2012.  Breaks from the RTCM and IMO Circular 289.
// "Area Notice Message Release Version: 1" 13 Aug 2012
//
// http://www.e-navigation.nl/content/geographic-notice
// http://www.e-navigation.nl/sites/default/files/asm_files/GN%20Release%20Version%201b.pdf

#include <cmath>

#include "ais.h"

namespace libais {

const size_t SUB_AREA_BITS = 96;

static int scale_multipliers[4] = {1, 10, 100, 1000};

Ais8_367_22_Circle::Ais8_367_22_Circle(const AisBitset &bs, const size_t offset)
    : precision(0), radius_m(0), spare(0) {
  const int scale_factor = bs.ToUnsignedInt(offset, 2);
  position = bs.ToAisPoint(offset + 2, 55);
  precision = bs.ToUnsignedInt(offset + 57, 3);
  radius_m =
      bs.ToUnsignedInt(offset + 60, 12) * scale_multipliers[scale_factor];
  spare = bs.ToUnsignedInt(offset + 72, 21);
}

Ais8_367_22_Rect::Ais8_367_22_Rect(const AisBitset &bs, const size_t offset)
    : precision(0), e_dim_m(0), n_dim_m(0), orient_deg(0), spare(0) {
  const int scale_factor = bs.ToUnsignedInt(offset, 2);
  position = bs.ToAisPoint(offset + 2, 55);
  precision = bs.ToUnsignedInt(offset + 57, 3);
  e_dim_m = bs.ToUnsignedInt(offset + 60, 8) * scale_multipliers[scale_factor];
  n_dim_m = bs.ToUnsignedInt(offset + 68, 8) * scale_multipliers[scale_factor];
  orient_deg = bs.ToUnsignedInt(offset + 76, 9);
  spare = bs.ToUnsignedInt(offset + 85, 8);
}

Ais8_367_22_Sector::Ais8_367_22_Sector(const AisBitset &bs, const size_t offset)
    : precision(0),
      radius_m(0),
      left_bound_deg(0),
      right_bound_deg(0),
      spare(0) {
  const int scale_factor = bs.ToUnsignedInt(offset, 2);
  position = bs.ToAisPoint(offset + 2, 55);
  precision = bs.ToUnsignedInt(offset + 57, 3);
  radius_m =
      bs.ToUnsignedInt(offset + 60, 12) * scale_multipliers[scale_factor];
  left_bound_deg = bs.ToUnsignedInt(offset + 72, 9);
  right_bound_deg = bs.ToUnsignedInt(offset + 81, 9);
  spare = bs.ToUnsignedInt(offset + 90, 3);
}

// Polyline or polygon.
Ais8_367_22_Poly::Ais8_367_22_Poly(const AisBitset &bs, const size_t offset,
                                   Ais8_366_22_AreaShapeEnum area_shape)
    : shape(area_shape), precision(0), spare(0) {
  const int scale_factor = bs.ToUnsignedInt(offset, 2);
  size_t poly_offset = offset + 2;
  for (size_t i = 0; i < 4; i++) {
    const int angle = bs.ToUnsignedInt(poly_offset, 10);
    poly_offset += 10;
    const int dist = bs.ToUnsignedInt(poly_offset, 11) *
                     scale_multipliers[scale_factor];
    poly_offset += 11;
    if (dist == 0) {
      break;
    }
    angles.push_back(angle);
    dists_m.push_back(dist);
  }
  spare = bs.ToUnsignedInt(offset + 86, 7);
}

Ais8_367_22_Text::Ais8_367_22_Text(const AisBitset &bs, const size_t offset) {
  text = string(bs.ToString(offset, 90));
  spare = bs.ToUnsignedInt(offset + 90, 3);
}

Ais8_367_22_SubArea *ais8_367_22_subarea_factory(const AisBitset &bs,
                                                 const size_t offset) {
  const Ais8_366_22_AreaShapeEnum area_shape =
      static_cast<Ais8_366_22_AreaShapeEnum>(bs.ToUnsignedInt(offset, 3));

  switch (area_shape) {
    case AIS8_366_22_SHAPE_CIRCLE:
      return new Ais8_367_22_Circle(bs, offset + 3);
    case AIS8_366_22_SHAPE_RECT:
      return new Ais8_367_22_Rect(bs, offset + 3);
    case AIS8_366_22_SHAPE_SECTOR:
      return new Ais8_367_22_Sector(bs, offset + 3);
    case AIS8_366_22_SHAPE_POLYLINE:  // FALLTHROUGH
    case AIS8_366_22_SHAPE_POLYGON:
      return new Ais8_367_22_Poly(bs, offset + 3, area_shape);
    case AIS8_366_22_SHAPE_TEXT:
      return new Ais8_367_22_Text(bs, offset + 3);
    case AIS8_366_22_SHAPE_RESERVED_6:  // FALLTHROUGH
    case AIS8_366_22_SHAPE_RESERVED_7:  // FALLTHROUGH
      // Leave area as 0 to indicate error.
      break;
    case AIS8_366_22_SHAPE_ERROR:
      break;
    default:
      assert(false);
  }
  return nullptr;
}

Ais8_367_22::Ais8_367_22(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), version(0), link_id(0), notice_type(0),
      month(0), day(0), hour(0), minute(0), duration_minutes(0), spare2(0) {
  assert(dac == 367);
  assert(fi == 22);

  if (num_bits < 216 || num_bits > 1016) {
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
  version = bs.ToUnsignedInt(56, 6);
  link_id = bs.ToUnsignedInt(62, 10);
  notice_type = bs.ToUnsignedInt(72, 7);
  month = bs.ToUnsignedInt(79, 4);
  day = bs.ToUnsignedInt(83, 5);
  hour = bs.ToUnsignedInt(88, 5);
  minute = bs.ToUnsignedInt(93, 6);
  duration_minutes = bs.ToUnsignedInt(99, 18);
  spare2 = bs.ToUnsignedInt(117, 3);

  const int num_sub_areas = static_cast<int>(
      floor((num_bits - 120) / static_cast<float>(SUB_AREA_BITS)));

  for (int area_idx = 0; area_idx < num_sub_areas; area_idx++) {
    const size_t start = 120 + area_idx * SUB_AREA_BITS;
    Ais8_367_22_SubArea *area = ais8_367_22_subarea_factory(bs, start);
    if (area != nullptr) {
      sub_areas.push_back(area);
    } else {
      status = AIS_ERR_BAD_SUB_SUB_MSG;
      return;
    }
  }

  // TODO(schwehr): Save the spare bits at the end of the message.
  assert(bs.GetRemaining() < 6);
  status = AIS_OK;
}

Ais8_367_22::~Ais8_367_22() {
  // Switch to unique_ptr.
  for (size_t i = 0; i < sub_areas.size(); i++) {
    delete sub_areas[i];
    sub_areas[i] = nullptr;
  }
}

}  // namespace libais
