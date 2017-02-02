// Defined sort of by Nav55, but being tweaked by the USCG and
// will hopefully become a RTCM standard.  Also hoping that this
// will be harmonized with the IMO Circ 289.

#include <cmath>

#include "ais.h"

namespace libais {

const char *shape_names[8] = {"Circle/Pt", "Rect", "Sector", "Polyline",
                              "Polygon", "Text", "Reserved_6", "Reserved_7"};

const char *ais8_366_22_notice_names[AIS8_366_22_NUM_NAMES] = {
  "Caution Area: Marine mammals habitat (no whales obs)",  // 0
  "Caution Area: Marine mammals in area – reduce speed",  // 1
  "Caution Area: Marine mammals in area – stay clear",  // 2
  "Caution Area: Marine mammals in area – report sightings",  // 3
  "Caution Area: Protected habitat – reduce speed",  // 4
  "Caution Area: Protected habitat – stay clear",  // 5
  "Caution Area: Protected habitat – no fishing or anchoring",  // 6
  "Caution Area: Derelicts (drifting objects)",  // 7
  "Caution Area: Traffic congestion",  // 8
  "Caution Area: Marine event",  // 9
  "Caution Area: Divers down",  // 10
  "Caution Area: Swim area",  // 11
  "Caution Area: Dredge operations",  // 12
  "Caution Area: Survey operations",  // 13
  "Caution Area: Underwater operation",  // 14
  "Caution Area: Seaplane operations",  // 15
  "Caution Area: Fishery – nets in water",  // 16
  "Caution Area: Cluster of fishing vessels",  // 17
  "Caution Area: Fairway closed",  // 18
  "Caution Area: Harbour closed",  // 19
  "Caution Area: Risk (define in Associated text field)",  // 20
  "Caution Area: Underwater vehicle operation",  // 21
  "(reserved for future use)",  // 22
  "Environmental Caution Area: Storm front (line squall)",  // 23
  "Environmental Caution Area: Hazardous sea ice",  // 24
  "Environmental Caution Area: Storm warning (cell or line of storms)",  // 25
  "Environmental Caution Area: High wind",  // 26
  "Environmental Caution Area: High waves",  // 27
  "Environmental Caution Area: Restricted visibility (fog, rain, etc.)",  // 28
  "Environmental Caution Area: Strong currents",  // 29
  "Environmental Caution Area: Heavy icing",  // 30
  "(reserved for future use)",  // 31
  "Restricted Area: Fishing prohibited",  // 32
  "Restricted Area: No anchoring.",  // 33
  "Restricted Area: Entry approval required prior to transit",  // 34
  "Restricted Area: Entry prohibited",  // 35
  "Restricted Area: Active military OPAREA",  // 36
  "Restricted Area: Firing – danger area.",  // 37
  "Restricted Area: Drifting Mines",  // 38
  "(reserved for future use)",  // 39
  "Anchorage Area: Anchorage open",  // 40
  "Anchorage Area: Anchorage closed",  // 41
  "Anchorage Area: Anchoring prohibited",  // 42
  "Anchorage Area: Deep draft anchorage",  // 43
  "Anchorage Area: Shallow draft anchorage",  // 44
  "Anchorage Area: Vessel transfer operations",  // 45
  "(reserved for future use)",  // 46
  "(reserved for future use)",  // 47
  "(reserved for future use)",  // 48
  "(reserved for future use)",  // 49
  "(reserved for future use)",  // 50
  "(reserved for future use)",  // 51
  "(reserved for future use)",  // 52
  "(reserved for future use)",  // 53
  "(reserved for future use)",  // 54
  "(reserved for future use)",  // 55
  "Security Alert – Level 1",  // 56
  "Security Alert – Level 2",  // 57
  "Security Alert – Level 3",  // 58
  "(reserved for future use)",  // 59
  "(reserved for future use)",  // 60
  "(reserved for future use)",  // 61
  "(reserved for future use)",  // 62
  "(reserved for future use)",  // 63
  "Distress Area: Vessel disabled and adrift",  // 64
  "Distress Area: Vessel sinking",  // 65
  "Distress Area: Vessel abandoning ship",  // 66
  "Distress Area: Vessel requests medical assistance",  // 67
  "Distress Area: Vessel flooding",  // 68
  "Distress Area: Vessel fire/explosion",  // 69
  "Distress Area: Vessel grounding",  // 70
  "Distress Area: Vessel collision",  // 71
  "Distress Area: Vessel listing/capsizing",  // 72
  "Distress Area: Vessel under assault",  // 73
  "Distress Area: Person overboard",  // 74
  "Distress Area: SAR area",  // 75
  "Distress Area: Pollution response area",  // 76
  "(reserved for future use)",  // 77
  "(reserved for future use)",  // 78
  "(reserved for future use)",  // 79
  "Instruction: Contact VTS at this point/juncture",  // 80
  "Instruction: Contact Port Administration at this point/juncture",  // 81
  "Instruction: Do not proceed beyond this point/juncture",  // 82
  "Instruction: Await instructions prior to proceeding beyond "
  "this point/juncture",  // 83
  "Proceed to this location – await instructions",  // 84
  "Clearance granted – proceed to berth",  // 85
  "(reserved for future use)",  // 86
  "(reserved for future use)",  // 87
  "Information: Pilot boarding position",  // 88
  "Information: Icebreaker waiting area",  // 89
  "Information: Places of refuge",  // 90
  "Information: Position of icebreakers",  // 91
  "Information: Location of response units",  // 92
  "VTS active target",  // 93
  "Rogue or suspicious vessel",  // 94
  "Vessel requesting non-distress assistance",  // 95
  "Chart Feature: Sunken vessel",  // 96
  "Chart Feature: Submerged object",  // 97
  "Chart Feature: Semi-submerged object",  // 98
  "Chart Feature: Shoal area",  // 99
  "Chart Feature: Shoal area due north",  // 100
  "Chart Feature: Shoal area due east",  // 101
  "Chart Feature: Shoal area due south",  // 102
  "Chart Feature: Shoal area due west",  // 103
  "Chart Feature: Channel obstruction",  // 104
  "Chart Feature: Reduced vertical clearance",  // 105
  "Chart Feature: Bridge closed",  // 106
  "Chart Feature: Bridge partially open",  // 107
  "Chart Feature: Bridge fully open",  // 108
  "(reserved for future use)",  // 109
  "(reserved for future use)",  // 110
  "(reserved for future use)",  // 111
  "Report from ship: Icing info",  // 112
  "(reserved for future use)",  // 113
  "Report from ship: Misc information - see associated text field",  // 114
  "(reserved for future use)",  // 115
  "(reserved for future use)",  // 116
  "(reserved for future use)",  // 117
  "(reserved for future use)",  // 118
  "(reserved for future use)",  // 119
  "Route: Recommended route",  // 120
  "Route: Alternative route",  // 121
  "Route: Recommended route through ice",  // 122
  "(reserved for future use)",  // 123
  "(reserved for future use)",  // 124
  "Other – Define in associated text field",  // 125
  "Cancellation – cancel area as identified by Message Linkage ID",  // 126
  "Undefined (default)"  // 127
};


Ais8_366_22::Ais8_366_22(const char *nmea_payload, const size_t pad)
    : Ais8(nmea_payload, pad), link_id(0), notice_type(0), month(0), day(0),
      utc_hour(0), utc_minute(0), duration_minutes(0) {
  assert(dac == 366);
  assert(fi == 22);

  if (num_bits <= 208 || num_bits >= 1020) {
    status = AIS_ERR_BAD_BIT_COUNT;
    return;
  }

  AisBitset bs;
  const AIS_STATUS r = bs.ParseNmeaPayload(nmea_payload, pad);
  if (r != AIS_OK) {
    status = r;
    return;
  }

  link_id = bs.ToUnsignedInt(56, 10);
  notice_type = bs.ToUnsignedInt(66, 7);
  month = bs.ToUnsignedInt(73, 4);
  day = bs.ToUnsignedInt(77, 5);
  utc_hour = bs.ToUnsignedInt(82, 5);
  utc_minute = bs.ToUnsignedInt(87, 6);

  duration_minutes = bs.ToUnsignedInt(93, 18);

  const int num_sub_areas = static_cast<int>(floor((num_bits - 111)/90.));
  for (int area_idx = 0; area_idx < num_sub_areas; area_idx++) {
    Ais8_366_22_SubArea *area =
        ais8_366_22_subarea_factory(bs, 111 + 90*area_idx);
    if (area) {
      sub_areas.push_back(area);
    } else {
      status = AIS_ERR_BAD_SUB_SUB_MSG;
      return;
    }
  }

  assert(bs.GetRemaining() == 0);
  status = AIS_OK;
}

Ais8_366_22::~Ais8_366_22() {
  // Switch to unique_ptr.
  for (size_t i = 0; i < sub_areas.size(); i++) {
    delete sub_areas[i];
#ifndef NDEBUG
    sub_areas[i] = NULL;
#endif
  }
}

// Lookup table for the scale factors to decode the length / distance fields.
// The index is the "Scale Factor"
static int scale_multipliers[4] = {1, 10, 100, 1000};

Ais8_366_22_Circle::Ais8_366_22_Circle(const AisBitset &bs,
                                       const size_t offset) {
  const int scale_factor = bs.ToUnsignedInt(offset + 3, 2);
  position = bs.ToAisPoint(offset + 5, 55);
  // TODO(schwehr): precision? And bit counts for radius  and spare?
  // TODO(schwehr): collapse these numbers
  radius_m =
      bs.ToUnsignedInt(offset + 60, 12) * scale_multipliers[scale_factor];
  spare = bs.ToUnsignedInt(offset + 72, 16);
}

Ais8_366_22_Rect::Ais8_366_22_Rect(const AisBitset &bs,
                                   const size_t offset) {
  const int scale_factor = bs.ToUnsignedInt(offset + 3, 2);
  position = bs.ToAisPoint(offset + 5, 55);
  e_dim_m = bs.ToUnsignedInt(offset + 60, 8) * scale_multipliers[scale_factor];
  n_dim_m = bs.ToUnsignedInt(offset + 68, 8) * scale_multipliers[scale_factor];
  orient_deg = bs.ToUnsignedInt(offset + 76, 9);
  spare = bs.ToUnsignedInt(offset + 85, 5);
}

Ais8_366_22_Sector::Ais8_366_22_Sector(const AisBitset &bs,
                                       const size_t offset) {
  const int scale_factor = bs.ToUnsignedInt(offset + 3, 2);
  position = bs.ToAisPoint(offset + 5, 55);
  radius_m =
      bs.ToUnsignedInt(offset + 60, 12) * scale_multipliers[scale_factor];
  left_bound_deg = bs.ToUnsignedInt(offset + 72, 9);
  right_bound_deg = bs.ToUnsignedInt(offset + 81, 9);
}

Ais8_366_22_Polyline::Ais8_366_22_Polyline(const AisBitset &bs,
                                           const size_t offset) {
  const int scale_factor = bs.ToUnsignedInt(offset + 3, 2);
  for (size_t i = 0; i < 4; i++) {
    const int angle = bs.ToUnsignedInt(offset + 5 + (i*21), 10);
    const int dist = bs.ToUnsignedInt(offset + 15 + (i*21), 11) *
                     scale_multipliers[scale_factor];
    if (dist == 0) {
      break;
    }
    angles.push_back(angle);
    dists_m.push_back(dist);
  }
  spare = bs[offset + 89];
}

// TODO(schwehr): Merge with Polyline.
Ais8_366_22_Polygon::Ais8_366_22_Polygon(const AisBitset &bs,
                                         const size_t offset) {
  const int scale_factor = bs.ToUnsignedInt(offset + 3, 2);
  for (size_t i = 0; i < 4; i++) {
    const int angle = bs.ToUnsignedInt(offset + 5 + (i*21), 10);
    const int dist = bs.ToUnsignedInt(offset + 15 + (i*21), 11) *
                     scale_multipliers[scale_factor];
    if (dist == 0) {
      break;
    }
    angles.push_back(angle);
    dists_m.push_back(dist);
  }
  spare = bs[offset + 89];
}

Ais8_366_22_Text::Ais8_366_22_Text(const AisBitset &bs,
                                   const size_t offset) {
  text = string(bs.ToString(offset + 3, 84));
  spare = bs.ToUnsignedInt(offset + 87, 3);
}

// Call the appropriate constructor
Ais8_366_22_SubArea *
ais8_366_22_subarea_factory(const AisBitset &bs,
                            const size_t offset) {
  const Ais8_366_22_AreaShapeEnum area_shape =
      (Ais8_366_22_AreaShapeEnum)bs.ToUnsignedInt(offset, 3);

  switch (area_shape) {
  case AIS8_366_22_SHAPE_CIRCLE:
    return new Ais8_366_22_Circle(bs, offset);
  case AIS8_366_22_SHAPE_RECT:
    return new Ais8_366_22_Rect(bs, offset);
  case AIS8_366_22_SHAPE_SECTOR:
    return new Ais8_366_22_Sector(bs, offset);
  case AIS8_366_22_SHAPE_POLYLINE:
    return new Ais8_366_22_Polyline(bs, offset);
  case AIS8_366_22_SHAPE_POLYGON:
    return new Ais8_366_22_Polygon(bs, offset);
  case AIS8_366_22_SHAPE_TEXT:
    return new Ais8_366_22_Text(bs, offset);
  case AIS8_366_22_SHAPE_RESERVED_6:  // FALLTHROUGH
  case AIS8_366_22_SHAPE_RESERVED_7:  // FALLTHROUGH
    // Leave area as 0 to indicate error
    break;
  case AIS8_366_22_SHAPE_ERROR:
    break;
  default:
    assert(false);
  }
  return nullptr;
}

}  // namespace libais
